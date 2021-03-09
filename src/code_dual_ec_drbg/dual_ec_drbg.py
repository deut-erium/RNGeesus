from tqdm import tqdm
from fastecdsa.curve import P256
from fastecdsa.point import Point
from gmpy2 import *
import random, os
from time import time


def xgcd(a, b):
    a1 = 1
    b1 = 0
    a2 = 0
    b2 = 1
    aneg = 1
    bneg = 1
    if(a < 0):
        a = -a
        aneg = -1
    if(b < 0):
        b = -b
        bneg = -1
    while (1):
        quot = -(a // b)
        a = a % b
        a1 = a1 + quot*a2
        b1 = b1 + quot*b2
        if(a == 0):
            return (b, a2*aneg, b2*bneg)
        quot = -(b // a)
        b = b % a
        a2 = a2 + quot*a1
        b2 = b2 + quot*b1
        if(b == 0):
            return (a, a1*aneg, b1*bneg)

def SqrRoots(a, n):
    def inverse_mod(a, n):
        (g, xa, xb) = xgcd(a, n)
        if(g != 1):
            raise ValueError(
                "***** Error *****: {0} has no inverse (mod {1}) as their gcd is {2}, not 1.".format(a, n, g))
        return xa % n

    def TSRsqrtmod(a, grpord, p):
        ordpow2 = 0
        non2 = grpord
        while(not ((non2 & 0x01) == 1)):
            ordpow2 += 1
            non2 //= 2
        for g in range(2, grpord-1):
            if (pow(g, grpord//2, p) != 1):
                break
        g = pow(g, non2, p)
        gpow = 0
        atweak = a
        for pow2 in range(0, ordpow2+1):
            if(pow(atweak, non2*2**(ordpow2-pow2), p) != 1):
                gpow += 2**(pow2-1)
                atweak = (atweak * pow(g, 2**(pow2-1), p)) % p
        d = inverse_mod(2, non2)
        tmp = pow(a*pow(g, gpow, p), d, p)
        return (tmp*inverse_mod(pow(g, gpow//2, p), p)) % p
    x1 = TSRsqrtmod(a, n-1, n)
    return x1, -x1 % n

def urandbits(n):
    return int.from_bytes(os.urandom(n//8),'big')

class Dual_EC:
    def __init__(self, seed, P, Q):
        self.state = seed
        self.P = P
        self.Q = Q
        # P and Q are not randomly chosen, but they're some numbers given to us!

    def next(self):
        r = (self.state*self.P).x
        self.state = (r*self.P).x
        t = (r * self.Q).x
        self.t = t

        return (t & ((1<<240)-1)) # get least signif 240 bits        


class Breaker(Dual_EC):
    def solution_exists(self, x):
        '''
            Checks if a solution exists for a given x-coordinate. Also outputs solutions in case they exist
            Returns (True, solutions) or (False, ())
        '''
        rhs = P256.evaluate(x)
        l = legendre(rhs, P256.p)
        if (l == -1):
            return (False, ())
        elif (l == 0):
            return (True, (0))
        else:
            p = tuple(set(SqrRoots(rhs, P256.p)))
            return (True, p)

    def get_random_point(self):
        '''
            Obtain a random point on the curve, given a generator
        '''
        x = random.randint(1, P256.q - 1)
        # Obtain point using generator
        return (x*P256.G)

    def __init__(self, seed):
        P = self.get_random_point()
        self.e = random.randint(1, P256.q - 1)
        Q = (P * self.e)
        self.d = int(invert(self.e,P256.q))
        # we ensure that P and Q are related, that allows us to exploit this possible backdoor
        # Q = e*P
        
        super().__init__(seed, P, Q)
        
    def possible_points(self, output):
        '''
            Given the output 240 bits, we want to obtain the possible points r*Q could be.
        '''
        l = []
        for i in tqdm(range(2 ** 16)):
            poss_x = (i << 240) + output
            a, b = self.solution_exists(poss_x)
            if a:
                for j in b:
                    p = Point(poss_x, j, curve=P256)
                    assert P256.is_point_on_curve((p.x,p.y)), "Point not on curve? How!"
                    l.append(p)
        return l
        
    def break_dec(self):
        '''
            Try to recover the current state of the Dual_EC_DRBG - can't recover older states!
        '''

        it = 240
        oup = self.next()
        l = self.possible_points(oup)

        # find all possible next states
        next_s = list(set([(p * self.d).x for p in l]))

        next_pts = [Dual_EC(i,self.P,self.Q).next() for i in next_s]
        # need first 240 bits to generate an initial list of possible states

        # the paper claims that no more than 256 bits are needed to break the RNG
        oup1 = self.next()
        p = 31

        print(f"Initial count of possible states : {len(next_s)}".format())

        inds = list(range(len(next_s)))
        
        while (p > 0) :

            if (len(inds) <= 1):
                break

            inds = list(filter(lambda x: ((oup1 & (1<<p)) == (next_pts[x] & (1<<p))),inds))
            
            it += 1
            p -= 1

        assert len(inds) == 1, [next_s[x] for x in inds]
        print(f"Old state - {next_s[inds[0]]}".format())
        print(f"New state - {((next_s[inds[0]] * self.P).x * self.P).x}".format())
        assert ((next_s[inds[0]] * self.P).x * self.P).x == self.state
        print(f"Took {it} bits to recover the seed".format())
        return next_s[inds[0]], it
                
if __name__ == '__main__':

    rand_seed = urandbits(256)
    d = Breaker(rand_seed)
    start_t = time()
    m = d.break_dec()
    print(m)
    print("Time taken: ", time()-start_t)