""" Some functions to implement Linear Feedback Shift Register. Berlekamp-Massey algo and Geffe Generator """

from functools import reduce
from z3 import *
import tqdm, itertools

def all_smt(s, var_list):
    """
    yielding all satisfying models over `var_list` on a 
    z3.Solver() instance `s` containing constraints
    """
    def block_term(s, m, t):
        s.add(t != m.eval(t))

    def fix_term(s, m, t):
        s.add(t == m.eval(t))

    def all_smt_rec(terms):
        if sat == s.check():
            m = s.model()
            yield m
            for i in range(len(terms)):
                s.push()
                block_term(s, m, terms[i])
                for j in range(i):
                    fix_term(s, m, terms[j])
                for m in all_smt_rec(terms[i:]):
                    yield m
                s.pop()
    for m in all_smt_rec(var_list):
        yield m

class LFSR:
    """ Normal LFSR impl with pythonic inputs. Everything is in `GF(2)`
    n-bit LFSR defined by given feedback polynomial
    seed = MSB to LSB list of bits
    feedback_polynomial = MSB to LSB list of bits
    """

    def __init__(self, seed, poly):
        assert len(seed) == len(poly), "Error: Seed and taps poly  should be of same length"
        self._seed = seed.copy()        # Sn, Sn-1, ..., S0
        self._comb_poly = poly          # C0, C1, ..., Cn
    
    def next_bit(self):
        """ Generate next output bit """
        if type(self._comb_poly[0]) == BitVecRef:
            tapped = self._comb_poly.copy()
        else:
            tapped = [self._seed[i] for i,j in enumerate(self._comb_poly) if j == 1]
        xored = reduce(lambda x,y: x^y, tapped)
        opt = self._seed.pop(0)
        self._seed.append(xored)
        return opt
    
    def get_lfsr(self, steps):
        """ Get next `steps` number of output bits """
        opt = [self.next_bit() for _ in range(steps)]
        return opt
    
    def set_seed(self, new_seed):
        """ Set the new seed to generate new LFSR using same polynomial """
        self._seed = new_seed.copy()

class Berlekamp_Massey:
    """ Berlekamp - Massey algo: PYTHON IMPLEMENTATION
    i/p:    S:  `list` of 0s and 1s, Sn, Sn-1, Sn-2, ... S1, S0.
    o/p:   min degree of C, Feedback Polynomial, anything else that we want 
    """

    def __init__(self, S):
        self.S = S
        C = [1]     # Connection polynomial. The one that generates next o/p bit. 1, C1, C2, ..., CL.
        L = 0       # Minimal size of LFSR at o/p
        m = -1      # num iterations since L and B were updated.
        B = [1]     # Previous value of C, before it was updated.
        n = 0       # counter. i.e. the iterator
        N = len(S)  # The length of i/p

        while(n < N):
            bit_calc = [i&j for i,j in zip(C,S[n-L:n+1][::-1])]
            d = reduce(lambda x, y: x^y, bit_calc,0)
            if d:
                c_temp = C.copy()
                lc = len(C)
                next_C = [0]*(n-m) + B + [0]*(lc - len(B) - n + m)
                C = [i^j for i,j in zip(C,next_C)] + next_C[lc:]
                # print(n, d, m,L, B, next_C, C)

                if L <= n>>1:
                    L = n + 1 - L
                    m = n
                    B = c_temp.copy()
            n += 1
        self._L = L
        self._C = C[::-1]
        self._seed = S[:L]
        assert len(self._seed) + 1 == len(self._C)

    def get_seed(self):
        return self._seed

    def get_taps(self):
        return self._C[1:]

    def get_degree(self):
        return self._L

class UnLFSR_Z3:
    """ Similar to berlekamp in the sense that it finds the seed and the comb poly using z3 solver. """

    def __init__(self, opt, degree_guess):
        """ opt is list of 0s and 1s. 1st bit 1st """
        self._opt = opt.copy()
        self._seed = [BitVec(f'k_{i}',1) for i in range(degree_guess)]
        self._poly = [BitVec(f'c_{i}',1) for i in range(degree_guess)]

    def solve(self):
        s = Solver()
        lfsr = LFSR(self._seed, self._poly)
        for i in range(len(self._opt)):
            s.add(lfsr.next_bit() == self._opt[i])
        if s.check() == sat:
            model = s.model()
            print(len(model))
            sd = ''.join(str(model[i]) for i in self._seed)
            taps = ''.join(str(model[i]) for i in self._poly)
            return sd,taps
        else:
            print("ERROR: unsolvable... unpossible!!")

class Geffe:
    """ Geffe Generator with solver in z3 as brute force as well. We need to know  the combination polynomial beforehand """

    def __init__(self, l1, l2, l3, c1, c2, c3):
        self._l1, self._l2, self._l3 = l1, l2, l3
        self._c1, self._c2, self._c3 = c1, c2, c3
        self._lfsrs = [LFSR(self._l1,c1), LFSR(self._l2, c2), LFSR(self._l3, c3)]
    
    def next_bit(self):
        bits = [lfsr.next_bit() for lfsr in self._lfsrs]
        # Equiv to if bits[0] then bits[1] else bits[2] in GF(2)
        return (bits[0] & bits[1]) | ((~bits[0]) & bits[2])
    
    def get_seqn(self, steps):
        """ Return steps geffe generated bits """
        return [self.next_bit() for _ in range(steps)]

    def solve(self, opt):
        """ opt is a list of output bits gen by jeff gen """
        s = Solver()
        for b in opt:
            s.add(self.next_bit() == b)
        for model in all_smt(s, self._l1 + self._l2 + self._l3):
            one = ''.join(str(model[k]) for k in self._l1)
            two = ''.join(str(model[k]) for k in self._l2)
            thr = ''.join(str(model[k]) for k in self._l3)
            yield (one,two,thr)
    
    def __temp_geffe(self, lfsr0: LFSR, lfsr1: LFSR, lfsr2: LFSR, steps: int):
        ans = []
        for _ in range(steps):
            bits = [lfsr0.next_bit(), lfsr1.next_bit(), lfsr2.next_bit()]
            ans.append((bits[0] & bits[1]) | ((~bits[0]) & bits[2]))
        return ans
    
    def solve_bruteforce(self, opt: list):
        n = len(opt)
        lfsr0 = self._lfsrs[0]
        lfsr1 = self._lfsrs[1]
        lfsr2 = self._lfsrs[2]

        # >75% match of opt with x3 i.e. lfsr2
        possible_seeds2 = []
        m2 = 0
        for seed2 in tqdm.tqdm(itertools.product('01', repeat=len(lfsr2._comb_poly)), total=pow(2, len(lfsr2._comb_poly))):
            lfsr2.set_seed(list(map(int,seed2)))
            x3 = lfsr2.get_lfsr(len(opt))
            corr = sum(x==y for x,y in zip(opt, x3))
            if corr >= int(0.70*n):
                possible_seeds2.append(''.join(seed2))
                break
        assert len(possible_seeds2) >=1, "Error: No x3 found, less data supplied."

        # > 75% match of opt with x2 i.e. lfsr1
        possible_seeds1 = []
        m1 = 0
        for seed1 in tqdm.tqdm(itertools.product('01', repeat=len(lfsr1._comb_poly)), total=pow(2, len(lfsr1._comb_poly))):
            lfsr1.set_seed(list(map(int,seed1)))
            x2 = lfsr1.get_lfsr(len(opt))
            corr = sum(x==y for x,y in zip(opt, x2))
            if corr >= int(0.70*n):
                possible_seeds1.append(''.join(seed1))
                break
        assert len(possible_seeds1) >=1, "Error: No x2 found, less data supplied"

        candidates = [(x, y) for x in possible_seeds1 for y in possible_seeds2]
        # print(candidates)

        # Now iterate through the remainig candidates and the all possible lfsr0 seeds
        # res = []
        # for s1, s2 in candidates:
        #     lfsr1.set_seed(list(map(int, s1)))
        #     lfsr2.set_seed(list(map(int, s2)))
        #     for s0 in tqdm.tqdm(itertools.product('01', repeat=len(lfsr0._comb_poly))):
        #         lfsr0.set_seed(list(map(int, s0)))
        #         geffe_opt = self.__temp_geffe(lfsr0, lfsr1, lfsr2, n)
        #         if geffe_opt == opt:
        #             res.append((s0, s1, ''.join(s2)))
        # print("Number of results: ", res)
        # return res[0]
        return (possible_seeds1[0], possible_seeds2[0])
