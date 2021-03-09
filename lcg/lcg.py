from z3 import *
from fpylll import IntegerMatrix, LLL
from time import time
from sympy import QQ, Matrix, invert
from functools import reduce
from sympy.polys.matrices import DomainMatrix
import sys, os, gmpy2
sys.setrecursionlimit(100000)

set_param('parallel.enable', True)
set_param('parallel.threads.max', 32)
set_param('sat.local_search_threads', 4)
set_param('sat.threads', 4)

def urandbits(n):
    return int.from_bytes(os.urandom(n//8),'big')

def gcd(*args):
    return reduce(gmpy2.gcd,args)

class lcg:
    def __init__(self, seed, a, b, m):
        self.state = seed
        self.a = a
        self.b = b
        self.m = m

    def next(self):
        self.state = (self.state * self.a + self.b) % self.m
        return self.state

class Breaker_lcg(lcg):
    def __init__(self, seed, a, b, m):
        super().__init__(seed,a,b,m)

    def break_lcg(self, ntimes):

        outputs = [self.next() for i in range(ntimes)]

        diffs = [(j-i) for i,j in zip(outputs,outputs[1:])]

        prods = [(b**2-a*c) for a,b,c in zip(diffs,diffs[1:],diffs[2:])]

        p = gcd(*prods)

        a = (diffs[1]*gmpy2.invert(diffs[0],p))%p

        b = (outputs[1]-a*outputs[0])%p

        assert all(j == (a*i + b)%p for i,j in zip(outputs,outputs[1:]))
        print(f"Recovered internal constants from {len(outputs)} outputs : p = {p} a = {a} b = {b}")
        return p,a,b

def all_smt(s, initial_terms):
    """
    Yield all satisfying models for solver s 
    for the given set of `initial_terms`
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
    for m in all_smt_rec(list(initial_terms)):
        yield m

class truncated_lcg:
    def __init__(self, seed, a, b, n,truncation):
        self.state = seed
        self.a = a
        self.b = b
        self.n = n
        self.truncation = truncation

    def next(self):
        self.state = ((self.a * self.state) + self.b) % self.n
        # print(self.state)
        return (self.state >> self.truncation)
        
        
class Breaker(truncated_lcg):
    def __init__(self, seed, a, b, n, truncation, **kwargs):
        super().__init__(seed, a, b, n, truncation)
        self.n_bitlen = n.bit_length()
        self.known_a: bool = kwargs.get('known_a', True)
        self.known_b: bool = kwargs.get('known_b', True)
        self.known_n: bool = kwargs.get('known_n', True)
        
    def break_sat(self, outputs):
        """
        Thought this wont suck
        well this sucks too XD
        """
        seed0 = BitVec('seed0', self.n_bitlen)
        seed = ZeroExt(self.n_bitlen,seed0)
        s = Solver()

        if (self.known_a):
            a = BitVecVal(self.a, self.n_bitlen)
        else:
            a = BitVec('a', self.n_bitlen)
            
        if (self.known_b):
            b = BitVecVal(self.b, self.n_bitlen)
        else:
            b = BitVec('b', self.n_bitlen)

        if (self.known_n):
            n = BitVecVal(self.n, self.n_bitlen)
        else:
            n = BitVec('n', self.n_bitlen)

        s.add(ULT(seed0,n),ULT(a,n),ULT(b,n),UGE(seed0,0),UGE(a,0),UGE(b,0))
        for v in outputs:
            seed = simplify(URem(ZeroExt(self.n_bitlen,a)*seed+ZeroExt(self.n_bitlen,b), ZeroExt(self.n_bitlen,n)))
            s.add(v == LShR(seed,self.truncation))

        start_time, last_time = time(), time()
        terms = [seed0,a,b,n]

        guess = []

        for m in all_smt(s,terms):
            SAT_guessed_seed = m[seed0]
            A = m.eval(a)
            B = m.eval(b)
            N = m.eval(n)
            print(f"{SAT_guessed_seed = } {A = } {B = } {N = }")
            guess.append((SAT_guessed_seed,A,B,N))
        print("Total time taken(SAT) :",time()-start_time)
        return guess

    def shorten(self,u):
        for i in range(u.nrows):
            u[i,0] %= self.n
            if 2*u[i,0] >=self.n:
                u[i,0]-=self.n

    def break_lattice(self, outputs):
        k = len(outputs)
        start_time = time()
        L = IntegerMatrix(k, k)
        v = IntegerMatrix(k, 1)
        U = IntegerMatrix.identity(k)
        for i in range(k):
            L[i, 0] = self.a**i
            L[i, i] = -1
            v[i, 0] = -(outputs[i] << self.truncation) % self.n
        L[0,0] = self.n
            
        v = L * v
        
        for i in range(k):
            v[i, 0] += ((1 - self.a ** i) // (self.a - 1)) * self.b
            v[i, 0] %= self.n
            
        _ = LLL.reduction(L, U)

        u = (U * v)
        self.shorten(u)

        A = DomainMatrix.from_Matrix(Matrix(k, k, lambda i, j: L[i, j])).convert_to(QQ)
        b = DomainMatrix.from_Matrix(Matrix(k, 1, lambda i, j: u[i, 0])).convert_to(QQ)
        M = (A.inv() * b).to_Matrix()

        next_st = (outputs[0] << self.truncation) | int(M[0, 0] % self.n)
        try:
            lattice_guessed_seed = (invert(self.a,self.n)*(next_st-self.b))%self.n
            print(f"{lattice_guessed_seed = }")
            print(f"Total time taken(LLL) : {time()-start_time}")
            return [lattice_guessed_seed]
        except:
            # if inverse of a in p does not exist
            seed = BitVec('seed', self.n_bitlen)
            
            s = Solver()
            s.add(ULT(seed,self.n),next_st == simplify(URem(self.a * ZeroExt(self.n_bitlen, seed) + self.b, self.n)))
            
            guess = []

            for m in all_smt(s, [seed]):
                lattice_guessed_seed = m[seed]
                print(f"{lattice_guessed_seed = }")
                guess.append(lattice_guessed_seed)
            print(f"Total time taken(LLL) : {time()-start_time}")
            return guess