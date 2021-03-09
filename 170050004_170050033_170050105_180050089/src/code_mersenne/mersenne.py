import random
from time import time
from collections import Counter
from statistics import mode
from functools import reduce
from z3 import *


def seed_arr_len(arr):
    """
    Figuring out the length of seed in words which was input to
    init_by_array
    """
    mode_vals = []
    for j in range(2,10):
        x = [i for i in range(624) if abs(arr[i]-arr[j])<624]
        mode_vals.append(
            mode([j-i for i,j in zip(x,x[1:])])
        )
    return mode(mode_vals)

def all_smt(s, initial_terms):
    """
    yielding all satisfying models over `initial_terms` on a 
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
    for m in all_smt_rec(list(initial_terms)):
        yield m


class MT19937:
    """
    Standard MT19937 instance for both 32 bit and 64 bit variants
    """
    def __init__(self, c_seed=0, bit_64=False):
        """
        initialize the mersenne twister with `c_seed`
        `bit_64` if True, would initialize the 64 variant of MT19937
        `c_seed` is 64-bit if `bit_64` set to True
        """
        # MT19937
        if bit_64:
            (self.w, self.n, self.m, self.r) = (64, 312, 156, 31)
            self.a = 0xB5026F5AA96619E9
            (self.u, self.d) = (29, 0x5555555555555555)
            (self.s, self.b) = (17, 0x71D67FFFEDA60000)
            (self.t, self.c) = (37, 0xFFF7EEE000000000)
            self.l = 43
            self.f = 6364136223846793005
        else:
            (self.w, self.n, self.m, self.r) = (32, 624, 397, 31)
            self.a = 0x9908B0DF
            (self.u, self.d) = (11, 0xFFFFFFFF)
            (self.s, self.b) = (7, 0x9D2C5680)
            (self.t, self.c) = (15, 0xEFC60000)
            self.l = 18
            self.f = 1812433253
        self.MT = [0 for i in range(self.n)]
        self.index = self.n + 1
        self.lower_mask = (1 << self.r) - 1  # 0x7FFFFFFF
        self.upper_mask = (1 << self.r)  # 0x80000000
        self.seed_mt(c_seed)

    def seed_mt(self, num):
        """initialize the generator from a seed"""
        self.MT[0] = num
        self.index = self.n
        for i in range(1, self.n):
            temp = self.f * (self.MT[i - 1] ^
                             (self.MT[i - 1] >> (self.w - 2))) + i
            self.MT[i] = temp & ((1 << self.w) - 1)

    def twist(self):
        """ Generate the next n values from the series x_i"""
        for i in range(0, self.n):
            x = (self.MT[i] & self.upper_mask) + \
                (self.MT[(i + 1) % self.n] & self.lower_mask)
            xA = x >> 1
            if (x % 2) != 0:
                xA = xA ^ self.a
            self.MT[i] = self.MT[(i + self.m) % self.n] ^ xA
        self.index = 0

    def extract_number(self):
        """
        extract tampered state at internal index i
        if index reaches end of state array, twist and set it to 0
        """
        if self.index >= self.n:
            self.twist()
        y = self.MT[self.index]
        y = y ^ ((y >> self.u) & self.d)
        y = y ^ ((y << self.s) & self.b)
        y = y ^ ((y << self.t) & self.c)
        y = y ^ (y >> self.l)
        self.index += 1
        return y & ((1 << self.w) - 1)

    def get_state(self):
        """
        returning python compatible state
        """
        return (3, tuple(self.MT + [self.index]), None)


class MTpython(MT19937):
    """
    Additional functionality offered by MT of python3, namely 
    better (non linear) initialization 
    """
    def __init__(self, seed=0):
        MT19937.__init__(self, 0)
        self.seed(seed) #python seed initialization

    def init_by_array(self, init_key):
        """
        Initialization with an `init_key` array of 32-bit words for
        better randomization properties
        """
        self.seed_mt(19650218)
        i, j = 1, 0
        for k in range(max(self.n, len(init_key))):
            self.MT[i] = (self.MT[i] ^ (
                (self.MT[i - 1] ^ (self.MT[i - 1] >> 30)) * 1664525)) + init_key[j] + j
            self.MT[i] &= 0xffffffff
            i += 1
            j += 1
            if i >= self.n:
                self.MT[0] = self.MT[self.n - 1]
                i = 1
            if j >= len(init_key):
                j = 0
        for k in range(self.n - 1):
            self.MT[i] = (self.MT[i] ^ (
                (self.MT[i - 1] ^ (self.MT[i - 1] >> 30)) * 1566083941)) - i
            self.MT[i] &= 0xffffffff
            i += 1
            if i >= self.n:
                self.MT[0] = self.MT[self.n - 1]
                i = 1
        self.MT[0] = 0x80000000

    def init_32bit_seed(self, seed_32):
        """
        Just an oversimplification of `init_by_array` for single element array
        of upto 32 bit number
        """
        self.seed_mt(19650218)
        i = 1
        for k in range(self.n):
            self.MT[i] = (self.MT[i] ^ (
                (self.MT[i - 1] ^ (self.MT[i - 1] >> 30)) * 1664525)) + seed_32
            self.MT[i] &= 0xffffffff
            i += 1
            if i >= self.n:
                self.MT[0] = self.MT[self.n - 1]
                i = 1
        for k in range(self.n - 1):
            self.MT[i] = (self.MT[i] ^ (
                (self.MT[i - 1] ^ (self.MT[i - 1] >> 30)) * 1566083941)) - i
            self.MT[i] &= 0xffffffff
            i += 1
            if i >= self.n:
                self.MT[0] = self.MT[self.n - 1]
                i = 1
        self.MT[0] = 0x80000000

    def seed(self,seed_int):
        """
        Replication of random.seed of cpython when seed is an integer
        """
        self.init_by_array(self.int_to_array(seed_int))

    def random(self):
        """
        python random.random() call which yeilds a uniformly random
        floating point between [0,1] employing two MT 32 bits calls
        """
        a = self.extract_number()>>5
        b = self.extract_number()>>6
        return (a*67108864.0+b)*(1.0/9007199254740992.0)

    def int_to_array(self,k):
        """
        converting a big integer to equivalent list of 32-bit integers
        as would be passed into python seed process
        """
        if k==0:
            return [0]
        k_byte = int.to_bytes(k,(k.bit_length()+7)//8,'little')
        k_arr = [k_byte[i:i+4] for i in range(0,len(k_byte),4)]
        return [int.from_bytes(i,'little') for i in k_arr ]

    def array_to_int(self,arr):
        """
        converting list of 32-bit integers back to a big integer
        """
        arr_bytes  = b"".join([int.to_bytes(i,4,'little') for i in arr])
        return int.from_bytes( arr_bytes ,'little')


class Breaker():
    """
    Class for breaking and seed recovery of standard mersenne twister
    """
    def __init__(self, bit_64=False):
        if bit_64:
            (self.w, self.n, self.m, self.r) = (64, 312, 156, 31)
            self.a = 0xB5026F5AA96619E9
            (self.u, self.d) = (29, 0x5555555555555555)
            (self.s, self.b) = (17, 0x71D67FFFEDA60000)
            (self.t, self.c) = (37, 0xFFF7EEE000000000)
            self.l = 43
            self.f = 6364136223846793005
            self.num_bits = 64
        else:
            (self.w, self.n, self.m, self.r) = (32, 624, 397, 31)
            self.a = 0x9908B0DF
            (self.u, self.d) = (11, 0xFFFFFFFF)
            (self.s, self.b) = (7, 0x9D2C5680)
            (self.t, self.c) = (15, 0xEFC60000)
            self.l = 18
            self.f = 1812433253
            self.num_bits = 32
        self.MT = [0 for i in range(self.n)]
        self.index = self.n + 1
        self.lower_mask = (1 << self.r) - 1
        self.upper_mask = (1 << self.r)

    def ut(self, num):
        """
        untamper a `num` to give back the internal state register
        """
        def get_bit(number, position):
            if position < 0 or position > self.num_bits - 1:
                return 0
            return (number >> (self.num_bits - 1 - position)) & 1

        def set_bit_to_one(number, position):
            return number | (1 << (self.num_bits - 1 - position))

        def undo_right_shift_xor_and(result, shift_len, andd=-1):
            original = 0
            for i in range(self.num_bits):
                if get_bit(result, i) ^ \
                   (get_bit(original, i - shift_len) &
                        get_bit(andd, i)):
                    original = set_bit_to_one(original, i)
            return original

        def undo_left_shift_xor_and(result, shift_len, andd):
            original = 0
            for i in range(self.num_bits):
                if get_bit(result, self.num_bits - 1 - i) ^ \
                   (get_bit(original, self.num_bits - 1 - (i - shift_len)) &
                        get_bit(andd, self.num_bits - 1 - i)):
                    original = set_bit_to_one(original, self.num_bits - 1 - i)
            return original
        num = undo_right_shift_xor_and(num, self.l)
        num = undo_left_shift_xor_and(num, self.t, self.c)
        num = undo_left_shift_xor_and(num, self.s, self.b)
        num = undo_right_shift_xor_and(num, self.u, self.d)
        return num

    def tamper_state(self, y):
        """
        Tamper the state of z3.BitVec(32) y and return the tampered
        thing
        """
        y = y ^ (LShR(y, self.u) & self.d)
        y = y ^ ((y << self.s) & self.b)
        y = y ^ ((y << self.t) & self.c)
        y = y ^ (LShR(y, self.l))
        return y

    def untamper_sat(self, num):
        """
        Same as self.ut but using z3 to find out the state (way slower
        than self.ut)
        """
        S = Solver()
        y = BitVec('y', self.num_bits)
        y = self.tamper_state(y)
        S.add(num == y)
        if S.check() == sat:
            m = S.model()
            return m[m.decls()[0]].as_long()

    def clone_state(self, outputs):
        """
        Clone the internal state given 624, 32-bit outputs
        """
        assert len(outputs) == 624, "To clone full state, 624 outputs needed"
        return list(map(self.ut, outputs))

    def get_seed_mt(self, outputs):
        """
        recovering the initializing knowing some `outputs`
        outputs: list of (output_num, output) pairs 
        (can recover seed with just three consecutive outputs)
        """
        STATE = [BitVec(f'MT[{i}]', self.num_bits) for i in range(self.n + 1)]
        SEED = BitVec('seed', self.num_bits)
        STATE[0] = SEED
        for i in range(1, self.n):
            temp = self.f * \
                (STATE[i - 1] ^ (LShR(STATE[i - 1], (self.w - 2)))) + i
            STATE[i] = temp & ((1 << self.w) - 1)
        self.twist_state(STATE)
        t_start = time()
        S = Solver()
        for index, value in outputs:
            S.add(STATE[index] == self.ut(value))
        if S.check() == sat:
            m = S.model()
            print("time taken :", time() - t_start)
            return m[m.decls()[0]].as_long()
        else:
            print(time() - t_start)

    def twist_state(self, MT):
        """
        Twist an array MT[z3.BitVec(32)] (or 64-bit) to give out the
        next set of internal state registers
        """
        for i in range(self.n):
            x = (MT[i] & self.upper_mask) + \
                (MT[(i + 1) % self.n] & self.lower_mask)
            xA = LShR(x, 1)
            xA = If(x & 1 == 1, xA ^ self.a, xA)
            MT[i] = simplify(MT[(i + self.m) % self.n] ^ xA)


    def untwist(self, outputs):
        """
        Recover the state post twisting
        (can get only the MSB of first element of the internal state)
        """
        MT = [BitVec(f'MT[{i}]', 32) for i in range(self.n)]
        self.twist_state(MT)
        s = Solver()
        for i in range(len(outputs)):
            s.add(outputs[i] == MT[i])
        if s.check() == sat:
            model = s.model()
            untwisted = {str(i): model[i].as_long() for i in model.decls()}
            untwisted = [untwisted[f'MT[{i}]'] for i in range(624)]
            return untwisted
        


class BreakerPy(Breaker):
    """
    Breaker for python functionalities
    """
    def __init__(self):
        Breaker.__init__(self)

    def get_ith(self, outputs, both=False):
        """
        Get i'th output given i-624, i-623 and i-227 th inputs
        if `both=True` then can be only i-623 and i-227 only as 
        we only need MSB of i-624 which can be two possibilities
        """
        if both:
            return self.get_ith([2**31]+outputs),self.get_ith([0]+outputs)
        i_min_624, i_min_623, i_min_227 = map(self.ut, outputs)
        x = (i_min_624 & self.upper_mask) + \
            (i_min_623 & self.lower_mask)
        xA = x >> 1
        if (x % 2) != 0:
            xA = xA ^ self.a
        y = i_min_227 ^ xA
        y = y ^ ((y >> self.u) & self.d)
        y = y ^ ((y << self.s) & self.b)
        y = y ^ ((y << self.t) & self.c)
        y = y ^ (y >> self.l)
        return y & ((1 << self.w) - 1)

    def get_32_bit_seed_python(self, outputs):
        """
        get the seed value if initialzed by a 32 bit seed using
        624 outputs
        """
        MT_init = MT19937(19650218).MT
        MT = [BitVec(f'MT[{i}]', 32) for i in range(self.n)]
        for i in range(self.n):
            MT[i] = BitVecVal(MT_init[i], 32)
        SEED = BitVec('seed', 32)
        i = 1
        for k in range(self.n):
            MT[i] = (MT[i] ^ ((MT[i - 1] ^ LShR(MT[i - 1], 30)) * 1664525)) + SEED
            i += 1
            if i >= self.n:
                MT[0] = MT[self.n - 1]
                i = 1
        for k in range(self.n - 1):
            MT[i] = (MT[i] ^ ((MT[i - 1] ^ LShR(MT[i - 1], 30)) * 1566083941)) - i
            i += 1
            if i >= self.n:
                MT[0] = MT[self.n - 1]
                i = 1
        MT[0] = BitVecVal(0x80000000, 32)
        untwisted = self.untwist(list(map(self.ut,outputs)))
        print(untwisted)
        t_start = time()
        S = Solver()
        S.add([i == j for i, j in zip(MT, untwisted)])
        if S.check() == sat:
            m = S.model()
            print("time taken :", time() - t_start)
            return m[m.decls()[0]].as_long()

    def get_seeds_python(self, outputs, num_seeds=5):
        """
        recovering seeds of python knowing the number of bits 
        in the seed, takes increasing time per increased number of words
        """
        MT_init = MT19937(19650218).MT
        MT = [BitVec(f'MT[{i}]', 32) for i in range(self.n)]
        for i in range(self.n):
            MT[i] = BitVecVal(MT_init[i], 32)
        SEEDS = [BitVec(f'seed[{i}]', 32) for i in range(num_seeds)]
        i,j = 1,0
        for k in range(self.n):
            MT[i] = (MT[i] ^ ((MT[i - 1] ^ LShR(MT[i - 1], 30)) * 1664525)) + SEEDS[j] + j
            i += 1
            j +=1
            if i >= self.n:
                MT[0] = MT[self.n - 1]
                i = 1
            if j==num_seeds:
                j=0
        for k in range(self.n - 1):
            MT[i] = (MT[i] ^ ((MT[i - 1] ^ LShR(MT[i - 1], 30)) * 1566083941)) - i
            i += 1
            if i >= self.n:
                MT[0] = MT[self.n - 1]
                i = 1
        MT[0] = BitVecVal(0x80000000, 32)
        untwisted = self.untwist(list(map(self.ut,outputs)))
        print(untwisted)
        t_start = time()
        S = Solver()
        S.add([i == j for i, j in zip(MT, untwisted)])
        if S.check() == sat:
            m = S.model()
            print("time taken :", time() - t_start)
            recovered = { str(i):m[i].as_long() for i in m.decls() }
            recovered = [ recovered[f'seed[{i}]'] for i in range(num_seeds) ]
            return recovered

    def get_seeds_python_fast(self,outputs):
        """
        recover seed of any size using 624 outputs in python
        runtime independent of seed size, takes anything around
        200-700 seconds
        """
        start_time = time()
        MT = [BitVec(f'MT[{i}]',32) for i in range(624)]
        i=2
        for k in range(self.n - 1):
            MT[i] = (MT[i] ^ ((MT[i - 1] ^ LShR(MT[i - 1], 30)) * 1566083941)) - i
            i += 1
            if i >= self.n:
                MT[0] = MT[self.n - 1]
                i = 1
        MT[0] = BitVecVal(0x80000000, 32)
        untwisted = self.untwist(list(map(self.ut,outputs)))
        S = Solver()
        for i in range(1,self.n):
            S.add(untwisted[i]==MT[i])
        if S.check()==sat:
            m = S.model()
            mt_vals = {str(i):m[i].as_long() for i in m.decls()}
            mt_intermediate = [mt_vals[f'MT[{i}]'] for i in range(1,624) ]

        MT_init = MT19937(19650218).MT
        MT = [BitVecVal(i,32) for i in MT_init]
        SEEDS = [BitVec(f'seed[{i}]', 32) for i in range(624)]
        i,j = 1,0
        for k in range(self.n):
            MT[i] = (MT[i] ^ ((MT[i - 1] ^ LShR(MT[i - 1], 30)) * 1664525)) + SEEDS[j] + j
            i += 1
            j +=1
            if i >= self.n:
                MT[0] = MT[self.n - 1]
                i = 1
            if j==self.n:
                j=0
        S = Solver()
        for i in range(1,self.n):
            S.add(mt_intermediate[i-1]==MT[i])
        if S.check() == sat:
            print("time taken :",time()-start_time)
            m = S.model()
            recovered = { str(i):m[i].as_long() for i in m.decls() }
            recovered = [ recovered[f'seed[{i}]'] for i in range(len(m.decls())) ]
            
            slen = seed_arr_len(recovered)
            seed_arr = [recovered[i]+ slen*(i//slen) for i in range(624)]
            if slen==1:
                return seed_arr[2]
            return seed_arr[slen:slen+2]+seed_arr[2:slen]

    def state_recovery_rand(self, outputs):
        """
        state recovery using python random.random() calls
        """
        MT = [BitVec(f'MT[{i}]',32) for i in range(624)]
        values = []
        for i in outputs:
            values.extend( divmod(int(i*2**53),2**26))
        start_time = time()
        S = Solver()
        for i in range(len(values)):
            if i%624==0:
                self.twist_state(MT)
            S.add(LShR(self.tamper_state(MT[i%624]),5+(i&1))==values[i])
        if S.check()==sat:
            print("time taken :",time()-start_time)
            model = S.model()
            mt = {str(i): model[i].as_long() for i in model.decls()}
            mt = [mt[f'MT[{i}]'] for i in range(len(model))]
            return mt
        
    def int_to_array(self,k):
        if k==0:
            return [0]
        k_byte = int.to_bytes(k,(k.bit_length()+7)//8,'little')
        k_arr = [k_byte[i:i+4] for i in range(0,len(k_byte),4)]
        return [int.from_bytes(i,'little') for i in k_arr ]

    def array_to_int(self,arr):
        arr_bytes  = b"".join([int.to_bytes(i,4,'little') for i in arr])
        return int.from_bytes( arr_bytes ,'little')


