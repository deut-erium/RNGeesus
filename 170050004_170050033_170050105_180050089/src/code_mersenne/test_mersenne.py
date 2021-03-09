import random
from time import time
from mersenne import *
import os

def urandbits(n):
    """
    get os true random bits
    """
    return int.from_bytes(os.urandom(n//8),'big')


def test_seed_mt():
    """
    test seed recovery using 3 outputs in standard MT (also used in numpy)
    takes around 200s
    """
    rand_seed = urandbits(32)
    r = MT19937(rand_seed)
    outputs = [(i,r.extract_number()) for i in range(3)]
    b = Breaker()
    recovered_seed = b.get_seed_mt(outputs)
    assert recovered_seed == rand_seed 
    print("success",recovered_seed)

def test_seed_mt_64():
    """
    test seed recovery using 3 outputs in MT19937-64
    also takes around 200s
    """
    rand_seed = urandbits(64)
    r = MT19937(rand_seed,bit_64=True)
    outputs = [(i,r.extract_number()) for i in range(3)]
    b = Breaker(bit_64=True)
    recovered_seed = b.get_seed_mt(outputs)
    assert recovered_seed == rand_seed 
    print("success",recovered_seed)


def twist(mt_orig):
    """
    helper - mersenne state twist operation
    """
    mt = [i for i in mt_orig]
    um = 0x80000000
    lm = 0x7fffffff
    a = 0x9908B0DF
    n = 624
    m = 397
    for i in range(n):
        x = (mt[i]&um) + (mt[(i+1)%n]&lm)
        xA = x>>1
        if x&1:
            xA=xA^a
        mt[i] = mt[(i+m)%n]^xA
    return mt

def tamper(num):
    """
    tamper operation for 32 bit MT
    """
    u,s,t,b,c,d,l,w,n,m = 11,7,15,0x9D2C5680,0xEFC60000,0xFFFFFFFF,18,32,624,397
    y = num
    y = y^((y>>u)&d)
    y = y^((y<<s)&b)
    y = y^((y<<t)&c)
    y = y^(y>>l)
    return y


def check_ut():
    """
    Check untwist operation to find the state before twist
    recovers all states except the first ( only 1 bit information)
    """
    rand_seed = urandbits(32)
    r = MT19937()
    r.seed_mt(rand_seed)
    untwisted_orig = [i for i in r.MT]
    outputs = [r.extract_number() for i in range(624)]
    b = Breaker()
    untampered = list(map(b.ut,outputs))
    assert list(map(tamper,untampered))==outputs
    assert twist(untwisted_orig)==untampered
    untwisted = b.untwist(untampered)
    assert twist(untwisted)==untampered
    #assert untwisted == untwisted_orig 
    #wont work we only know a single from the first element
    assert untwisted_orig[1:]==untwisted[1:]
    print("success")

def test_recover_32bit():
    rand_seed = urandbits(32)
    r = MTpython(rand_seed)
    b = BreakerPy()
    outputs = [r.extract_number() for i in range(624)]
    recovered_seed = b.get_32_bit_seed_python(outputs)
    print(rand_seed,recovered_seed)
    assert recovered_seed == rand_seed
    print("success")

def test_recover_init_by_array(x):
    rand_seeds = [urandbits(32) for i in range(x)]
    r = MTpython(0)
    r.init_by_array(rand_seeds)
    b = BreakerPy()
    outputs = [r.extract_number() for i in range(624)]
    recovered_seeds = b.get_seeds_python(outputs,x)
    print(rand_seeds)
    print(recovered_seeds)
    assert rand_seeds==recovered_seeds
    print("success")

def int_to_array(k):
    k_byte = int.to_bytes(k,(k.bit_length()+7)//8,'little')
    k_arr = [k_byte[i:i+4] for i in range(0,len(k_byte),4)]
    return [int.from_bytes(i,'little') for i in k_arr ]

def array_to_int(arr):
    return int.from_bytes( b"".join([int.to_bytes(i,4,'little') for i in arr]) ,'little')

def test_python_int_seeds():
    """
    checking init_by_array works as intended with random.seed(integer)
    """
    r = MTpython(0)
    for i in range(1,1000):
        int_seed = urandbits(8*i)
        array_seed = int_to_array(int_seed)
        assert array_to_int(array_seed)==int_seed
        r.init_by_array(array_seed)
        random.seed(int_seed)
        assert r.get_state()==random.getstate()

def test_python_int_seeds2():
    """
    checking random.seed(integer) works as intended
    """
    r = MTpython(0)
    for i in range(1,1000):
        int_seed = urandbits(8*i)
        r.seed(int_seed)
        random.seed(int_seed)
        assert r.get_state()==random.getstate()

def test_python_seed_recovery_fast():
    """
    testing the integer seed recovery in python
    should take anywhere from 200 - 800 s
    """
    seed_len = random.randint(1,624)*32
    rand_seed = urandbits(seed_len)
    rand_seed_arr = int_to_array(rand_seed)
    random.seed(rand_seed)
    outputs = [random.getrandbits(32) for i in range(624)]
    b = BreakerPy()
    seed_arr = b.get_seeds_python_fast(outputs)
    assert seed_arr==rand_seed_arr
    print("success")

def state_recovery_rand():
    """
    state recovery using 
    """
    rand_seed = urandbits(1234)
    random.seed(rand_seed)
    state_orig = list(random.getstate()[1][:-1])
    outputs = [random.random() for i in range(624)]
    b = BreakerPy()
    recovered_state = b.state_recovery_rand(outputs)
    print(sum(i==j for i,j in zip(recovered_state,state_orig)))
    assert recovered_state == state_orig
    print("success")

def compare_ut_sat_ut(num_outs=10000):
    """
    Comparing the performance of Breaker.untamper_sat and
    Breaker.ut i.e. speed comparison of sat and direct
    algorighm
    """
    outputs = [random.getrandbits(32) for i in range(num_outs)]
    b = BreakerPy()
    start_time = time()
    ut1 = list(map(b.ut,outputs))
    time_ut = time()-start_time
    start_time = time()
    ut2 = list(map(b.untamper_sat,outputs))
    time_ut_sat = time()-start_time
    print(f'time taken by ut: {time_ut}, time taken by sat: {time_ut_sat}')
    return time_ut_sat/time_ut




#rand_seed = urandbits(96)
#random.seed(rand_seed)
outputs = [random.getrandbits(32) for i in range(624)]
b = BreakerPy()


#r.init_by_array([0x44434241,0x48474645,0x49])
#random.seed(0x494847464544434241)
#print(r.get_state()==random.getstate())
#r.init_by_array([0x44434241])
#init_state = r.MT.copy()
#r.init_32bit_seed(0x44434241)
#random.seed(0x44434241)
#print(r.get_state()==random.getstate())
#outputs = [r.extract_number() for i in range(3000)]
#untampered = list(map(b.ut,outputs))

#r2 = random.Random()
#r2.setstate(r.get_state())
#initial = r.get_state()
#assert r.get_state() == r2.getstate()
#o1 = [r.extract_number() for i in range(1000)]
#o2 = [r2.getrandbits(32) for i in range(1000)]
#common = list(set(o1)&set(o2))
#test()
#test_64()
#outputs = [random.getrandbits(32) for i in range(3)]
#b = Breaker()
#print(b.get_seed(outputs))
