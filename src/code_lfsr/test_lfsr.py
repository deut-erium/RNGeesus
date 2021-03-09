import time, random
from lfsr import *

def str_to_lst_int(string):
    return list(map(int, string))

def test_n_bit_k_steps(n: int, k: int):
    # Generate seed and combination polynomial and generate some eandom bits
    rndm_seed = bin(random.getrandbits(n))[2:]
    seed = rndm_seed + '0'*(n-len(rndm_seed))
    rndm_poly = bin(random.getrandbits(n))[2:]
    feedback_poly = rndm_poly + '0'*(n - len(rndm_poly))
    lfsr = LFSR(list(map(int,seed)), list(map(int,feedback_poly[::-1])))
    gen_opt = lfsr.get_lfsr(2*k)

    # Test bm algo
    start_time = time.time()
    bm = Berlekamp_Massey(gen_opt)
    print("Time taken to recover LFSR seed: ", time.time() - start_time)
    sd = bm.get_seed()
    taps = bm.get_taps()
    # print(bm.get_degree())
    # UnLFSR test
    # unlfsr = UnLFSR_Z3(gen_opt[:len(gen_opt)//2], len(gen_opt)//4)
    # start_time = time.time()
    # seed_unlfsr, taps_unlfsr = unlfsr.solve()
    # print("Unlfsr solved in ", time.time() - start_time)
    print(feedback_poly)
    print(''.join(map(str,taps)))
    print(seed)
    print(''.join(map(str, sd)), len(sd), bm.get_degree())
    lfsr_new = LFSR(sd, taps[::-1])
    bm_opt = lfsr_new.get_lfsr(2*k)

    if bm_opt == gen_opt:
        print(f"No mismatch for {n} bit seed. Matched {2*k} (random) output bits")
        print("Success!")
    else:
        for i, j in enumerate(zip(bm_opt, gen_opt)):
            if j[0] != j[1]:
                print(f"For {2*k} bits, 1st mismatch at index: {i}")
                print("Partial Success. Need more output bits")
                break
    return

# test_n_bit_k_steps(256,2000)

# Test Geffes generator
def test_geffe_generator(num_opt_bits, size_taps):
    """ Given n output bits and taps of all the 3 LFSRs, find the actual seeds of LFSRs """
    # c1 = bin(random.getrandbits(size_taps))[2:]
    # c1 = c1 + '0'*(size_taps - len(c1))
    c1 = '11011' + '0'*(size_taps - 5)
    # c2 = bin(random.getrandbits(size_taps))[2:]
    # c2 = c2 + '0'*(size_taps - len(c2))
    c2 = '11100001' + '0'*(size_taps - 8)
    # c3 = bin(random.getrandbits(size_taps))[2:]
    # c3 = c3 + '0'*(size_taps - len(c3))
    c3 = '1110000111' + '0'*(size_taps - 10)

    # Z3 Ls
    l1 = [BitVec(f'l1_{i}',1) for i in range(len(c1))]
    l2 = [BitVec(f'l2_{i}',1) for i in range(len(c2))]
    l3 = [BitVec(f'l3_{i}',1) for i in range(len(c3))]

    # Normal l1, l2 l3
    ll1 = bin(random.getrandbits(size_taps))[2:]
    ll1 = ll1 + '0'*(size_taps - len(ll1))
    ll2 = bin(random.getrandbits(size_taps))[2:]
    ll2 = ll2 + '0'*(size_taps - len(ll2))
    ll3 = bin(random.getrandbits(size_taps))[2:]
    ll3 = ll3 + '0'*(size_taps - len(ll3))

    geffe_normal = Geffe(str_to_lst_int(ll1), str_to_lst_int(ll2), str_to_lst_int(ll3),  str_to_lst_int(c1), str_to_lst_int(c2), str_to_lst_int(c3))

    opt = geffe_normal.get_seqn(num_opt_bits)
    # opt = opt + '0'*(num_opt_bits - len(opt))
    geffe = Geffe(l1, l2, l3, list(map(int,c1)), list(map(int,c2)), list(map(int,c3)))
    start_time = time.time()
    for l1_z3, l2_z3, l3_z3 in geffe.solve(opt):
        print("Time taken Geffe using z3: " , time.time() - start_time)
        print(ll1 == ''.join(map(str,l1_z3)), ll2 == ''.join(map(str,l2_z3)), ll3 == ''.join(map(str,l3_z3)))
        start_time = time.time()

    start_time = time.time()
    l2_normal, l3_normal = geffe_normal.solve_bruteforce(opt)
    print("Time taken to break Geffe using bruteforce: ", time.time() - start_time)
    # print(ll2, l2_normal, ll3, l3_normal)
    print((ll2 == l2_normal) & (ll3 == l3_normal))

test_geffe_generator(200, 16)
