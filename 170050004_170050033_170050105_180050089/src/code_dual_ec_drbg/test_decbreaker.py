from dual_ec_drbg import *
import os

def test_dec(times):
    for i in range(times):
        bts = urandbits(256)
        d = Breaker(bts)
        _, b = d.break_dec()
        print(b)

test_dec(4)