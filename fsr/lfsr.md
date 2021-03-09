## Linear Feedback Shit Registers - LFSRs
Linear Feeback shift registers are one of the easiest way to generate seemingly random bits from known bits. The word linear suggests that the algorithm is linear in nature, as in, the next output bit depends linearly on previous bit(s).

## Algorithmic Details

LFSR were used for stream ciphers and are still used today in algorithms like A5/1, A5/2 for providing over-the-air communication in GSM phones. They are also used for pattern generation for exhaustive testing.

### LFSR algorithm:
Linear Feedback Shift Register is a shift register whose input bit is a function of it's previous state.\
Let's say we want to generate random bits using n-bit LFSR with seed $S = S_1, S_2, ..., S_n$. Some fixed values lets say $a_1, a_2,... , a_n$ are used to compute the next input bit (these $a_i$'s are also called feedback polynomial or combination coefficients). Next input bit is given by the formula: 
$$S_{n} = \sum_{j=1}^{j=n}a_j*S_{n-j}$$
where\
$n$ is total numebr of input bits\
$S_i$ is the input bits being used\
$a_j$ is the value of the cofficient, which in our case is 0 or 1 as all the calculation is modulo 2 

After computing the next bit ($S_N$), the rightmost bit of the input is considered as *output bit* and the newly computed bit is attached to the left and the process is repated **k** number of times (called steps or cycles). This gives the *$k$* bits of output which is than used in cryptographic systems or for some other purpose.
<!-- ![](16bit-lfsr.png) -->

During every such step, the $a_i$'s remain the same, as the next input bit is linearly dependent on the previous bits thorough $a_i$'s. Given sufficient number of output bits, we can successfully guess a minimal feedback polynomial and minimal size of LFSR such that we can successfully generate the same random bits.

So we need to find the $a_i$'s and the minimal value $n$. This is where Berlekamp - Massey algorithm comes into play. It's helps in finding $a_i$'s in $O(n^2)$.

## Background
LFSRs were used in stream ciphers in early years of internet. Later on, Berlekamp published a paper that talked about an algorithm to decode [BCH codes](https://en.wikipedia.org/wiki/BCH_code). Later on, James Massey recognized its application to LFSRs and simplified the algorithm.

### Berlekamp – Massey Algorithm
The Berlekamp–Massey algorithm is an algorithm that will find the shortest linear feedback shift register (LFSR) for a given binary output sequence. \
This algorithm starts with the assumption that the length of the LSFR is $l = 1$, and then *iteratively* tries to generate the known sequence and if it succeeds, everything is well, if not, $l$ must be *increased*. 

To solve a set of linear equations of the form $S_i+v+\Lambda_1S_i+ v-1 + ... + \Lambda_{v-1}S_{i + 1} + \Lambda_vS_i=0$, a potential instance of $\Lambda$ is constructed step by step. Let $C$ denote such a potential candidate, it is sometimes also called the "connection polynomial" or the "error locator polynomial" for L errors, and defined as $C = c_LX_L + c_{L-1}X_{L-1} + ...+ c_2X_2 + c_1X + 1$. The goal of the Berlekemp - Massey algorithm is to now determine the minimal degree $L$ and to construct $C$ such, that $S_n+c_1S_{n-1} + ... + c_LS_{n-L}= 0, \forall L\leq n \leq N-1$, where $N$ is the total number of syndromes, and $n$ is the index variable used to loop through the syndromes from $0\ to\ N-1$.

With each iteration, the algorithm calculates the **discrepancy** $d$ between the candidate and the actual feedback polynomial: 
$$ d = S_k+c_1S_{k-1}+ ... + c_LS_{k-L} $$
If the discrepancy is zero, the algorithm assumes that $C$ is a valid candidate and continues. Else, if $d\neq0$, the candidate $C$ must be adjusted such, that a recalculation of the discrepancy would lead to $d = 0$. This re-adjustments is done as follows: 
$$ C= C - (d/b)X^mB $$
where,\
$B$ is a copy of the *last candidate* $C$ since $L$ was updated,\
$b$ a copy of the *last discrepancy* since $L$ was updated,\
and the multiplication by X^m is but an index shift. \
The new discrepancy can now easily be computed as $d = d-(d/b)b = d-d = 0$. This above algorithm can further be simplified for modulo 2 case. See Wiki.

### Geffe Generator
The Geffe generator consists of three LFSRs: LFSR-1, LFSR-2 and LFSR-3 using primitive feedback polynomials. If we denote the outputs of these registers by $x_1$, $x_2$ and $x_3$, respectively, then the Boolean function that combines the three registers to provide the generator output is given by
$$ F(x_1, x_2, x_3) = (x_1 \land x_2) \oplus (\lnot x_1 \land x_2) $$
There are $2^3 = 8$ possible values for the outputs of the three registers, and the value of this combining function for each of them is shown in the table below: 
|  |  |  |  |  |  |  |  |  | 
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| $x_1$ | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 | 
| $x_2$ | 0 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 
| $x_3$ | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 
| $F(x_1, x_2, x_3)$ | 0 | 1 | 0 | 1 | 0 | 0 | 1 | 1 | 

Consider the output of the third register, $x_3$. The table above makes it clear that of the 8 possible outputs of $x_3$, 6 of them are equal to the corresponding value of the generator output, $F(x_1, x_2, x_3)$, i.e. $x_3 = F( x_1, x_2,x_3 )$ in **75%** of all possible cases. This correlation can be exploited to have a brute force search on the key bits of LSFR-3, since on correct key we would observe an advantage of 1/4 over any other key.
This will reduce the brute forcing of 3 LFSRs to only 2 LFSRs!!\
Similarly LFSR-2 has 75% of correlation with the output. After recovering key bits of LFSR-2 and LFSR-3, we can brute force over the bits of LFSR-1 to get the correct key.\
Thus, burte search of $2^{n_1 + n_2 + n_3}$ is reduced to $2^{n1} + 2^{n2} + 2^{n3}$.

## Our work
### Modeling
Though Geffe generator is susceptible to correlation attack, we found a faster and much more efficient attack by modeling the Geffe generator as a Boolean formula over the key-bits and solving the satisfiability problem over generated output bits.  
We encoded the **seed** finding problem into [Z3Prover](https://github.com/Z3Prover/z3).  
We abstracted the Geffe Generated in a manner that it would return a boolean function when the input key bits are a boolean function and and python bool when the input bits are python bool.
[See solve function](#geffe-generator-implementation) and [Test cases](#testing-the-implementation).

### Results
We observed significantly faster runtimes using the Z3 boolean model as compared to brute force correlation attack.

| Specifications | Time taken using brute-force| Time taken using Z3 solver |
|:---| :---: | :---: |
| 10-bit seed each, 128 bit output | 01.50s | 0.25s |
| 12-bit seed each, 128 bit output | 06.25s | 0.26s |
| 12-bit seed each, 256 bit output | 12.50s | 0.41s |
| 12-bit seed each, 350 bit output | 16.62s | 0.54s |
| 13-bit seed each, 256 bit output | 25.08s | 0.74s |
| 14-bit seed each, 256 bit output | 52.30s | 1.12s |
| 16-bit seed each, 256 bit output | 222.66s | 4.53s |
| 16-bit seed each, 512 bit output | 449.29s | 5.99s |
| 18-bit seed each, 256 bit output | 936.59s | 29.33s |
<!-- | 24-bit seed each, 2048 bit output | - Timout - | 400.45s | -->
While the runtime of discovered correlation attack is observably *exponential* in the number of bits of LFSRs whereas, observed runtime of our approach is *subexponential/polynomial*, since boolean constraints are relatively sparse and SAT solvers are highly optimized in solving such boolean constraints.  

### Limitations
#### Berlekamp-Massey VS SAT modeling
For finding the minimum degree feedback polynomial using SAT encoding, we ran into the problem of not knowing the degree of the polynomial, thus we need to enumerate over the possible guesses of the degree and checking the satisfiability of the generated boolean formula over LSFR state bits and output bits.

Since we have no expected bounds on runtimes, we could not conclude termination while recovering the minimal polynomial using the SAT encoding approach.

### Future Scope
We explored a known weak combiner generator where the correlations between various LFSR bits and the generated output bits is obvious, the solver might me internally exploiting some higher order correlation which might be difficult to discover.

This approach can be extended to different combiner generators and seemingly undiscovered correlations can be expoilted in a similar efficient way.

## References: 
- [Wikipedia - Berlekamp Massey Algorithm](https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm)
- [Bell0bytes blog](https://bell0bytes.eu/linear-feedback-shift-registers/)
- [BMA - IIT Kharagpur](https://cse.iitkgp.ac.in/~debdeep/courses_iitkgp/Crypto/slides/BMA.pdf)
- [Wikipedia - Correlation Attack](https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm)

## Code:

### LFSR: 
- This class is an abstarct implementation handeling `pythonic` as well as `Z3-type` input of seed and feedback/combination polynomial.
```python
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
```
### Geffe Generator Implementation:
- This class is also an abstraction. Depending on our utility we can use the same instance to guess the seeds generated using same combination polynomial.
```python
class Geffe:
    """ Geffe Generator's with solver in z3 and brute force as well. We need to know the combination polynomial beforehand """

    def __init__(self, l1, l2, l3, c1, c2, c3):
        self._l1, self._l2, self._l3 = l1, l2, l3
        self._c1, self._c2, self._c3 = c1, c2, c3
        self._lfsrs = [LFSR(self._l1,c1), LFSR(self._l2, c2), LFSR(self._l3, c3)]
    
    def next_bit(self):
        bits = [lfsr.next_bit() for lfsr in self._lfsrs]
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
        for seed2 in tqdm.tqdm(itertools.product('01', repeat=len(lfsr2._comb_poly))):
            lfsr2.set_seed(list(map(int,seed2)))
            x3 = lfsr2.get_lfsr(len(opt))
            corr = sum(x==y for x,y in zip(opt, x3))
            if corr >= int(0.70*n):
                possible_seeds2.append(''.join(seed2))
            if m2 < corr:
                m2 = corr
        assert len(possible_seeds2) >=1, "Error: No x3 found, less data supplied."

        # > 75% match of opt with x2 i.e. lfsr1
        possible_seeds1 = []
        m1 = 0
        for seed1 in tqdm.tqdm(itertools.product('01', repeat=len(lfsr1._comb_poly))):
            lfsr1.set_seed(list(map(int,seed1)))
            x2 = lfsr1.get_lfsr(len(opt))
            corr = sum(x==y for x,y in zip(opt, x2))
            if corr >= int(0.70*n):
                possible_seeds1.append(''.join(seed1))
            if m1 < corr:
                m1 = corr
        assert len(possible_seeds1) >=1, "Error: No x2 found, less data supplied"
        
        candidates = [(x, y) for x in possible_seeds1 for y in possible_seeds2]

        return (possible_seeds1[0], possible_seeds2[0])
```

### UnLFSR (Z3):
- Finding the seed of LFSR given some output bits. We need to guess the number of bits in the seeds though.
```python
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
            print("ERROR: unsolvable... provide more output or bits guessed is wrong")
```
### Testing the implementation
- This code implemtes certain test cases for testing the Z3 modeling
```python
import time, random
from lfsr import *

def str_to_lst_int(string):
    return list(map(int, string))

def test_n_bit_k_steps(n: int, k: int):
    """ Generate seed and combination polynomial and generate some eandom bits """
    rndm_seed = bin(random.getrandbits(n))[2:]
    seed = rndm_seed + '0'*(n-len(rndm_seed))
    rndm_poly = bin(random.getrandbits(n))[2:]
    feedback_poly = rndm_poly + '0'*(n - len(rndm_poly))
    lfsr = LFSR(list(map(int,seed)), list(map(int,feedback_poly)))
    gen_opt = lfsr.get_lfsr(2*k)

    # Test bm algo
    start_time = time.time()
    bm = Berlekamp_Massey(gen_opt[:len(gen_opt)//2])
    print("Time taken to recover LFSR seed: ", time.time() - start_time)
    sd = bm.get_seed()
    taps = bm.get_taps()

    # UnLFSR test
    unlfsr = UnLFSR_Z3(gen_opt[:len(gen_opt)//2], len(gen_opt)//4)
    start_time = time.time()
    seed_unlfsr, taps_unlfsr = unlfsr.solve()
    print("Unlfsr solved in ", time.time() - start_time)
    
    lfsr_new = LFSR(sd, taps)
    bm_opt = lfsr_new.get_lfsr(2*k)

    if bm_opt == gen_opt:
        print(f"No mismatch for {n} bit seed. Matched {2*k} (random) output bits")
        print("Success!")
    else:
        for i, j in enumerate(zip(bm_opt[k:], gen_opt[k:])):
            if j[0] != j[1]:
                print(f"For {2*k} bits, 1st mismatch at index: {i}")
                print("Partial Success. Need more output bits")
                break
    return

test_n_bit_k_steps(2048,4096)

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
    geffe = Geffe(l1, l2, l3, list(map(int,c1)), list(map(int,c2)), list(map(int,c3)))
    start_time = time.time()
    for l1_z3, l2_z3, l3_z3 in geffe.solve(opt):
        print("Time taken Geffe using Z3: " , time.time() - start_time)
        print(ll1 == ''.join(map(str,l1_z3)), ll2 == ''.join(map(str,l2_z3)), ll3 == ''.join(map(str,l3_z3)))
        start_time = time.time()

    start_time = time.time()
    l2_normal, l3_normal = geffe_normal.solve_bruteforce(opt)
    print("Time taken to break Geffe using bruteforce: ", time.time() - start_time)
    print((ll2 == l2_normal) & (ll3 == l3_normal))

test_geffe_generator(2048, 24)
```