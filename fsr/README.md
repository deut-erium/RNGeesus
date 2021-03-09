# Linear Feedback Shit Registers - LFSRs
Linear Feeback shift registers are one of the easiest way to generate seemingly random bits from known bits. The word linear suggests that the algorithm is linear in nature, as in, the next output bit depends linearly on the current bits.

A simple way to write LFSR is:\
$$S_{n} = \sum_{j=1}^{j=n}a_j*S_{n-j}$$
where\
$n$ is total numebr of input bits\
$S_i$ is the input bits being used\
$a_j$ is the value of the cofficient, which in our case is 0 or 1 as all the calculation is modulo 2

This type of approach does generate seeming random bits faster, but due to its **linear** nature, its cryptanalysis becomes easy.

### LFSR algorithm:
Let's say we are using n-bit LFSR $S_1, S_2, ..., S_n$ called seed. There are some fixed values of $a_1, a_2,... , a_n$ which are used to compute the next input bit (these $A_i$'s are also called feedback polynomial). After computing the next bit, the rightmost bit (MSB) is popped of input is considered as *output* and the newly computed bit is attached to the left (LSB) and the process is reated **k** number of times (or steps or cycles). This gives the *$k$* bits of output which is than used in cryptographic systems.

During every such step, the $a_i$'s remain the same, as the next out is linearly dependent on the inputs, given sufficient numeber of output bits, we can successfully guess the initial n-bit and the feedback polynomial and hence obtain all the subsequent output bits.

So we need to find the $a_i$'s and the value $n$. This is where Berlekamp - Massey algorithm comes into play. It's helps in finding $S_i$'s in $O(n^2)$ time. 

## Berlekamp – Massey Algorithm
This algorithm starts with the assumption that the length of the LSFR is $l = 1$, and then *iteratively* tries to generate the known sequence and if it succeeds, everything is well, if not, $l$ must be *increased*. 

To solve a set of linear equations of the form $S_i+ν+Λ_1S_i+ ν−1 + ... + Λ_{ν−1}S_{i + 1} + Λ_νS_i=0$, a potential instance of $Λ$ is constructed step by step. Let $C$ denote such a potential candidate, it is sometimes also called the "connection polynomial" or the "error locator polynomial" for L errors, and defined as $C = c_LX_L + c_{L−1}X_{L−1} + ...+ c_2X_2 + c_1X + 1$. The goal of the Berlekemp - Massey algorithm is to now determine the minimal degree $L$ and to construct $C$ such, that $S_n+c_1S_{n−1} + ... + c_LS_{n−L}= 0, \forall L≤n≤N−1$, where $N$ is the total number of syndromes, and $n$ is the index variable used to loop through the syndromes from 0 to N−1.

With each iteration, the algorithm calculates the **discrepancy** $d$ between the candidate and the actual feedback polynomial: 
$$ d = S_k+c_1S_{k−1}+ ... + c_LS_{k−L} $$
If the discrepancy is zero, the algorithm assumes that $C$ is a valid candidate and continues. Else, if $d≠0$, the candidate $C$ must be adjusted such, that a recalculation of the discrepancy would lead to $d = 0$. This re-adjustments is done as follows: 
$$ C= C− (d/b)X^mB $$
where,\
$B$ is a copy of the *last candidate* $C$ since $L$ was updated,\
$b$ a copy of the *last discrepancy* since $L$ was updated,\
and the multiplication by X^m is but an index shift. \
The new discrepancy can now easily be computed as $d = d−(d/b)b = d−d = 0$. This above algorithm can further be simplified for modulo 2 case. See Wiki.

Sauce: 
- [Wikipedia](https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm)
- [Blog Post, Seems correct idk](https://bell0bytes.eu/linear-feedback-shift-registers/)
- [BMA - IIT Kharagpur](https://cse.iitkgp.ac.in/~debdeep/courses_iitkgp/Crypto/slides/BMA.pdf)