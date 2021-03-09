# Dual EC DRBG - Kleptographic Backdoor

Dual EC(Elliptic Curve) DRBG(Deterministic Random Bit Generator) was a pseudo-random number generator which was a one of its kind of PRNGS in that it relied on elliptic curves, while most earlier PRNGs relied on symmetric components like block ciphers or hash functions. As a result, it was quite slow than it's competitive RNGs. However, since it was based on elliptic curves, it was believed that this was a **Cryptographically Secure** PRNG; however it later turned out that it suffered from many flaws. Despite all its shortcomings, it was highly standardized by NSA and was one of the four PRNGs proposed by NIST as **approved CSRNGs**. The main issue was that **it had no security proofs**!

## Algorithmic Details

![](./decdrbg.png)

A particular curve $C$ and two points on the curve $P$ and $Q$ are chosen apriori and remain fixed throughout.

The RNG is initially seeded with a random seed $s_0$. Whenever the RNG is asked for an output, assuming that the current seed(state) is $s_i$; $r_i \leftarrow (s_i * P) |_{x}$ is computed where $M |_{x}$ denotes the $x$-coordinate of the point $M$.

The seed is then updated as $s_{i+1} \leftarrow (r_i * P) |_{x}$, $t_i$ is computed as $t_i \leftarrow (r_i * Q) |_{x}$ and then the output generated is $LSB_{bitlen - 16}(t_i)$ where $LSB_{k}$ denotes the Least Significant $k$ bits.

We work with the case where $C$ is the NIST $P-256$ curve, in which case $bitlen = 256$ and therefore, in a single round, $240$ bits are obtained.


### Backdoor

Since the RNG only throws out the top $16$ bits, which is a very small amount compared to the total bit-length being $256$, these $16$ bits can be bruteforced to obtain list of points on the curve admitting the mentioned `240` least significant bits ensure. This would usually be the end of any possible attack, however this RNG had another issue!  

If $P$ and $Q$ were chosen to be totally random points, being able to deduce $r_i * Q$ doesn't help us in recovering the internal state at all. However, if $P$ and $Q$ were deliberately chosen such that $Q = e*P$ or $P = d * Q$ where either of $e$ or $d$ is known to the attacker; given the generated point $r_i * Q$, one could obtain the next state by computing the internal point as $r_i * P = r_i * d * Q = (r_i * Q) * d$. This not only compromises the security of the RNG, but also allows the attacker to be able to predict all future outputs of the RNG.

## Background

It was [shown](https://eprint.iacr.org/2006/190.pdf) that the generated `240` bits are not indistinguishable from truly random bits but actually admit non-negligible bias; thus demonstrating that the generator is insecure! Even though [some papers](https://eprint.iacr.org/2007/048.pdf) tried to show that DEC-DRBG was based on `cryptographically hard` problems; the possibility of a [kleptographic backdoor](https://rump2007.cr.yp.to/15-shumow.pdf) was later demonstrated. With this backdoor, one could essentially break down TLS just by monitoring one encryption connection.

## Our Work

In our Proof-of-Concept, we demonstrate that choosing $P$ and $Q$ of our own accord(that is, if we by chance know $e$ or $d$) allows us to recover the internal state of the RNG in approximately 32 bytes(256 bits) of the output. This confirms the possibility of there being a backdoor in the RNG, and hence, allowing the attacker to compromise encryptions which relied on Dual-EC-DRBG as the random number generator.

## Challenges/Difficulties

How to perform mathematical operations on elliptic curves efficiently in python was one of the small challenges we encountered.  
Another challenge is that one can not demonstrate the backdoor on the original values of $P$ and $Q$ chosen by the NSA implementation, since it's an extremely hard problem to find $e$ given points on curve $P$ and $Q$.

## Assumptions

We demonstrated the backdoor by choosing our own random `multiplicative relation` between the "random points" $P$ and $Q$.

## Critique/Comparison

If only one output(`240` bits) can be obtained from the RNG, $~2^15$ possible states exist; thus, the attack doesn't work in such a case(atleast around `256` bits need to be seen, which essentially means two runs of the RNG).

The values of $P$ and $Q$ which were used in the actual implementation had been publicised by NSA to be the ones which allowed "fast computations", nobody knows where these values actually came from! Since the values were chosen by themselves, it is unknown whether they actually had utilized this backdoor but the existence of a backdoor in a popular PRNG is very troublesome to the cryptographic community in itself.

## References
 - [A blog post](https://blog.cryptographyengineering.com/2013/09/18/the-many-flaws-of-dualecdrbg/)
 - https://eprint.iacr.org/2006/190.pdf  
 Cryptanalysis of the Dual Elliptic Curve Pseudorandom Generator  
 Berry Schoenmakers and Andrey Sidorenko  
 Cryptology ePrint Archive, Report 2006/190
 - [Possibility of Backdoor](https://rump2007.cr.yp.to/15-shumow.pdf)
 - https://eprint.iacr.org/2007/048.pdf  
 A Security Analysis of the NIST SP 800-90 Elliptic Curve Random Number Generator  
Daniel R. L. Brown and Kristian Gjøsteen  
Cryptology ePrint Archive, Report 2007/048