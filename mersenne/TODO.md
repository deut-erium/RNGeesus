1. recovering state of python mersenne twister for random.random() calls (skips 5 or 6 bits from each element of MT array)
    - How much outputs do we need to recover the MT array (state)
2. recovering MT array for getrandbits(n) when n!=32 same as point 1
3. python seed recovery, possible? trying it too
4. [python random](https://github.com/python/cpython/blob/b8fde8b5418b75d2935d0ff93b20d45d5350f206/Modules/_randommodule.c)


5. DUAL_EC_DRBG kleptographic backdoor

6. LCG
7. Truncated LCG
8. LSFR

https://dspace.cvut.cz/bitstream/handle/10467/69409/F8-BP-2017-Molnar-Richard-thesis.pdf?sequence=-1&isAllowed=y

https://www.ambionics.io/blog/php-mt-rand-prediction