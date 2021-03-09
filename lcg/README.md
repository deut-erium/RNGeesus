# Linear Congruential Generator

Linear Congruential Generator(LCG) is a method of generating a sequence of pseudo-randomized numbers calculated using modular arithmetic. This method has seen quite widespread usage since the theory behind this method is pretty easy to understand and is also easy to implement as well as fast and require minimal memory, especially on computer hardware. However, having seen widespread usage does not imply security in any manner. In fact, it's recommended that LCGs not be used as random generators at all!

Donald Knuth suggested the usage of Truncated LCG, where only some bits of the internal state were outputted(say upper half bits). These turned out to have much better statistical properties than the original LCGs. However, these are not cryptographically secure either; and indeed there exist attacks which can find out the internal state given a few outputs!

## Algorithmic Details
A linear congruential generator maintains an internal state $s$, which is updated on every call to the generator as $s := (a*s + b) \% m$. The updated state is the generated pseudo-random number. Therefore, the generated numbers $X_i$ follow the recurrence relation $X_{i+1} = (a*X_i + b) \% m$, or, equivalently $X_i = (a^{i} * X_0 + b*(1 + a + \dots + a^{i-1})) \% m$. Note that $0 \le X_i, s, a, b < m$.

For a truncated LCG which outputs certain most significant bits of the internal state, generated number $X$ can be written as $X = (s \gg trunc)$ where $\gg$ denotes logical right-shift and $trunc$ is the number of lower bits to be truncated.

## Background

### LCG
It has been shown that given sufficient number of outputs, the parameters of a secret LCG can be recovered. The algorithm for the same is described next.

Let's first assume that $b = 0$. Then we have $X_{i+1} = a*X_i \% m$ \& $X_{i+2} = a^2 * X_i \% m \implies m | (X_{i+2} * X_i - X_{i+1}^2)$. Thus, $m$ would be a divisor of the GCD of $X_{i+2} * X_i - X_{i+1}^2 , \forall i$. Given more and more outputs, the probability of $m$ being the GCD itself rises. In this manner $m$ can be recovered. Knowing $m$ and $X_i$, $a$ can be recovered as $a = X_{i+1} * X_i^{-1} \% m$.  

Suppose now that $b \ne 0$. Given $X_i$, we define $Y_i = (X_{i+1} - X_i) \% m$. Then $Y_{i+1} = (X_{i+2} - X_{i+1}) \% m = [(a * X_{i+1} + b) - (a * X_i + b)] \% m = [a * (X_{i+1} - X_i)] \% m = (a * Y_i) \% m$. $Y_i$ are therefore LCG with $b = 0$ with the same $a$ and $m$, which can be recovered using the above method. Further $b$ can be recovered as $b = (X_2 - a * X_1) \% m$. Thus, any of the 3 parameters, if unknown, can be recovered as mentioned.

Note that here, $X_1$ denotes the first output and so on; we can then recover the seed as $X_0 = a^{-1} * (X_1 - b) \% m$.

### Truncated LCG
Lattice-based attacks have been described on truncated LCGs to recover the internal states. Here we describe one form of attack in which we're provided some number(say $k$) of generated (truncated) outputs $x_i$, $a$, $b$, $m$ and the truncation(say $t$).

$x_i =  X_i \gg t, X_{i+1} = (a*X_i + b) \% m \implies X_i = 2^t * x_i + y_i$, where $0 \le y_i < 2^t$. Here, $x_i$ are known to us while $y_i$ are the unknowns.

The forthcoming attack is borrowed from this [paper](https://www.math.cmu.edu/~af1p/Texfiles/RECONTRUNC.pdf) on reconstructing truncated integer variables satisfying linear congruences. 


Consider the matrix $L$ defined for some $k$ as -
$$\begin{bmatrix}
    a & -1 & 0 & \dots & 0\\
    a^2 & 0 & -1 & \dots & 0\\
    \vdots & \vdots & \vdots & \ddots & \vdots\\
    a^{k-1} & 0 & 0 & \dots & -1\\
    M & 0 & 0 & \dots & 0\\
    0 & M & 0 & \dots & 0\\
    0 & 0 & M & \dots & 0\\
    \vdots & \vdots & \vdots & \ddots & \vdots\\
    0 & 0 & 0 & \dots & M\\
\end{bmatrix} \implies L \begin{bmatrix}
    X_1\\
    X_2\\
    X_3\\
    \vdots\\
    X_k\\
\end{bmatrix} = \begin{bmatrix}
    b + M \alpha_1\\
    b(1+a) + M \alpha_2\\
    \vdots\\
    b(1+a+\dots+a^{k-2}) + M \alpha_{k-1}\\
    M X_1\\
    M X_2\\
    M X_3\\
    \vdots\\
    M X_k\\
\end{bmatrix} = \begin{bmatrix}
   
    0\\
    0\\
    \vdots\\
    0\\
\end{bmatrix} (\% M)$$
since $X_i = [a^{i-1} * X_1 + b(1 + a + \dots + a^{i-2})] \% M = a^{i-1} * X_1 + b \frac{a^{i-1}-1}{a-1} + M \alpha_{i-1}$ for some $\alpha_{i-1} \in \Z$. Note that here $L$ is a $2k-1 \times k$ lattice, and we also observe that the bottom $k-1$ rows can all be written as linear combinations of the top $k$ rows, and therefore, the top $k$ rows form a basis for this lattice. Thus, we can write it as -
$$L' = \begin{bmatrix}
    a & -1 & \dots & 0\\
    a^2 & 0 & \dots & 0\\
    \vdots & \vdots & \ddots & \vdots\\
    a^{k-1} & 0 & \dots & -1\\
    M & 0 & \dots & 0\\
    \end{bmatrix} \implies L' \begin{bmatrix}
    X_1\\
    X_2\\
    X_3\\
    \vdots\\
    X_k\\
\end{bmatrix} = \begin{bmatrix}
    b\\
    b\frac{a^2-1}{a-1}\\
    \vdots\\
    b\frac{a^{k-1}-1}{a-1}\\
    0\\
\end{bmatrix} (\% M)$$
Also since, $X_i = 2^t * x_i + y_i$, the above equation can be re-written as -
$$L' \begin{bmatrix}
    y_1\\
    y_2\\
    y_3\\
    \vdots\\
    y_k\\
\end{bmatrix} = \left( b * \begin{bmatrix}
    1\\
    \frac{a^2-1}{a-1}\\
    \vdots\\
    \frac{a^{k-1}-1}{a-1}\\
    0\\
\end{bmatrix} + 2^t * L' * \begin{bmatrix}
    x_1\\
    x_2\\
    \vdots\\
    x_k\\
\end{bmatrix} \right) (\% M) = \text{(let) } c (\% M)$$

Consider the LLL reduced basis for $L'$ denoted by $L'_r$, and consider $c_r$ such that each element of $c_r$ is $\le \frac{M}{2}$ in absolute value(ensuring $c_r (\% M) = c (\% M)$). Then, the mentioned paper shows that there exists **atmost one integral "small" solution** to the (non-modular) linear equation $L'_r \cdot y = c_r$, where $y$ denotes the vector consisting of entries $y_1$ upto $y_k$! Thus, we can solve for $y$ by computing the inverse of $L'_r$. Thus, the obtained first coordinate of $y$ would be our $y_1$; and we can then obtain $X_1$ as $X_1 = 2^t * x_1 + y_1$. The only catch here is whether this "small" solution indeed is the correct solution, that is whether our expected $y$ indeed satisfies the mentioned norm bounds. The paper proves that for random LCGs this holds true with a good probability, given sufficient number of output-bits and sufficient information to be guessed.

## Our Work
We have implemented both the aforementioned attacks in python3. The attack on LCG allows us to recover the seed easily. However, the lattice-based attack on truncated LCG is somewhat different in the sense that the method only allows us to recover $X_1$, which is not the seed we seek.  
One solution is to modify our original algorithm to include $X_0$ in the unknown vector $y$; however since $X_0$ does not have a known $x_0$ part, this modification may actually yield results much worse, since $X_0$ is a possibly large vector, and hence the norm bounds may now be violated! We had tried this method earlier, but it couldn't correctly recover the seed in many cases, especially the cases in which $a$ was not co-prime with $M$!  
Another possible solution is to realize that $X_0$ in most cases needn't be unique, since the only outputs we obtain start from $X_1$! Thus, there could be multiple possible $X_0$'s which could yield the same sequence. We rely on our algorithm to obtain $X_1$, and then a SAT solver is incorporated to find out all possible $X_0$ which could yield the expected $X_1$. This way, we do not have to rely on the existence of the modular inverse of $a$, and several possible existing seeds can be recovered easily.

Another attack on truncated LCGs was implemented which doesn't rely on the knowledge or absence thereof of the parameters of the (truncated) LCG. This attack proceeds by modelling the parameter recovery problem for LCG as a [Satisfiability Modulo Theories](https://en.wikipedia.org/wiki/Satisfiability_modulo_theories) (SMT) decision problem.  
We again used the SMT solver [Z3Prover](https://github.com/Z3Prover/z3) and encoded our parameter recovery problem by constructing the unknown parameters as bit-vectors of same length as our modulus(that is, even though the parameters might be unknown, their maximum possible bit-lengths are assumed to be known and so is the truncation!).

### Observations
We observed that in cases where $a$ was not co-prime to $M$, say for example when $M$ is a power of $2$ and $a$ is even, multiple solutions always exist! Both the mentioned attacks on truncated LCGs were able to recover these multiple solutions.  
It was observed that Lattice-based attack is much more faster than SAT-based attack. However, since the lattice attack requires a bit more information than the information-theoretic lower bound, it doesn't perform very well in cases where the number of outputs are just enough in bit-sizes! In these cases, however, a SAT-based attack still performs very well, and as is usual, multiple solutions if possible, are reported. 

## Challenges
Finding a good library for lattice algorithms as well as which allows matrix operations to be performed on Integer or Rational Matrices was a challenge we faced.  
Another challenge was that the original paper on lattice-based attack didn't describe the algorithm in much detail, neither was the intuition of the algorithm obvious, so we spent a lot of time understanding the working of the same, and had also faced several issues, such as the $X_0$ issue mentioned above.  

## Limitations
Lattice-based attack relies on knowing the parameters of the truncated LCG. Had we not known those parameters, the above attack wouldn't work at all.  
All attacks rely on knowing the truncation as well as the bit-length of the modulus, none of the attacks are so general as to be able to figure these out on the fly!  
Another major limitation is that since SAT solvers work in exponential or subexponential time, they're much slower compared to the lattice attack which works in polynomial time.


## References
- [LCG](https://en.wikipedia.org/wiki/Linear_congruential_generator)
- [Truncated LCG Reconstruction](https://www.math.cmu.edu/~af1p/Texfiles/RECONTRUNC.pdf)  
Reconstructing Truncated Integer Variables Satisfying Linear Congruences  
Alan M. Frieze, Johan Hastad, Ravi Kannan, Jeffrey C. Lagarias, and Adi Shamir  
SIAM Journal on Computing 1988 17:2, 262-280 
- [Z3Prover](https://github.com/Z3Prover/z3)