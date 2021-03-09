from collections import Counter
from itertools import combinations
from math import log2
from z3 import *
from tqdm import tqdm
from fractions import Fraction

def parity(x):
    res = 0
    while x:
        res^=1
        x=x&(x-1)
    return res

def calc_bias(sbox,no_sign=True,fraction=False):
    n = len(sbox)
    bias = Counter({(i,j):-(n//2) for i in range(n) for j in range(n)})
    for imask in tqdm(range(n),desc='calculating sbox bias'):
        for omask in range(n):
            for i in range(n):
                bias[(imask,omask)]+= parity((sbox[i]&omask)^(i&imask))^1
    if no_sign:
        for i in bias:
            bias[i]=abs(bias[i])
    if fraction:
        for i in bias:
            bias[i]=Fraction(bias[i]/n)
    return bias

def print_bitrelations(inp_masks,out_masks,n,s):
    """
    s = num bits in sbox
    n = block size
    """
    def bin_sep(val):
        v = bin(val)[2:].zfill(n)
        return "|".join(v[i:i+s] for i in range(0,n,s))

    rounds = len(out_masks)
    for i in range(rounds):
        imask,omask = inp_masks[i],out_masks[i]
        print(bin_sep(imask))
        print(' '.join(['-'*s]*(n//s)))
        print(bin_sep(omask))
        print()
    print(bin_sep(inp_masks[-1]))

def gen_pbox(s,n):
    '''n s-bit sboxes'''
    return [ (s*i+j)%(n*s) for j in range(s) for i in range(n) ]

def get_optimal_masks(sbox,pbox,num_rounds,bias=None,non_zero=[0],prune_level=0):
    n = int(log2(len(sbox)))
    num_blocks = len(pbox)//n
    if not bias:
        bias = calc_bias(sbox)
    sboxf = Function('sbox',BitVecSort(n),BitVecSort(n),RealSort())
    def permutation(inp,oup,pbox):
        pn = len(pbox)
        constraints = []
        for i,v in enumerate(pbox):
            constraints.append(
                Extract(pn-1-i,pn-1-i,inp)==Extract(pn-1-v,pn-1-v,oup)
            )
        return constraints
    constraints = []
    for i in range(2**n):
        for j in range(2**n):
            # just some pruning of very small biases
            if bias[(i,j)]>=2**(prune_level):
                constraints.append(sboxf(i,j)==bias[(i,j)])
            else:
                constraints.append(sboxf(i,j)==0)
    s = Optimize()
    inps = [[BitVec('r{}_i{}'.format(r,i),n) for i in range(num_blocks)] 
            for r in range(num_rounds+1)]
    oups = [[BitVec('r{}_o{}'.format(r,i),n) for i in range(num_blocks)] 
            for r in range(num_rounds)]
    objectives = [
        # the actual objective, which is just product of bias [0,1/2]
        2**(num_blocks*num_rounds-1)*Product([
            sboxf(
                inps[i//num_blocks][i%num_blocks],
                oups[i//num_blocks][i%num_blocks]) 
            for i in range(num_blocks*num_rounds)
        ]),        
        # reducing optimization problem of product to sums 
        Sum([
            sboxf(
                inps[i//num_blocks][i%num_blocks],
                oups[i//num_blocks][i%num_blocks]) 
            for i in range(num_blocks*num_rounds)
        ]),
        # objective when the input biases are [0,2**n] just the final 
        # division
        2**(num_blocks*num_rounds-1)*Product([
            sboxf(
                inps[i//num_blocks][i%num_blocks],
                oups[i//num_blocks][i%num_blocks]) 
            for i in range(num_blocks*num_rounds)
        ])/((2**n)**(num_blocks*num_rounds))
    ]
    #for objective in objectives:
    s.maximize(objectives[1])
    s.add(constraints)
    s.add(Not(And( *[inps[0][i]==0 for i in range(num_blocks)])))
    # the last layer is input, which we would like to be
    # reduced to as few sboxes as possible
    if non_zero: #if specified which boxes to look for
        for i in range(num_blocks):
            if i in non_zero:
                s.add(inps[-1][i]!=0)
            else:
                s.add(inps[-1][i]==0)
    #s.add(PbEq([(i!=0,1) for i in inps[-1]],1))dd
    for r in range(num_rounds):
        for i in range(num_blocks):
            # if sbox has input, it should have ouput
            s.add(Implies(inps[r][i]!=0,oups[r][i]!=0))
            # if sbox has no input it should not have any output
            s.add(Implies(inps[r][i]==0,oups[r][i]==0))
            # skip taking input/outputs with no bias
            s.add(
                Implies(
                    And(inps[r][i]!=0,oups[r][i]!=0),
                    sboxf(inps[r][i],oups[r][i])!=0
                )
            )
    # permutation of output of sboxes are inputs of next round
    for i in range(num_rounds):
        s.add(permutation(Concat(oups[i]),Concat(inps[i+1]),pbox))
    results = []
    #print("began searching")
    if s.check()==sat:
        m = s.model()
        inp_masks = [ m.eval( Concat(inps[i])).as_long() 
                     for i in range(num_rounds+1) ]
        oup_masks = [ m.eval( Concat(oups[i])).as_long() 
                         for i in range(num_rounds) ]
        print_bitrelations(inp_masks,oup_masks,len(pbox),n)
        total_bias = m.eval(objectives[2]).as_fraction()
        print("total bias:",total_bias)
        return inp_masks,oup_masks,total_bias

if __name__=="__main__":
    num_rounds = int(input())
    pbox_size = int(input())
    pbox = list(map(int,input().strip().split()))
    sbox_size = int(input())
    sbox = list(map(int,input().strip().split()))
    x = get_optimal_masks(sbox,pbox,num_rounds-1,None,non_zero=[])
    print('Bit Masks: ',x[0],x[1])
