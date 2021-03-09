from collections import Counter
from itertools import combinations
from math import log2
from z3 import *
from tqdm import tqdm
from fractions import Fraction
import random

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


def get_optimal_masks(sbox,pbox,num_rounds,bias=None,non_zero=[0],prune_level=0):
    n = int(log2(len(sbox)))
    num_blocks = len(pbox)//n
    #print(n,num_blocks)
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

def get_all_pos_masks(sbox,pbox,num,num_key_blocks=1,prune=0):
    n = int(log2(len(sbox)))
    num_blocks = len(pbox)//n
    bias = calc_bias(sbox)
    round_masks = []
    try:
        for num_rounds in range(1,num+1):
            masks_this_round = []
            print("depth:",num_rounds)
            for pos in combinations(range(num_blocks),num_key_blocks):
                print("block positions:",pos)
                io_masks = get_optimal_masks(sbox,pbox,num_rounds,bias,pos,prune_level=prune)
                masks_this_round.append(io_masks)
            round_masks.append(masks_this_round)
    except KeyboardInterrupt:
        return round_masks
    return round_masks


def gen_pbox(s,n):
    '''n s-bit sboxes'''
    return [ (s*i+j)%(n*s) for j in range(s) for i in range(n) ]

def ROTL(value, bits, size): return \
            ((value % (1 << bits)) << (size - bits)) | (value >> bits)
class SPN:
    def __init__(self,SBOX,PBOX,key,rounds):
        self.SBOX = SBOX
        self.PBOX = PBOX
        self.SINV = [SBOX.index(i) for i in range(len(SBOX))]
        self.PINV = [PBOX.index(i) for i in range(len(PBOX))]
        self.BLOCK_SIZE = len(PBOX)
        self.BOX_SIZE = int(log2(len(SBOX)))
        self.NUM_SBOX = len(PBOX)//self.BOX_SIZE
        self.rounds = rounds
        self.round_keys = self.expand_key(key,rounds)

    def perm(self,inp:int) -> int:
        ct = 0
        for i,v in enumerate(self.PBOX):
            ct |= (inp>>(self.BLOCK_SIZE-1-i)&1)<<(self.BLOCK_SIZE-1-v)
        return ct

    def inv_perm(self,inp:int) -> int:
        ct = 0
        for i,v in enumerate(self.PINV):
            ct |= (inp>>(self.BLOCK_SIZE-1-i)&1)<<(self.BLOCK_SIZE-1-v)
        return ct

    def sbox(self,inp:int) -> int:
        ct,BS = 0,self.BOX_SIZE
        for i in range(self.NUM_SBOX):
            ct |= self.SBOX[(inp>>(i*BS))&((1<<BS)-1)]<<(BS*i)
        return ct

    def inv_sbox(self,inp:int) -> int:
        ct,BS = 0,self.BOX_SIZE
        for i in range(self.NUM_SBOX):
            ct |= self.SINV[(inp>>(i*BS))&((1<<BS)-1)]<<(BS*i)
        return ct
         
    def int_to_list(self,inp):
        BS = self.BOX_SIZE
        return [ (inp>>(i*BS))&((1<<BS)-1) 
                for i in range(self.NUM_SBOX-1,-1,-1) ]

    def list_to_int(self,lst):
        res = 0
        for i,v in enumerate(lst[::-1]):
            res |= v<<(i*self.BOX_SIZE)
        return res

    def expand_key(self,key,rounds):
        if isinstance(key,(list)):
            key = self.list_to_int(key)
        elif isinstance(key,(bytes,bytearray)):
            key = int.from_bytes(key,'big')
        key = key&((1<<self.BLOCK_SIZE)-1)
        keys = [key]
        for _ in range(rounds):
            key = ROTL(key,self.BOX_SIZE,self.BLOCK_SIZE)
            keys.append(key)
        return [random.randint(0,(1<<self.BLOCK_SIZE)-1) for _ in  range(rounds+1)]

    def enc(self,ct:int,round_keys=None) -> int:
        if round_keys is None:
            round_keys = self.round_keys
        for round_key in round_keys[:-2]:
            ct^= round_key
            ct = self.sbox(ct)
            ct = self.perm(ct)
        ct^= round_keys[-2]
        ct = self.sbox(ct)
        ct^= round_keys[-1]
        return ct
    
    def dec(self,ct,round_keys=None):
        if round_keys is None:
            round_keys = self.round_keys[::-1]
        # partial decryption
        # round keys in reverse order
        ct^= round_keys[0]
        ct = self.inv_sbox(ct)
        if len(round_keys)==1:
            return ct
        ct^= round_keys[1]
        for round_key in round_keys[2:]:
            ct = self.inv_perm(ct)
            ct = self.inv_sbox(ct)
            ct^= round_key
        return ct
 
    def parity_bias(self,inp_mask,out_mask,pt_ct_pairs,round_keys_rev):
        key_guesses = []
        key_blocks = [(1<<self.BOX_SIZE)-1 if i else 0 for i in self.int_to_list(out_mask)]
        num_key_blocks = sum(i!=0 for i in key_blocks)
        key_bits = self.list_to_int(key_blocks)
        if round_keys_rev: #not the last round
            key_bits = self.perm(key_bits)
        for i in tqdm(range((1<<self.BOX_SIZE)**num_key_blocks)):
            key = to_key(i,key_bits,self.BLOCK_SIZE)
            bias = Counter()
            for pt,ct in pt_ct_pairs:
                ct_partial = self.dec(ct,round_keys_rev + [key])
                bias[parity((pt&inp_mask)^(ct_partial&out_mask))]+=1
            key_guesses.append(bias)
        score,key_i = max((abs(key_guesses[i][1]-len(pt_ct_pairs)/2),i) 
            for i in range( (1<<self.BOX_SIZE)**num_key_blocks ))
        print(sorted((abs(key_guesses[i][1]-len(pt_ct_pairs)/2),i) 
            for i in range( (1<<self.BOX_SIZE)**num_key_blocks )))
        print("bias:",score/(len(pt_ct_pairs)))
        return to_key(key_i,key_bits,self.BLOCK_SIZE)
    
    def recover_round_keys(self,pt_ct_pairs,prune=0):
        all_round_masks = get_all_pos_masks(self.SBOX,self.PBOX,self.rounds-1,prune=prune)
        recovered_keys = []
        for round_masks in all_round_masks[::-1]: #current round depth
            key=0
            for mask in round_masks: #block positions in current round
                inp_mask,oup_mask = mask[0][0],mask[0][-1]
                key |= self.parity_bias(inp_mask,oup_mask,pt_ct_pairs,recovered_keys)
            recovered_keys.append(key)
        return recovered_keys

def to_key(i,key_bits,bs):
    key = 0
    for bit in range(bs):
        if (key_bits>>bit)&1:
            key |= (i&1)<<bit
            i>>=1
    return key

if __name__ == "__main__":
    SBOX = [ 8, 15, 2, 10, 12, 11, 6, 13, 14, 5, 4, 0, 1, 3, 9, 7 ]
    PBOX = [ 9, 8, 2, 6, 1, 7, 11, 3, 10, 4, 0, 5 ]
    s = SPN(SBOX,PBOX,0,2) #2 rounds
    pairs = []
    for i in range(8192):
        x = random.randint(0,2**16-1)
        pairs.append((x,s.enc(x)))
    recovered = s.recover_round_keys(pairs)
    print('actual last round key:',s.round_keys[-1],'{:018_b}'.format(s.round_keys[-1]))
    print('recovered last round key:',recovered[0],'{:018_b}'.format(recovered[0]))