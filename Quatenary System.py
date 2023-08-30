# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 22:08:13 2023

@author: jakob
"""

import numpy as np


# DNA als quatenary schreiben und nicht als string

generator = np.random.default_rng()

def generate_ID(seq):
    ID = 0
    for index, number in enumerate(reversed(seq)):
        ID += number*4**index
    
    return ID


def branch_or_not(p_branch):
    branching = False
    p_random = generator.rand()
    if p_random < p_branch:
        branching = True
    return branching


def mutate(seq):
    #mutates one base in the DNA sequence
    
    nucs = [0, 1, 2, 3]
    position = generator.integers(0, len(seq))
    old_nuc = seq[position]
    nucs.remove(old_nuc)
    new_nuc = generator.choice(nucs)
    seq[position] = new_nuc
    
    return seq
        

def generate_DNA_seq(DNA_length):
    nucs = [0, 1, 2, 3]
    seq = []
    
    for i in range(0, DNA_length):
        seq.append(generator.choice(nucs))

    return seq
