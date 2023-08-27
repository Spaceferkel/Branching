
import random

    
def generate_DNA_seq(DNA_length):
    nucs = ["A", "T", "C", "G"]
    seq = ""
    
    for i in range(0, DNA_length):
        seq += random.choice(nucs)

    return seq



def mutate(seq):
    #mutates one base in the DNA sequence
    nucs = ["A", "T", "C", "G"]
    
    position = random.randint(0, len(seq)-1)
    old_nuc = seq[position]
    nucs.remove(old_nuc)
    new_nuc = random.choice(nucs)
    
    seq_mutated = seq[:position] + new_nuc + seq[position+1:]
    return seq_mutated
        


def count_diff_seq(DNA_dic):
    count_per_gen = {}
    for i in DNA_dic:
        diff_seq = len(set(DNA_dic[i]))
        count_per_gen[i] = diff_seq
        
    return count_per_gen
 
