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
