import numpy as np
import hashlib


# Set seed for RNG
seed = 444
generator = np.random.default_rng(seed=seed)


def branch_or_not(p_branch):
    # IDEA: checks if branching event occurs or not 
    
    branching = False
    p_random = generator.random()
    if p_random < p_branch:
        branching = True
    return branching


def generate_DNA_seq(DNA_length):
    # IDEA: generates a random DNA sequence with given length
    
    nucs = ["A", "T", "C", "G"]
    seq = ""
    for i in range(0, DNA_length):
        seq += generator.choice(nucs)

    return seq



def mutate(seq):
    # IDEA: mutates one base in the DNA sequence
    
    nucs = ["A", "T", "C", "G"]
    position = generator.integers(0, len(seq))
    old_nuc = seq[position]
    nucs.remove(old_nuc)
    new_nuc = generator.choice(nucs)
    
    seq_mutated = seq[:position] + new_nuc + seq[position+1:]
    return seq_mutated



def generate_ID(DNA_seq):
    # IDEA: generates a unique ID for the given DNA sequence
    
    sha256 = hashlib.sha256()
    sha256.update(DNA_seq.encode('utf-8'))
    ID = sha256.hexdigest()
        
    return ID
 

def BBP_DNA(p_branch, max_time, time_steps, p_mut, DNA_seq, rate=1):
    # IDEA: returns an array of the number of different DNA sequences, a population array and a time array for a binary branching process
    # NOTE: Each different DNA sequence has its own ID (and not every cell)
    
    generator = np.random.default_rng()
    N_array = [1]
    seq_count_array = [1] # Stores the number of different sequences per time
    time_array = np.linspace(0, max_time, time_steps)
    N = 1
    ID_ancestor = generate_ID(DNA_seq)
    ID_seq = {ID_ancestor:DNA_seq} # Contains IDs and sequences of all DNA sequences which have ever existed
    ID_count_current_alive = {ID_ancestor:1} # Contains IDs and number of those sequences which are currently alive
    current_time = 0
    index = 1
    
    while current_time <= max_time and N > 0:
        waiting_time = generator.exponential(1/(N*rate))
        counts_DNA_alive = np.array(list(ID_count_current_alive.values()))
        p_distribution = counts_DNA_alive/np.sum(counts_DNA_alive)
        ID_parent = generator.choice(list(ID_count_current_alive), p=p_distribution) # Choose ID of one living DNA sequence with given probability distribution

        if branch_or_not(p_branch) == True:
            N += 1
            p_random_mut = generator.random()
            
            if p_random_mut < p_mut: 
                # The sequence of the parent gets passed on and an additional mutant sequence is added 
                new_seq = mutate(ID_seq[ID_parent])
                ID_new = generate_ID(new_seq)
                
                if ID_new in ID_count_current_alive:
                    ID_count_current_alive[ID_new] = ID_count_current_alive[ID_new] + 1
                else:
                    ID_count_current_alive[ID_new] = 1
                    ID_seq[ID_new] = new_seq   
            else:
                ID_count_current_alive[ID_parent] = ID_count_current_alive[ID_parent] + 1 # If it doesn't mutate, one identical sequence is added
        else:
            # If it doesn't branch, DNA sequence dies
            N -= 1
            if ID_count_current_alive[ID_parent] > 1: #if multiple identical sequences exist, just lower their count
                ID_count_current_alive[ID_parent] = ID_count_current_alive[ID_parent] - 1
            else:
                del ID_count_current_alive[ID_parent]

        current_time += waiting_time
        
        # Only get data for wanted time points
        for position, t in enumerate(time_array[index:]):
            if current_time > time_array[-1]: 
                N = 0 # Break while loop
                N_array = N_array + [N_array[-1]]*(len(time_array)-index) 
                seq_count_array = seq_count_array + [seq_count_array[-1]]*(len(time_array)-index) 
                break
            elif current_time < t and N == 0:
                if N == 0:
                    N_array.append(N)
                    seq_count_array.append(len(ID_count_current_alive))
                break
            elif current_time >= t and current_time < time_array[index+position+1]:
                N_array.append(N)
                seq_count_array.append(len(ID_count_current_alive))
                index += position+1
                break
            elif current_time > t and current_time > time_array[index+position+1]:
                N_array.append(N_array[-1])
                seq_count_array.append(seq_count_array[-1])

    return time_array, N_array, seq_count_array
