import numpy as np
from Mutation import *

def branch_or_not(p_branch):
    branching = False
    p_random = np.random.rand()
    if p_random < p_branch:
        branching = True
    return branching




def BBP_N(p_branch, max_time, rate=1):
    generator = np.random.default_rng()
    N_array = [1]
    time_array = [0]
    N = 1
    current_time = 0
    
    while current_time <= max_time and N > 0:
        waiting_time = generator.exponential(1/(N*rate))
        if branch_or_not(p_branch) == True:
            N += 1
        else:
            N -= 1
        current_time += waiting_time
        N_array.append(N)
        time_array.append(current_time)

    return time_array, N_array



def BBP_S(p_branch, max_time, rate=1):
    generator = np.random.default_rng()
    time_array = [0]
    N = 1
    S = 1
    S_array = [1]
    current_time = 0
    while current_time <= max_time and N > 0:
        waiting_time = generator.exponential(1/(N*rate))
        if branch_or_not(p_branch) == True:
            N += 1
            S += 2
        else:
            N -= 1
            
        current_time += waiting_time
        time_array.append(current_time)
        S_array.append(S)
        
    return time_array, S_array



#cell_ID festlegen sepereat festlegen? 
#bessere datentruktur als dicitonaries? 

#wie viele sequenzen insgesamt mutiert (und ausgestorben)
#was ist die max-distance der mutierten Sequenzen? 
#wie ist die verteilung - wie viel nur eine mutation wie viel wie viele andere mutationen?
#Sequenz finden und schauen, ob die schon mutiert ist 
#Most recent common ancestor finden 


#!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Habe ich das jetzt drinnen, das auch wirklich die Sequnez vom Vater genommen wird?? 
#diese Mutationsrate biologisch sinnvoll und wauf was anwendbar?
#!!!!!!!!!!!!!!!!!!!!!

def BBP_DNA(p_branch, max_time, p_mut, DNA_seq, rate=1):
    #Each different DNA sequence has its own ID
    
    generator = np.random.default_rng()
    N_array = [1]
    seq_count_array = [1] # Stores the number if different sequences per time
    time_array = [0]
    N = 1
    ID_seq = {0:DNA_seq} # Contains IDs and sequences of all DNA sequences which have ever existed
    ID_count_current_alive = {0:1} # Contains IDs and number of those sequences which are currently alive
    current_time = 0
    ID_count = 0
    
    while current_time <= max_time and N > 0:
        waiting_time = generator.exponential(1/(N*rate))
        ID_parent = generator.choice(list(ID_count_current_alive)) # Randomly choose ID of one living DNA sequence
        
        if branch_or_not(p_branch) == True:
            N += 1
            p_random_mut =  np.random.rand()
            
            if p_random_mut < p_mut: 
                # The sequence of the parent gets passed on and an additional mutant sequence is added 
                new_seq = mutate(ID_seq[ID_parent])
                
                if new_seq in ID_seq.values(): # Check if mutated sequence has already existed (i.e. has same mutation or is a back mutation)
                    for key, value in ID_seq.items(): # Get old ID
                        if value == new_seq:
                            ID_new = key
                            
                            if ID_new in ID_count_current_alive: # Check if ID is currently alive
                                ID_count_current_alive[ID_new] = ID_count_current_alive[ID_new] + 1
                            else:
                                ID_count += 1
                                ID_new = ID_count
                                ID_seq[ID_new] = new_seq
                                ID_count_current_alive[ID_new] = 1
                            break
                else:
                    ID_count += 1
                    ID_new = ID_count
                    ID_seq[ID_new] = new_seq
                    ID_count_current_alive[ID_new] = 1
                    
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
        N_array.append(N)
        time_array.append(current_time)
        seq_count_array.append(len(ID_count_current_alive))
        
    return time_array, seq_count_array #, N_array



"""import matplotlib.pyplot as plt
p_branch = 0.8
max_time = 20
p_mut = 0.003
DNA_seq = generate_DNA_seq(100)
rate = 1


time, seq, N = BBP_DNA(p_branch, max_time, p_mut, DNA_seq, rate)
import matplotlib.pyplot as plt
plt.plot(time, np.log(seq))
plt.plot(time, np.log(N))

print(seq[-1])"""
















