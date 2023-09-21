import numpy as np


# Define the seed of the RNG
seed = 10
generator = np.random.default_rng(seed=seed)
    

def branch_or_not(p_branch):
    # IDEA: checks if branching event occurs or not 
    
    branching = False
    p_random = generator.random()
    if p_random < p_branch:
        branching = True
        
    return branching



def BBP_N(p_branch, max_time, time_steps, rate=1):
    # IDEA: returns population array and time array for a binary branching process

    N_array = [1]
    time_array = np.linspace(0, max_time, time_steps)
    N = 1
    current_time = 0
    index = 1

    while current_time <= max_time and N > 0:
        waiting_time = generator.exponential(1/(N*rate))
        if branch_or_not(p_branch) == True:
            N += 1
        else:
            N -= 1
            
        current_time += waiting_time
        
        # Only get population data for wanted time points
        for position, t in enumerate(time_array[index:]):
            if current_time > time_array[-1]: 
                N = 0 # Break while loop
                N_array = N_array + [N_array[-1]]*(len(time_array)-index) 
                break
            elif current_time < t and N == 0:
                if N == 0:
                    N_array.append(N)
                break
            elif current_time >= t and current_time < time_array[index+position+1]:
                N_array.append(N)
                index += position+1
                break
            elif current_time > t and current_time > time_array[index+position+1]:
                N_array.append(N_array[-1])
                
    return time_array, N_array


def BBP_S(p_branch, max_time, time_steps, rate=1):
    # IDEA: returns size array and time array for a binary branching process
    
    time_array = np.linspace(0, max_time, time_steps)
    N = 1
    S = 1
    S_array = [1]
    current_time = 0
    index = 1
    
    while current_time <= max_time and N > 0:
        waiting_time = generator.exponential(1/(N*rate))
        if branch_or_not(p_branch) == True:
            N += 1
            S += 2
        else:
            N -= 1
            
        current_time += waiting_time
        
        # Only get population data for wanted time points
        for position, t in enumerate(time_array[index:]):
            if current_time > time_array[-1]: 
                N = 0 # Break while loop
                S_array = S_array + [S_array[-1]]*(len(time_array)-index) 
                break
            elif current_time < t and N == 0:
                if N == 0:
                    S_array.append(S_array[-1])
                break
            elif current_time >= t and current_time < time_array[index+position+1]:
                S_array.append(S)
                index += position+1
                break
            elif current_time > t and current_time > time_array[index+position+1]:
                S_array.append(S_array[-1])
            
    return time_array, S_array
