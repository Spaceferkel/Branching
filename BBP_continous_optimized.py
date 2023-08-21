import numpy as np
import heapq as hq


def branch_or_not(p_branch):
    #decides whether to branch or not
    branching = False
    p_random = np.random.rand()
    if p_random < p_branch:
        branching = True        
    return branching



def BBP_continous(p_branch, max_time, rate):
    #OPTIMIZATIONS: heap instead of list; rounded int instead of float
    
    max_time = max_time * 1000 #We calculte in integeter instead of float
    current_alive = []
    times = [0]
    population = [1]
    generation = 0
    gr = np.random.default_rng()
    life_time_ancestor = int(round(gr.exponential(1/rate), 3) * 1000)
    hq.heappush(current_alive, (life_time_ancestor, generation))
    current_birth_time = 0
    
    
    while current_birth_time <= max_time and len(current_alive) > 0:
        offspring = 0
        parent = hq.heappop(current_alive)
        current_generation = parent[1]
        current_birth_time = parent[0] #birth_time of offspring = death_time of parent
        if branch_or_not(p_branch) == True:
            offspring = 2
          
        for i in range(0, offspring):
            life_time = int(round(gr.exponential(1/rate),3) * 1000)
            hq.heappush(current_alive, (current_birth_time + life_time, current_generation+1))
        
        times.append(current_birth_time)
        population.append(len(current_alive))
        
    times_array = np.array(times) / 1000
    population_array = np.array(population)
    return (times_array, population_array)
