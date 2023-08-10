import random
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def BP(p_branch, generations):
    #IDEA: start with one element in 0th generation. Each element branches with p_branch in 2 or zero offspring
    #NOTE: When an element branches, it is also it's own offspring (does not increase size)
    
    N = 1 #one element in 0th generation
    N_dic = {0:1} #initialize N and S (total size) for 0th generation; Store N(g), S(g) here
    S_dic = {0:1}
    g = 1   #start with second generation
    
    while g <= generations: #go through each generation
    
        N_next = 0  #reset number of elements in generation every new generation (can't use N because that would inflict with the 1st generation)
        S = 0
        for element in range(1, N+1): #go through each element in generation
            p_e = random.random()
            if p_e <= p_branch:
                offspring = 2 
                S +=1
                p_branching_occurred = True
            else: 
                offspring = 0
            N_next += offspring
            
        S_dic[g] = S_dic[g-1] + S  
        N = N_next
        N_dic[g] = N
        g += 1
    
    Ns = np.array(list(N_dic.values()), dtype=np.double)
    #Ns[Ns==0]=np.nan #set 0-values to nan so that they won't be plotted (if I do not want to plot a dead process any further)
    S_array = np.array(list(S_dic.values()), dtype=np.int)
    
    return (Ns, S_array)

