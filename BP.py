# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 21:06:13 2023

@author: jakob
"""

import random
import numpy as np


def BP(p_branch, generations):
    #Idea: start with one element in 0th generation. Each element branches with p_branch in 2 or zero offspring
    
    N = 1 #one element in 0th generation
    N_dic = {0:1} #initialize N for 0th generation; Store N(g) here
    g = 1   #start with second generation
    S = N   #size of the colony
    
    while g <= generations: #go through generations 
        N_next = 0  #reset number of elements in generation every new generation (can't use N because that would inflict with the 1st generation)
        for element in range(1, N+1): #go through each element in generation
            p_e = random.random()
            if p_e <= p_branch:
                offspring = 2  
                S +=1
            else: 
                offspring = 0
            N_next += offspring
        
        N = N_next
        N_dic[g] = N
        g += 1
    
    Ns = np.array(list(N_dic.values()), dtype=np.double)
    Ns[Ns==0]=np.nan #set 0-values to nan so that they won't be plotted (don't want to plot a dead process any further)
    
    return Ns #return size S if wanted

