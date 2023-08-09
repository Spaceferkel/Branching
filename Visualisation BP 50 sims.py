# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 17:36:49 2023

@author: jakob
"""

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from BP import BP

def plot_nsims_multiplots(nsims, generations, p_values):
    sns.set_style("white")
    population_data_list = []  # List to store population data
    
    for p_branch in p_values:
        temp_list = []
        for j in range(nsims):
            population_data = BP(p_branch, generations)  
            temp_list.append(population_data)
        population_data_list.append(temp_list)
    
    # Create a 2x2 grid of subplots
    plt.figure(figsize=(12, 8))
    
    for idx, p_branch in enumerate(p_values, 1):
        plt.subplot(2, 2, idx)
        sns.lineplot(data=population_data_list[idx - 1], dashes=False, palette="rocket_r")
        
        plt.xlabel('Generation', fontsize=12)
        plt.ylabel('Population', fontsize=12)
        plt.title(f'p_branch = {p_branch}', fontsize=14)
        
        plt.legend().set_visible(False)
        sns.despine(right=True, top=True)
        for x in np.arange(0, generations + 1, step=5):
            plt.axvline(x=x, color='gray', linestyle='-', alpha=0.3)
    
    # Adjust layout and save/show the plot
    plt.suptitle(f"Binary Branching Process (Number of simulations = {nsims})", fontsize=16)
    plt.tight_layout(rect=[0, 0.07, 1, 0.95])
    plt.savefig(f"BP_multiplot_{nsims}_sims.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0

# Call the function with desired parameters
nsims = 50
generations = 30
p_values = [0.3, 0.5, 0.55, 0.8]
plot_nsims_multiplots(nsims, generations, p_values)