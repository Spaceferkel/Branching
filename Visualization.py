import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd
#from Binary_BP import BP
from Binary_BP_mut import BP


nsims = 50
generations = 30
p_values = [0.3, 0.5, 0.55, 0.8]


def multiplot_nsims(nsims, generations, p_values):
    
    #Set plot parameters
    sns.set_style("white")
    plt.figure(figsize=(12, 8))
    
    #Get BP data for every wanted probabilitiy
    data_list = []  # List to store data
    for p_branch in p_values:
        temp_list = []
        for j in range(nsims):
            data = BP(p_branch, generations)[0]  
            temp_list.append(data)
        data_list.append(temp_list)
      
    #Generate subplots
    for idx, p_branch in enumerate(p_values, 1):
        
        plt.subplot(2, 2, idx)
        sns.lineplot(data=data_list[idx - 1], dashes=False, palette="rocket_r")
        plt.xlabel('Generation', fontsize=14)
        plt.ylabel('Population', fontsize=14)
        plt.title(f'p_branch = {p_branch}', fontsize=16)
        plt.legend().set_visible(False)
        sns.despine(right=True, top=True)
        
        for x in np.arange(0, generations + 1, step=5):
            plt.axvline(x=x, color='gray', linestyle='-', alpha=0.3)
    
    # Adjust layout and show the plot
    plt.suptitle(f"Binary Branching Process (Number of simulations = {nsims})", fontsize=18)
    plt.tight_layout(rect=[0, 0.07, 1, 0.95])
    plt.savefig(f"BP_multiplot_{nsims}_sims_size.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0





p_values = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
nsims = 1000
generations = 20
p_mut = 0.1
DNA_length = 1000


def plot_avg(p_values, generations, nsims, p_mut, DNA_length, var="N"):
    #IDEA: plot average Population or Size vs generations (var=N/S)
    
    sns.set_style("whitegrid")
    sns.set_palette("rocket")
    plt.figure(figsize=(20, 8))  
    
    # Check wheter to plot Population or Size
    if var == "N":
        i = 0
    elif var == "S":
        i = 1
    elif var =="Number_diff_seq":
        i = 2
    else: 
        raise ValueError("Not permitted variable")
        
    for p_branch in p_values:
        
        # Perform nsims sumulations and generate table with data
        plots_df = pd.DataFrame({'Generation': np.arange(generations + 1)})  
        for j in range(nsims):
            data = BP(p_branch, generations, p_mut, DNA_length)[i]
            df = pd.DataFrame({f'Population_{j}': data})  
            plots_df = pd.concat([plots_df, df], axis=1)  # Each simulation's data as a new column
        
        
        mean_array = np.array(plots_df.iloc[:, 1:].mean(axis=1))
        
        if var == "N":
            if p_branch < 0.5:
                zero_indices = np.where(mean_array==0)[0]
                if zero_indices.size > 0:
                    mean_gen_death = zero_indices[0] - 1
                mean_array[mean_array==0]=np.nan
                plt.plot(range(0, generations+1), mean_array, label=f"p_branch={p_branch}, mean death={mean_gen_death}")
            else:
                plt.plot(range(0, generations+1), mean_array, label=f"p_branch={p_branch}")
        elif var == "Number_diff_seq":
           plt.plot(range(0, generations+1), mean_array, label=f"p_branch={p_branch}") 
            
        else:
            if p_branch == 0.5:
                slope = round(linregress(range(0, generations+1), mean_array)[0], 3)
                plt.plot(range(0, generations+1), mean_array, label=f"p_branch={p_branch}, slope={slope}")
            else:
                plt.plot(range(0, generations+1), mean_array, label=f"p_branch={p_branch}")
  
    #Plot
    plt.xlabel('Generation', fontsize=16)
    plt.ylabel(f'<{var}>', fontsize=16)
    plt.ylim(-1, 10)
    plt.xlim(0, 15)
    plt.title(f'<{var}(g)> for {nsims} simulations', fontsize=16)
    plt.xticks(np.arange(0, generations + 1, step=5))
    plt.legend(fontsize=18)
    plt.savefig(f"BP.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0
