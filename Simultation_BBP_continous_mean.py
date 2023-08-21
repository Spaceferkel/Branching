from BBP_continous_optimized import BBP_continous
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


p_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9]
nsims = 1000
max_time = 100
rate = 1
time_steps = 1000


def nsims_mean(p_values, max_time, time_steps, nsims, rate):
    #IDEA: simulate the BBP nsims times and plot mean for every branching probability
    #NOTE: the larger the time_steps, the more accurate is the result
    
    sns.set_style("whitegrid")
    sns.set_palette("rocket")
    plt.figure(figsize=(20, 8))  
        
    time = np.linspace(0, max_time, time_steps)
    
    for p_branch in p_values:
        plots_df = pd.DataFrame({'Time': time})  
        for j in range(nsims):
            N_array = [1]
            time_data, population_data = BBP_continous(p_branch, max_time, rate)
            
            #get N for same times (Originally, time_data for different BBP_continous calls is different)
            for t in time:
                index = next((i for i, td in enumerate(time_data) if td >= t), None)
                if index:
                    N_array.append(population_data[index - 1])
            
            #Add zeros to dead processes so the calulation of the mean is correct        
            if len(N_array) < time_steps:
                N_array = N_array + [0]*(time_steps - len(N_array))
                
            #Add the data for every simulation to a DataFrame
            df = pd.DataFrame({f'Population_{j}': N_array}) 
            plots_df = pd.concat([plots_df, df], axis=1)
        
        mean_array = np.array(plots_df.iloc[:, 1:].mean(axis=1))
        
        #Do not plot dead processes further
        if p_branch < 0.5:
            zero_indices = np.where(mean_array==0)[0]
            if zero_indices.size > 0:
                index_mean_time_death = zero_indices[0]
                mean_time_death = time[index_mean_time_death]
                mean_array[mean_array==0]=np.nan
                mean_array[index_mean_time_death] = 0
                plots_df = pd.concat([plots_df, pd.DataFrame(mean_array)], axis=1)
                plt.plot(time, mean_array, label=f"p_branch={p_branch}, mean time of death: {round(mean_time_death, 2)}")       
            else:
                plots_df = pd.concat([plots_df, pd.DataFrame(mean_array)], axis=1)
                plt.plot(time, mean_array, label=f"p_branch={p_branch}") 
        else:
            plots_df = pd.concat([plots_df, pd.DataFrame(mean_array)], axis=1)
            plt.plot(time, mean_array, label=f"p_branch={p_branch}")                
            
    #Plot
    plt.xlabel('time', fontsize=16)
    plt.ylabel('<N>', fontsize=16)
    plt.ylim(-1, 10)
    plt.xlim(0, 42)
    plt.title(f'<N(t)> for {nsims} simulations', fontsize=16)
    plt.legend(fontsize=18)
    plt.savefig("BBP_continous_mean.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0

nsims_mean(p_values, max_time, time_steps, nsims, rate)

