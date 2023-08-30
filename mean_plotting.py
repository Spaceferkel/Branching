from BBP_N_new import BBP_N, BBP_S
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from DNA_new import *
from scipy.stats import linregress


p_values = [0.1]
nsims = 1000
max_time = 20
rate = 1
time_steps = 100
p_mut = 0.003
DNA_seq = generate_DNA_seq(100)


def nsims_mean_N(p_values, max_time, time_steps, nsims, rate):
    # IDEA: simulate the BBP nsims times and plot the mean of the population over time for every branching probability
    # NOTE: the larger the time_steps, the more accurate the result
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(20, 8))  

    time = np.linspace(0, max_time, time_steps)

    for p_branch in p_values:
        plots_df = pd.DataFrame({'Time': time})  
        
        for j in range(nsims):
            time_data, population_data = BBP_N(p_branch, max_time, time_steps, rate)
            
            #Add zeros to dead processes so that the calculation of the mean is correct
            if len(population_data) < time_steps:
                population_data = population_data + [0]*(time_steps - len(population_data))
            
            #Add the data for every simulation to a DataFrame
            df = pd.DataFrame({f'Population_{j}': population_data}) 
            plots_df = pd.concat([plots_df, df], axis=1)
            print(f"simulation {j} done")
            
        mean_array = np.array(plots_df.iloc[:, 1:].mean(axis=1))

        #Do not plot dead processes further
        if p_branch < 0.5:
            zero_indices = np.where(mean_array==0)[0]
            if zero_indices.size > 0:
                index_mean_time_death = zero_indices[0]
                mean_time_death = time[index_mean_time_death]
                mean_array[mean_array==0]=np.nan
                mean_array[index_mean_time_death] = 0
                plt.plot(time, mean_array, label=f"p_branch={p_branch}, mean time of death: {round(mean_time_death, 2)}")       
            else:
                plt.plot(time, mean_array, label=f"p_branch={p_branch}") 
        else:
            plt.plot(time, mean_array, label=f"p_branch={p_branch}") 
            
        plots_df = pd.concat([plots_df, pd.DataFrame(mean_array)], axis=1)
        print(f"branching probability {p_branch} done")
            
    #Plot
    plt.xlabel('time', fontsize=16)
    plt.ylabel('<N>', fontsize=16)
    #plt.ylim(0, 100)
    #plt.xlim(0, 25)
    plt.title(f'<N>(t) for {nsims} simulations', fontsize=16)
    plt.legend(fontsize=18)
    #plt.savefig("BBP_continous_mean.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0



def nsims_mean_size(p_values, max_time, time_steps, nsims, rate):
    # IDEA: simulate the BBP nsims times and plot the mean of the size over time for every branching probability
    # NOTE: the larger the time_steps, the more accurate the result
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(20, 8))  

    time = np.linspace(0, max_time, time_steps)

    for p_branch in p_values:
        plots_df = pd.DataFrame({'Time': time})  
        
        for j in range(nsims):
            time_data, size_data = BBP_S(p_branch, max_time, time_steps, rate)
            
            #Add last element of size to a dead process for the rest of the time     
            if len(size_data) < time_steps:
                size_data = size_data + [size_data[-1]]*(time_steps - len(size_data))
            
            #Add the data for every simulation to a DataFrame
            df = pd.DataFrame({f'Population_{j}': size_data}) 
            plots_df = pd.concat([plots_df, df], axis=1)
            print(f"simulation {j} done")
            
        mean_array = np.array(plots_df.iloc[:, 1:].mean(axis=1))
        plots_df = pd.concat([plots_df, pd.DataFrame(mean_array)], axis=1)

        if p_branch == 0.5:
            slope, intercept, r_value, p_value, std_err = linregress(time, mean_array)
            plt.plot(time, mean_array, label=f"p_branch={p_branch}, slope: {round(slope, 2)}") 
        else:
            plt.plot(time, mean_array, label=f"p_branch={p_branch}")                

        print(f"branching probability {p_branch} done")
            
    #Plot
    plt.xlabel('time', fontsize=16)
    plt.ylabel('<S>', fontsize=16)
    #plt.ylim(0, 100)
    #plt.xlim(0, 25)
    plt.title(f"Average Size over time <S>(t) for {nsims} simulations", fontsize=16)
    plt.legend(fontsize=18)
    #plt.savefig("BBP_continous_mean.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0




def nsims_mean_DNA(p_values, max_time, time_steps, nsims, p_mut, DNA_seq, rate):
    # IDEA: simulate the BBP nsims times and plot the mean of the number of different DNA sequences over time for every branching probability
    # NOTE: the larger the time_steps, the more accurate is the result
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(20, 8))  
    
    time = np.linspace(0, max_time, time_steps)
    
    for p_branch in p_values:
        df_N_all = pd.DataFrame({'Time': time})  
        df_seq_all = pd.DataFrame({'Time': time})  
        
        for j in range(nsims):
            time_data, population_data, sequences = BBP_DNA(p_branch, max_time, time_steps, p_mut, DNA_seq, rate)
            
            #Add zeros to dead processes so that the calculation of the mean is correct
            if len(population_data) < time_steps:
                population_data = population_data + [0]*(time_steps - len(population_data))
                sequences = sequences + [0]*(time_steps - len(sequences))
            
            #Add the data for every simulation to a DataFrame
            df_N = pd.DataFrame({f'Population_{j}': population_data}) 
            df_N_all = pd.concat([df_N_all, df_N], axis=1)
            df_seq = pd.DataFrame({f'Number_different_seq_{j}': sequences})
            df_seq_all = pd.concat([df_seq_all, df_seq], axis=1)
            print(f"simulation {j} done")
        
        mean_N_array = np.array(df_N_all.iloc[:, 1:].mean(axis=1))
        mean_seq_array = np.array(df_seq_all.iloc[:, 1:].mean(axis=1))
        
        
        #Do not plot dead processes further
        if p_branch < 0.5:
            zero_indices = np.where(mean_N_array==0)[0]
            if zero_indices.size > 0:
                index_mean_time_death = zero_indices[0]
                mean_time_death = time_data[index_mean_time_death]
                mean_N_array[mean_N_array==0]=np.nan
                mean_N_array[index_mean_time_death] = 0
                mean_seq_array[index_mean_time_death] = 0
                mean_seq_array[index_mean_time_death+1:] = np.nan
                plt.plot(time, mean_N_array, label=f"<N>(t), p_branch={p_branch}, mean time of death: {round(mean_time_death, 2)}") 
                plt.plot(time, mean_seq_array, label=f"<S>(t), p_branch={p_branch}") 
            else:
                plt.plot(time, mean_N_array, label=f"<N>(t), p_branch={p_branch}") 
                plt.plot(time, mean_seq_array, label=f"<S>(t), p_branch={p_branch}") 

        else:
            plt.plot(time, mean_N_array, label=f"<N>(t), p_branch={p_branch}") 
            plt.plot(time, mean_seq_array, label=f"<S>(t), p_branch={p_branch}") 
        
        print(f"{p_branch} done")
        
    print(mean_N_array[-2], mean_seq_array[-2])
    #Plot
    plt.xlabel('time', fontsize=16)
    plt.ylabel('<N>, <S>', fontsize=16)
    plt.ylim(0, 2)
    #plt.xlim(0, 25)
    plt.title(f'Average number of Population and different DNA sequences over time for {nsims} simulations, p_mut={p_mut}', fontsize=16)
    plt.legend(fontsize=18)
    #plt.savefig("BBP_continous_mean.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0
