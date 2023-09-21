from BBP_N_new import BBP_N, BBP_S
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from DNA_new import *
from scipy.stats import linregress



def nsims_mean_N(p_values, max_time, time_steps, nsims, rate):
    # IDEA: simulate the BBP nsims times and plot the mean of the population over time for every branching probability
    # NOTE: the larger the time_steps, the more accurate the result
    
    #sns.set_style("whitegrid")
    plt.style.use('classic')
    plt.figure(figsize=(20, 8))  
    plt.grid(True)
    
    time = np.linspace(0, max_time, time_steps)
    all_data = {'Time': time}

    for p_branch in p_values:
        data = {'Time': time}
        
        for j in range(nsims):
            time_array, population_data = BBP_N(p_branch, max_time, time_steps, rate)
            
            # Add zeros to dead processes so that the calculation of the mean is correct
            if len(population_data) < time_steps:
                population_data = population_data + [0]*(time_steps - len(population_data))
            
            # Store data of each simulation run
            data[f'Population_{j}'] = population_data
            print(f"simulation {j} done")
            
        df = pd.DataFrame.from_dict(data)
        mean_array = np.array(df.iloc[:, 1:].mean(axis=1))
        
        # Do not plot dead processes further
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
            
        all_data[f"p_branch {p_branch}"] = mean_array
        print(f"branching probability {p_branch} done")
    
    # Save data in a file
    p_values_string = "\t\t".join(["p: " + str(i) for i in p_values])
    column_names = "time\t\t" + p_values_string
    df_all_data = pd.DataFrame.from_dict(all_data)
    header = f"MEAN POPULATION <N>(t)\n\nNumber of simulations: {nsims}\nTime: {max_time}\nRNG seed: {seed}\
        \n\n{column_names}"
    np.savetxt(f"Mean_N_{nsims}_simulations.txt", df_all_data, header = header, delimiter="\t", fmt="%e")
            
    #Plot
    plt.xlabel('time', fontsize=16)
    plt.ylabel('<N>', fontsize=16)
    #plt.ylim(0, 100)
    #plt.xlim(0, 25)
    plt.title(f'<N>(t) for {nsims} simulations', fontsize=16)
    plt.legend(fontsize=18)
    #plt.savefig("33.png", dpi=300, bbox_inches='tight')
    plt.show()
  
    return 0



def nsims_mean_size(p_values, max_time, time_steps, nsims, rate):
    # IDEA: simulate the BBP nsims times and plot the mean of the size over time for every branching probability
    # NOTE: the larger the time_steps, the more accurate the result
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(20, 8))  

    time = np.linspace(0, max_time, time_steps)
    all_data = {'Time': time}

    for p_branch in p_values:
        data = {'Time': time}

        for j in range(nsims):
            time_data, size_data = BBP_S(p_branch, max_time, time_steps, rate)
            
            # Add last element of size to a dead process for the rest of the time     
            if len(size_data) < time_steps:
                size_data = size_data + [size_data[-1]]*(time_steps - len(size_data))
            
            # Add the data for every simulation to a DataFrame
            data[f'Size_{j}'] = size_data
            print(f"simulation {j} done")
            
        df = pd.DataFrame.from_dict(data)
        mean_array = np.array(df.iloc[:, 1:].mean(axis=1))
        all_data[f"p_branch {p_branch}"] = mean_array

        if p_branch == 0.5:
            slope, intercept, r_value, p_value, std_err = linregress(time, mean_array)
            plt.plot(time, mean_array, label=f"p_branch={p_branch}, slope: {round(slope, 2)}") 
        else:
            plt.plot(time, mean_array, label=f"p_branch={p_branch}")                

        print(f"branching probability {p_branch} done")
            
    # Save data in a file
    p_values_string = "\t\t".join(["p: " + str(i) for i in p_values])
    column_names = "time\t\t" + p_values_string
    df_all_data = pd.DataFrame.from_dict(all_data)
    header = f"MEAN SIZE <S>(t)\n\nNumber of simulations: {nsims}\nTime: {max_time}\nRNG seed: {seed}\
        \n\n{column_names}"
    np.savetxt(f"Mean_S_{nsims}_simulations.txt", df_all_data, header = header, delimiter="\t", fmt="%e")
    
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
    
    sns.set_style("white")
    plt.figure(figsize=(20, 8))  
    
    time = np.linspace(0, max_time, time_steps)
    all_data = {'Time': time}
    
    for p_branch in p_values:
        N_data = {'Time': time} 
        Seq_data = {'Time': time}
        
        for j in range(nsims):
            time_data, population_data, sequences = BBP_DNA(p_branch, max_time, time_steps, p_mut, DNA_seq, rate)
            
            # Add zeros to dead processes so that the calculation of the mean is correct
            if len(population_data) < time_steps:
                population_data = population_data + [0]*(time_steps - len(population_data))
                sequences = sequences + [0]*(time_steps - len(sequences))
            
            # Add the data for every simulation to a DataFrame
            N_data[f'Population_{j}'] =  population_data
            Seq_data[f'Number_different_seq_{j}'] = sequences
            print(f"simulation {j} done")
        
        df_N = pd.DataFrame.from_dict(N_data)
        df_Seq = pd.DataFrame.from_dict(Seq_data)
        mean_N_array = np.array(df_N.iloc[:, 1:].mean(axis=1))
        mean_seq_array = np.array(df_Seq.iloc[:, 1:].mean(axis=1))
        
        # Do not plot dead processes further
        if p_branch < 0.5:
            zero_indices = np.where(mean_N_array==0)[0]
            if zero_indices.size > 0:
                index_mean_time_death = zero_indices[0]
                mean_time_death = time_data[index_mean_time_death]
                mean_N_array[mean_N_array==0]=np.nan
                mean_N_array[index_mean_time_death] = 0
                mean_seq_array[index_mean_time_death] = 0
                mean_seq_array[index_mean_time_death+1:] = np.nan
                #plt.plot(time, mean_N_array, label=f"<N>(t), p_branch={p_branch}, mean time of death: {round(mean_time_death, 2)}") 
                plt.plot(time, mean_seq_array, label=f"<Seq>(t), p_branch={p_branch}") 
            else:
                #plt.plot(time, mean_N_array, label=f"<N>(t), p_branch={p_branch}") 
                plt.plot(time, mean_seq_array, label=f"<Seq>(t), p_branch={p_branch}") 
        else:
            #plt.plot(time, mean_N_array, label=f"<N>(t), p_branch={p_branch}") 
            plt.plot(time, mean_seq_array, label=f"<Seq>(t), p_branch={p_branch}")
            
        all_data[f"p_branch {p_branch}"] = mean_seq_array
        print(f"{p_branch} done")
        
    # Save data in a file
    p_values_string = "\t\t".join(["p: " + str(i) for i in p_values])
    column_names = "time\t\t" + p_values_string
    df_all_data = pd.DataFrame.from_dict(all_data)
    header = f"MEAN NUMBER OF DIFFERENT DNA SEQUENCES <Seq>(t)\n\nNumber of simulations: {nsims}\nTime: {max_time}\nRNG seed: {seed}\
        \n\n{column_names}"
    np.savetxt(f"Mean_Seq_{nsims}_simulations.txt", df_all_data, header = header, delimiter="\t", fmt="%e")
    
    #Plot
    plt.xlabel('time', fontsize=16)
    plt.ylabel('<Seq>', fontsize=16)
    #plt.ylim(0, 8)
    #plt.xlim(0, 25)
    main_title = "Average number of different DNA sequences over time <Seq>(t)"
    sub_title = f"{nsims} simulations, p_mut={p_mut}"    
    plt.suptitle(main_title, fontsize=22, color='black', y=1.05)
    plt.text(0.5, 0.98, sub_title, fontsize=18, color='#333333', ha='center', va='center', transform=plt.gcf().transFigure)
    plt.legend(fontsize=18)
    plt.savefig("33.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    return 0
