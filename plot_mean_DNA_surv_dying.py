import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from DNA_new import *



def nsims_mean_DNA_surv_dying(p_branch, max_time, time_steps, nsims, p_mut, DNA_seq, rate):
    # IDEA: Plot Number of different DNA sequences for all, dying and surviving processes for a given branching probability
    # NOTE: the larger the time_steps, the more accurate is the result
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(20, 8))  
    
    time = np.linspace(0, max_time, time_steps)
    dying_processes = {'Time': time}
    surviving_processes = {'Time': time}
    all_processes = {'Time': time}
        
    for j in range(nsims):
        time_data, _, sequences = BBP_DNA(p_branch, max_time, time_steps, p_mut, DNA_seq, rate)
                
        # Store dead processes in a different DataFrame
        if len(sequences) < time_steps:
            sequences = sequences + [0]*(time_steps - len(sequences)) # Add zeros to dead processes so that the calculation of the mean is correct
            dying_processes[f"Sim. {j}"] = sequences
            all_processes[f"Sim. {j}"] = sequences
            
        else: 
            all_processes[f"Sim. {j}"] = sequences
            surviving_processes[f"Sim. {j}"] = sequences
        print(f"simulation {j} done")
       
    df_dying = pd.DataFrame.from_dict(dying_processes)  
    df_surviving = pd.DataFrame.from_dict(surviving_processes)  
    df_all = pd.DataFrame.from_dict(all_processes)  
     
    mean_seq_array_die = np.array(df_dying.iloc[:, 1:].mean(axis=1))
    mean_seq_array_survive = np.array(df_surviving.iloc[:, 1:].mean(axis=1))
    mean_seq_array_all = np.array(df_all.iloc[:, 1:].mean(axis=1))

    #Do not plot dead processes further
    zero_indices = np.where(mean_seq_array_die==0)[0]
    if zero_indices.size > 0:
        index_mean_time_death = zero_indices[0]
        mean_time_death = time_data[index_mean_time_death]
        mean_seq_array_die[mean_seq_array_die==0]=np.nan
        mean_seq_array_die[index_mean_time_death] = 0
        plt.plot(time, mean_seq_array_die, c = "#FF3333", label=f"Dying processes <Seq>(t) mean time of death: {round(mean_time_death, 2)}") 
    else:
        plt.plot(time, mean_seq_array_die, c = "#FF3333", label="Dying processes <Seq>(t)") 

    plt.plot(time, mean_seq_array_survive, c="#008000", label="Surviving processes <Seq>(t)") 
    plt.plot(time, mean_seq_array_all, label="All processes <Seq>(t)") 
    
    # Save data in a file
    all_data = {'Time': time}
    all_data["all processes"] = mean_seq_array_all
    all_data["surviving processes"] = mean_seq_array_survive
    all_data["dying processes"] = mean_seq_array_die
    df_all_data = pd.DataFrame.from_dict(all_data)
    
    column_names = "time\t\tall proc.\tsurv proc.\tdying proc."
    header = f"MEAN NUMBER OF DIFFERENT DNA SEQUENCES <Seq>(t) FOR ALL, SURVIVING AND DYING BRANCHING PROCESSES\
    \n\nNumber of simulations: {nsims}\nTime: {max_time}\nRNG seed: {seed}\
        \n\n{column_names}"
    np.savetxt(f"Mean_Seq_{nsims}_simulations_surv_dying.txt", df_all_data, header = header, delimiter="\t", fmt="%e")
    
    # Finish the plot
    plt.xlabel('time', fontsize=16)
    plt.ylabel("<Seq>", fontsize=16)
    plt.ylim(0, 7)
    #plt.xlim(0, 25)
    main_title = "Average number of different DNA sequences of all, surviving and dying processes"
    sub_title = f"{nsims} simulations, p_branch={p_branch}, p_mut={p_mut}, DNA_length={DNA_length}"    
    plt.suptitle(main_title, fontsize=22, color='black', y=1.05)
    plt.text(0.5, 0.98, sub_title, fontsize=18, color='#333333', ha='center', va='center', transform=plt.gcf().transFigure)
    plt.legend(fontsize=18)
    plt.savefig("33.png", dpi=300, bbox_inches='tight')
    plt.show()
    return 0

