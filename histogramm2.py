from DNA_new import *
import matplotlib.pyplot as plt
import pandas as pd


def BBP_DNA_seq_number_at_t(p_branch, max_time, p_mut, DNA_seq, rate=1):
    # IDEA: returns a dictionary of the living IDs (with counts) at max_time and a DataFrame which contains the associated DNA sequences
    # NOTE: Each different DNA sequence has its own ID (and not every cell)
    
    generator = np.random.default_rng()
    N = 1
    ID_ancestor = generate_ID(DNA_seq)
    ID_seq = {ID_ancestor:DNA_seq} # Contains IDs and sequences of all DNA sequences which have ever existed
    ID_count_current_alive = {ID_ancestor:1} # Contains IDs and number of those sequences which are currently alive
    current_time = 0

    while current_time <= max_time and N > 0:
        waiting_time = generator.exponential(1/(N*rate))
        counts_DNA_alive = np.array(list(ID_count_current_alive.values()))
        p_distribution = counts_DNA_alive/np.sum(counts_DNA_alive)
        ID_parent = generator.choice(list(ID_count_current_alive), p=p_distribution) # Choose ID of one living DNA sequence with given probability distribution

        if branch_or_not(p_branch) == True:
            N += 1
            p_random_mut = generator.random()
            
            if p_random_mut < p_mut: 
                # The sequence of the parent gets passed on and an additional mutant sequence is added 
                new_seq = mutate(ID_seq[ID_parent])
                ID_new = generate_ID(new_seq)
                
                if ID_new in ID_count_current_alive:
                    ID_count_current_alive[ID_new] = ID_count_current_alive[ID_new] + 1
                else:
                    ID_count_current_alive[ID_new] = 1
                    ID_seq[ID_new] = new_seq   
            else:
                ID_count_current_alive[ID_parent] = ID_count_current_alive[ID_parent] + 1 # If it doesn't mutate, one identical sequence is added
        else:
            # If it doesn't branch, DNA sequence dies
            N -= 1
            if ID_count_current_alive[ID_parent] > 1: #if multiple identical sequences exist, just lower their count
                ID_count_current_alive[ID_parent] = ID_count_current_alive[ID_parent] - 1
            else:
                del ID_count_current_alive[ID_parent]
                
        current_time += waiting_time

    return ID_count_current_alive, ID_ancestor


def nsims_barplot(p_branch, max_time, nsims, p_mut, DNA_seq, rate=1):
    # IDEA: Plot a barplot for the average numbers of occurrence of different mutated DNA sequences for nsims simulations
    
    # Generate a DataFrame which contains the data from each simulation run
    data_all = {}
    for j in range(nsims):
        dic, ID_ancestor = BBP_DNA_seq_number_at_t(p_branch, max_time, p_mut, DNA_seq, rate)
        if ID_ancestor in dic: 
            del dic[ID_ancestor]
        if dic: 
            data = np.array(list(dic.values()))/(sum(dic.values())) # Get a sorted array of normalized occurence of each sequence
            data_all[j] = np.sort(data)
        print(f"Simulation {j} done")
        
    df = pd.DataFrame.from_dict(data_all, orient='index').transpose()
    df.fillna(0, inplace=True)
    
    # Calculate the sorted mean frequency
    mean_array = np.sort(np.array(df.mean(axis=1)))[::-1] # Do not take ancestor DNA sequence
    ascending_list = np.array(list(range(1, len(mean_array)+1)))
  
    # Plot Barplot
    plt.figure(figsize=(20, 8)) 
    ax = plt.gca() 
    ax.set_facecolor('#F0F0F0')
    plt.xlabel('frequency rank', fontsize=16)
    plt.ylabel('mean frequency', fontsize=16)
    #plt.yscale("log")
    #plt.xscale("log")
    title_line1 = f'Mean frequency of the i-th frequent mutated DNA sequence (snapshot)'
    plt.text(0.5, 1.15, title_line1, fontsize=18, color='black', ha='center', va='center', transform=plt.gca().transAxes)
    title_line2 = f'for {nsims} simulations, p_branch = {p_branch}, p_mut = {p_mut}, DNA_length = {len(DNA_seq)} bp at time {max_time}'
    plt.text(0.5, 1.08, title_line2, fontsize=14, color='gray', ha='center', va='center', transform=plt.gca().transAxes)
    plt.grid(True, linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
        
    plt.bar(ascending_list, mean_array, alpha=0.7) 
    plt.legend(fontsize='xx-large')   
    plt.savefig("distribution12.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    return 0
