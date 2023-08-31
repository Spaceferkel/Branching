from DNA_new import *
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit


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

    return ID_count_current_alive, ID_seq



def polynomial_function(x, a, b, c):
    return a*x**b+c


def linear_function(x, slope, intercept):
    return slope * x + intercept



def nsims_barplot(p_branch, max_time, nsims, p_mut, DNA_seq, rate=1):
    # IDEA: Plot a barplot for the average numbers of occurence of different mutated DNA sequneces for nsims simulations
    
    # Generate a DataFrame with contains the IDs and counts for every run
    df = pd.DataFrame(columns=['ID'])
    for j in range(nsims):
        dic, ID_seq = BBP_DNA_seq_number_at_t(p_branch, max_time, p_mut, DNA_seq, rate)
        for ID, count in dic.items():
            if ID not in df['ID'].values:
                df = df.append({'ID': ID}, ignore_index=True)
            row = df[df['ID'] == ID].index[0]
            df.at[row, j] = count
    df.fillna(0, inplace=True)
    
    # Calculate the mean number of occurence per ID
    mean_array = np.array(df.iloc[:, 1:].mean(axis=1))
    ID_array = df.iloc[:, 0]
    
    # Sort the arrays
    sorted_indices = np.argsort(mean_array)[::-1]
    sorted_mean_array = np.array(mean_array)[sorted_indices]
    sorted_ID_array = np.array(ID_array)[sorted_indices]
    ascending_list = np.array(list(range(1, len(sorted_mean_array))))

    
    # Plot 1 (Barplot)
    plt.figure(figsize=(20, 8)) 
    ax = plt.gca() 
    ax.set_facecolor('#F0F0F0')
    plt.xlabel('Different DNA sequences', fontsize=16)
    plt.ylabel('Average occurence <DNA_seq>', fontsize=16)
    title_line1 = f'Average occurrences of different mutated DNA sequences'
    plt.text(0.5, 1.15, title_line1, fontsize=18, color='black', ha='center', va='center', transform=plt.gca().transAxes)
    title_line2 = f'for {nsims} simulations, p_branch = {p_branch}, p_mut = {p_mut}, DNA_length = {len(DNA_seq)} bp at time {max_time}'
    plt.text(0.5, 1.08, title_line2, fontsize=14, color='gray', ha='center', va='center', transform=plt.gca().transAxes)
    plt.grid(True, linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Fit a function
    init = [1, -1.7, -1]
    params, covariance = curve_fit(polynomial_function, ascending_list, sorted_mean_array[1:], p0=init)
    a, b, c = params
    y_fit = polynomial_function(ascending_list, a, b, c)

    plt.bar(ascending_list, sorted_mean_array[1:], alpha=0.7) 
    plt.plot(ascending_list, y_fit, color='red', alpha=0.5, label=f"Fitted function: f(x)=({round(a, 2)})*x^({round(b, 2)})-{round(abs(c), 2)}")  # Angepasste Funktion
    plt.legend(fontsize='xx-large')    
    plt.show()
    
    
    # Plot 2 (double logarithmic)
    plt.figure(figsize=(20, 8)) 
    ax = plt.gca() 
    ax.set_facecolor('#F0F0F0')
    plt.xlabel('Logarithm of different DNA sequences', fontsize=16)
    plt.ylabel('Logarithm of average occurence <DNA_seq>', fontsize=16)
    title_line1 = f'Double logarithmic plot of average occurrences of different mutated DNA sequences'
    plt.text(0.5, 1.15, title_line1, fontsize=18, color='black', ha='center', va='center', transform=plt.gca().transAxes)
    title_line2 = f'for {nsims} simulations, p_branch = {p_branch}, p_mut = {p_mut}, DNA_length = {len(DNA_seq)} bp at time {max_time}'
    plt.text(0.5, 1.08, title_line2, fontsize=14, color='gray', ha='center', va='center', transform=plt.gca().transAxes)
    plt.grid(True, linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.scatter(np.log(ascending_list), np.log(y_fit), edgecolors='red', marker='o', facecolors='none', alpha=0.5, label=f"log(f(x))")  # Angepasste Funktion
    
    # Fit linear function 
    end_value = 3.5
    selected_ascending_list = np.log(ascending_list)[np.log(ascending_list) <= end_value]
    selected_sorted_mean_array = np.log(sorted_mean_array[1:len(selected_ascending_list)+1])
    p0 = (0, selected_sorted_mean_array[0]) 
    params, covariance = curve_fit(linear_function, selected_ascending_list, selected_sorted_mean_array, p0=p0)
    slope_fit, intercept_fit = params
    linear_fit = linear_function(selected_ascending_list, slope_fit, intercept_fit)
    plt.plot(selected_ascending_list, linear_fit, color='black', alpha=0.5, label=f"Linear Fit with slope {round(slope_fit, 2)}")

    plt.legend(fontsize='xx-large')    
    plt.show()
    
    return 0
