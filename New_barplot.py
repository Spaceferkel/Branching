from DNA_new import *
import pandas as pd
from datetime import datetime
import os 
import matplotlib.pyplot as plt



def BBP_DNA_seq_number_at_t(p_branch, max_time, p_mut, DNA_seq, rate=1):
    # IDEA: returns a dictionary of the living IDs (with counts) at max_time and a DataFrame which contains the associated DNA sequences
    # NOTE: Each different DNA sequence has its own ID (and not every cell)
    
    N = 1
    ID_ancestor = generate_ID(DNA_seq)
    ID_seq = {ID_ancestor:DNA_seq}           # Contains IDs and sequences of all DNA sequences which have ever existed
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


def BBP_DNA_seq_number_size(p_branch, max_time, p_mut, DNA_seq, rate=1):
    # IDEA: returns a dictionary with the ID and count of DNA sequences which have ever existed
    # NOTE: Each different DNA sequence has its own ID (and not every cell)
    
    N = 1
    ID_ancestor = generate_ID(DNA_seq)
    ID_seq = {ID_ancestor:DNA_seq}           # Contains IDs and sequences of all DNA sequences which have ever existed
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
            # If it doesn't branch, cell dies
            N -= 1

        current_time += waiting_time
        
    return ID_count_current_alive, ID_ancestor




def nsims_distribution(p_branch, max_time, nsims, p_mut, DNA_seq, rate, mode="snapshot"):
    # IDEA: Generate data for the average numbers of occurence of different mutated DNA sequences for nsims simulations
    
    # Define type of the generated data
    if mode == "snapshot":
        generating_function = BBP_DNA_seq_number_at_t
    elif mode == "size":
        generating_function = BBP_DNA_seq_number_size
    else:
        raise ValueError("Invalid mode entered")
        
    # Get time of the start of program execution
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%H-%M-%S")

    # Generate a DataFrame wich contains the data from each simulation run
    data_all = {}
    for j in range(nsims):
        dic, ID_ancestor = generating_function(p_branch, max_time, p_mut, DNA_seq, rate)
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

    # Save the data in a file 
    raw_data = np.column_stack([ascending_list, mean_array])
    folder = "Data" 
    file_name = f"barplot_{nsims}_{p_branch}_{dt_string}.txt"
    file_path = os.path.join(os.getcwd(), folder, file_name)
    column_labels = "\t".join(["Rank", "Mean frequency"])
    time_format = now.strftime("%d/%m/%Y, %H:%M:%S")
    header = f" MEAN FREQUENCY OF THE I-TH FREQUENT MUTATED DNA SEQUENCE\
        \n\n Program execution time: \t {time_format} \n Mode: \t\t\t {mode}\
        \n Number of simulations: \t {nsims} \n Branching probability: \t {p_branch}\
        \n At time: \t\t\t {max_time}\n Rate: \t\t\t {rate}\n DNA length: \t\t\t {len(DNA_seq)} \n Mutation probability: \t {p_mut}\
        \n RNG seed: \t\t\t {seed}\n\n\n{column_labels}" 
    np.savetxt(file_path, raw_data, fmt = ['%d', '%e'], header=header, delimiter="\t")
    
    return 0


def exponential_function(x, a, b):
    return a*np.exp(b*x)

def polynomial_function(x, a, b):
    return a*x**b


def plot_barplot(file_name):
    # IDEA: Plot the data in a barplot
    
    # Get the data from the file
    file_path = 'Data\\' + file_name + ".txt"
    data = np.genfromtxt(file_path, delimiter='\t', skip_header=13)
    ascending_list = data[:, 0]
    mean_array = data[:, 1]
    
    # Get the parameters from the file
    parameters = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines[3:10]:
        key, value = line.split(":")
        key = key.strip("# ").strip()
        value = value.strip()
        parameters[key] = value

    mode = parameters["Mode"]
    nsims = int(parameters["Number of simulations"])
    p_branch = float(parameters["Branching probability"])
    max_time = int(parameters["At time"])
    rate = int(parameters["Rate"])
    DNA_length = int(parameters["DNA length"])
    p_mut = float(parameters["Mutation probability"])
    
    
    
    """
    # Fit exponential
    init = [0.1, -10]
    params, covariance = curve_fit(exponential_function, ascending_list, mean_array, p0=init, maxfev=2000)
    a, b = params
    y_exp_fit = exponential_function(ascending_list, a, b)
    total_sum_of_squares = np.sum((mean_array - np.mean(mean_array))**2)
    residual_sum_of_squares = np.sum((mean_array - y_exp_fit)**2)
    r_squared = 1 - (residual_sum_of_squares / total_sum_of_squares)
    plt.plot(ascending_list, y_exp_fit, color='green', alpha=0.5, label=f"Fitted exponential function: f(x)=({round(a, 2)})*exp({round(b, 2)}*x), R^2={round(r_squared,2)}")  # Angepasste Funktion

    # Fit power function 
    init = [2, -1]
    params, covariance = curve_fit(polynomial_function, ascending_list, mean_array, p0=init, maxfev=2000)
    c, d = params
    y_polynomial_fit = polynomial_function(ascending_list, c, d)
    residual_sum_of_squares = np.sum((mean_array - y_polynomial_fit)**2)
    r_squared = 1 - (residual_sum_of_squares / total_sum_of_squares)
    plt.plot(ascending_list, y_polynomial_fit, color='red', alpha=0.5, label=f"Fitted power function: f(x)=({round(c, 2)})*x^({round(d, 2)}), R^2={round(r_squared,2)}")  # Angepasste Funktion
    """
    
    plt.plot(ascending_list, mean_array, label=f"Mode: {mode}", alpha=0.7) 
    
    
    
    
    return 0

