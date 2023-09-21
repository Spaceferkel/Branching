import seaborn as sns
import matplotlib.pyplot as plt
from BBP_N_new import *
import pandas as pd
import numpy as np

def plot_nsims_multiplots(nsims, max_time, time_steps, p_values):
    
    # Prepare the plot
    sns.set_style("white")
    plt.figure(figsize=(12, 8))
    plt.subplots_adjust(wspace=0.2, hspace=0.5)
    with open("Quadruplot.txt", 'w') as file:
        file.write(f"SIMULATIONS OF BINARY BRANCHING PROCESSES \n\n Number of Simulations: {nsims} \n RNG seed: {seed}\n")
    
    # Create subplots
    rows = 2
    cols = 2
    axes = []
    time_data = np.linspace(0, max_time, time_steps)
    
    for idx, p_branch in enumerate(p_values, 1):
        ax = plt.subplot(rows, cols, idx)
        
        data = {"time": time_data}
        
        for j in range(nsims):
            time_data, population_data = BBP_N(p_branch, max_time, time_steps)
            if len(population_data) < time_steps:
                population_data = population_data + [np.nan]*(time_steps - len(population_data))
            data[f"Simulation {j}"] = population_data
            sns.lineplot(x=time_data, y=population_data, dashes=False, palette="rocket_r", ci=None)
            print(f"Simulation {j} done")

        df = pd.DataFrame.from_dict(data)
        with open("Quadruplot.txt", 'a') as file:
            np.savetxt(file, X=df, header = f"\n\nBranching probability: {p_branch}\n", delimiter="\t", fmt="%e")
 
        plt.xlabel('time', fontsize=13)
        plt.ylabel('Population', fontsize=13)
        plt.title(f'p_branch = {p_branch}', fontsize=14)
        plt.legend().set_visible(False)
        axes.append(ax)

    # Finish the plot
    main_title = "Binary Branching Processes"
    sub_title = f"Number of simulations: {nsims}"    
    plt.suptitle(main_title, fontsize=22, color='black', y=1.05)
    plt.text(0.5, 0.98, sub_title, fontsize=18, color='#333333', ha='center', va='center', transform=plt.gcf().transFigure)
    plt.savefig(f"BP_multiplot_{nsims}.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    return 0
    
