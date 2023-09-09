import sys
import os
import re
import os.path
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from matplotlib.ticker import MultipleLocator

meta_gene_scan_matrix = []


meta_gene_scan = [0] * 61
label = list(range(61))
numbersite = 0

input_data_wave = sys.argv[1]
input_path = sys.argv[2]
with open(input_data_wave) as f:
    for line in f.readlines():
        line = line.split()
        geneid = line[0]
        position = []
        amino_acid = []
        charges = []
        reads_counts = []

        my_file = Path("/home/bbian/Data_all_calibrate/result/merge_charge_reads_each_gene/"+input_path+"/{}.txt".format(geneid))
        if not my_file.is_file():
            continue
        with open(my_file, "r") as t:
            for lines in t.readlines():
                line2 = lines.split()
                position.append(int(line2[0]))
                amino_acid.append(line2[1])
                charges.append(line2[2])
                reads_counts.append(float(line2[3]))

            for i in range(len(reads_counts)):
                index1 = position[i] - 30
                index2 = position[i] + 30
                if amino_acid[i] == "P" and index1 > 0 and index2 < len(reads_counts):
                    new_reads = reads_counts[index1 - 1:index2]
                    numbersite += 1 
                    meta_gene_scan_matrix.append(new_reads[:])
                    for j in range(len(meta_gene_scan)):
                        meta_gene_scan[j] += new_reads[j]
                    
    print(numbersite)
                    
    meta_gene_scan_matrix=np.array(meta_gene_scan_matrix)
    #print(meta_gene_scan_matrix[:,0])

    mean_values = [value / numbersite for value in meta_gene_scan]
    print(mean_values)
    
    lower_ci=[]
    upper_ci=[]
    for i in range (61):
        std_error_new=st.t.interval(alpha=0.95, df=len(meta_gene_scan_matrix[:,i])-1,
                  loc=np.mean(meta_gene_scan_matrix[:,i]),
                  scale=st.sem(meta_gene_scan_matrix[:,i]))
        lower_ci.append(std_error_new[0])
        upper_ci.append(std_error_new[1])
    print(std_error_new)

    #lower_ci = [mean - 1.96 * std_err for mean, std_err in zip(mean_values, std_errors)]
    #upper_ci = [mean + 1.96 * std_err for mean, std_err in zip(mean_values, std_errors)]
    print(lower_ci)
    print(upper_ci)

    # Rest of the code for plotting and exporting results
    # Create the plot

plt.figure(figsize=(10, 6),dpi=300)
plt.plot(label, mean_values, color='green', label='Mean')
plt.fill_between(label, lower_ci, upper_ci, color='green', alpha=0.2, label='95% CI')
plt.xlabel('Position')
plt.ylabel('Mean Read Counts')
plt.title('Metagene Plot')
#x_ticks = np.arange(min(label), max(label)+1,5)
#plt.locator_params(axis='x', nbins=4)
x_ticks = np.arange(0, 61, 10)  # -30 to 30, with interval of 5
plt.xticks(x_ticks, [str(i-30) for i in x_ticks])  # Convert to strings for labeling

plt.axvline(x=30, color='gray', linestyle='--')
#plt.xticks(label, [str(i-30) for i in label])  # Change x-axis ticks to -30 to 30
#plt.axvline(x=30, color='gray', linestyle='--')  # Add a vertical line at the center (position 0)
plt.legend()

# Save or show the plot
plt.savefig("/home/bbian/Data_all_calibrate/result/metagene_plot/"+input_path+".proline.metagene_plot.png")
plt.show()



