import pandas as pd
import scipy.stats as stats
import os
import csv


folder_path = "/Your/work/path/result/partial_correlation_analysis/"


all_results_df = pd.DataFrame(columns=["Column1", "Column2", "P_value", "Mean_value"])


for file_name in os.listdir(folder_path):
    if file_name.endswith(".txt"):
        file_parts = file_name.split(".")
        prefix1 = file_parts[0]
        prefix2 = file_parts[-2]

     
        file_path = os.path.join(folder_path, file_name)
        corre_distri_new1 = pd.read_csv(file_path, header=None, sep="\t", names=["Gene_name", "spearman_coefficient"])
        corre_distri_new1 = corre_distri_new1.dropna()
     
        coeff_new1 = corre_distri_new1["spearman_coefficient"]

   
        t_test_result = stats.ttest_1samp(coeff_new1, 0)

   
        p_value = t_test_result.pvalue
        mean_value = coeff_new1.mean()

    
        temp_result_df = pd.DataFrame({"Column1": [prefix1], "Column2": [prefix2], "P_value": [p_value], "Mean_value": [mean_value]})

  
        all_results_df = pd.concat([all_results_df, temp_result_df], ignore_index=True)


result_file = "/Your/work/path/result/partial_correlation_analysis/ttest/all_results.txt"
all_results_df.to_csv(result_file, sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC, quotechar='"')

