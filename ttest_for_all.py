import pandas as pd
import scipy.stats as stats
import os
import csv

# 文件夹路径，包含多个 txt 文件
folder_path = "/home/bbian/Data_all_calibrate/result/correlation_analysis_idrs_NC_term/"

# 创建一个空的数据框来存储所有文件的结果
all_results_df = pd.DataFrame(columns=["Column1", "Column2", "P_value", "Mean_value"])

# 遍历文件夹中的所有文件
for file_name in os.listdir(folder_path):
    if file_name.endswith(".txt"):
        # 提取文件名中的 "T.brucei" 和 "C_termi"
        file_parts = file_name.split(".")
        prefix1 = file_parts[0]
        prefix2 = file_parts[-2]

        # 读取数据
        file_path = os.path.join(folder_path, file_name)
        corre_distri_new1 = pd.read_csv(file_path, header=None, sep="\t", names=["Gene_name", "spearman_coefficient"])
        corre_distri_new1 = corre_distri_new1.dropna()
        # 提取 spearman_coefficient 列
        coeff_new1 = corre_distri_new1["spearman_coefficient"]

        # 执行 t-检验
        t_test_result = stats.ttest_1samp(coeff_new1, 0)

        # 获取 p-value 和均值
        p_value = t_test_result.pvalue
        mean_value = coeff_new1.mean()

        # 创建临时数据框来存储当前文件的结果
        temp_result_df = pd.DataFrame({"Column1": [prefix1], "Column2": [prefix2], "P_value": [p_value], "Mean_value": [mean_value]})

        # 将当前文件的结果添加到总结果数据框
        all_results_df = pd.concat([all_results_df, temp_result_df], ignore_index=True)

# 保存所有结果到一个结果文件
result_file = "/home/bbian/Data_all_calibrate/result/correlation_analysis_idrs_NC_term/ttest/all_results.txt"
all_results_df.to_csv(result_file, sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC, quotechar='"')

