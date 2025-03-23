import pandas as pd


import statsmodels.stats.multitest as stats
import math
from tqdm import tqdm
import numpy as np

# def detect_zero_sum_rows(input_dataframe):
#     # Extract only the columns corresponding to allele counts (A, T, C, G)
#     allele_columns = input_dataframe.iloc[:, 4:]  # Ignoring the first 4 columns: chr, pos, index, ref_n
#     # Group the data in sets of 4 columns (A, T, C, G) per sample
#     num_samples = allele_columns.shape[1] // 4
    
#     # Create 3D numpy array with rows (genomic positions), columns (samples), allele counts (ATCG) as dimensions
#     reshaped_data = allele_columns.to_numpy().reshape(allele_columns.shape[0], num_samples, 4)
    
#     # Sum ATCG values for each sample
#     allele_sums = reshaped_data.sum(axis=2)
    
#     # Extract the positions and rows where any sample has an allele sum of zero
#     zero_sum_positions = input_dataframe[np.any(allele_sums == 0, axis=1)]
#     zero_sum_indices = zero_sum_positions.index.tolist()    
    
#     # Display first few rows with zero allele sum for inspection
#     # print(zero_sum_positions.head())
    
#     all_values = allele_columns.to_numpy()
#     # -----------------------------------------------------------
#     # For logging purposes:
#     # Note: This is the total count of ATCG sums equaling zero. 
#     # This means that we might have many zero sum instances in the same position. 
    
    
#     zero_count = np.sum(allele_sums == 0) 
#     nan_count = np.sum(np.isnan(all_values))
#     print(f'nan count: {nan_count}')
    
    
#     # Create the log folder if it doesn't already exist
#     # log_folder_path = '\zero_sum_logs'
#     # os.makedirs(log_folder_path, exist_ok=True)
#     # # # Define log file path
#     # log_path = os.path.join(log_folder_path, f'{maf_filename}.txt')
    
#     # # Write the total ATCG zero sum instances count and total nan count to the combined log file
#     # with open(log_path, 'w') as log_file:
#     #     log_file.write(f"Number of rows with zero allele sum: {zero_sum_positions.shape[0]}")
#     #     log_file.write(f"Total ATCG zero sum instances count: {zero_count}")
#     #     log_file.write(f"Total nan count: {nan_count}")
    
#     print(f"Number of rows with zero allele sum: {zero_sum_positions.shape[0]}")
#     print(f"Total ATCG zero sum instances count: {zero_count}")
#     # print(f"Total nan count: {nan_count}")
#     # TODO: Create output log files
#     # -----------------------------------------------------------
    
#     return zero_sum_indices
# TEST BELOW
chromosome=rf'LR537121'

input_baypass_path = rf'D:/final_baypass_results_resorted/Final_Sparus_{chromosome}_baypass_summary_pi_xtx.out'
input_maf_path = rf'final_maf_files/maf_{chromosome}.csv'
output_path = f'D:/final_baypass_results_resorted/Sparus_{chromosome}_baypass_result.csv'
# input_baypass_path = rf'Sorted_Results/Sparus_{chromosome}_baypass_summary_pi_xtx.out'
# input_maf_path = rf'maf_{chromosome}.csv'
# output_path = f'Sparus_{chromosome}_baypass_result.csv'

BayPass_sparus_result = pd.read_csv(input_baypass_path,  sep='\s+')
BayPass_sparus_result
BayPass_sparus_result.rename(columns = {0:'MRK', 1:'M_P', 2:'SD_P', 3:'M_XtX', 4:'SD_XtX', 5:'XtXst', 6:'log10(1/pval)'}, inplace=True)

#---------------------
# Constants
for chromosome_num in tqdm(range(21,45), desc="Processing Chromosomes", unit="chromosome", colour="blue"):
    chromosome=rf'LR5371{chromosome_num}'
    # adding these lines to sort baypass based on maf pos order, check maf_prep.py for details
    # Sparus_maf = pd.read_csv(rf'maf_{chromosome}.csv', sep=",", header=1)
    # Sparus_maf.drop(index = detect_zero_sum_rows(Sparus_maf), inplace = True)
    # Sparus_maf.reset_index(drop=True, inplace=True)
    # post_list = Sparus_maf['pos'].tolist()
    # Sparus_maf.sort_values('pos', inplace=True)
    # Sparus_maf.reset_index(drop=True, inplace=True)
    # Find way to remove popoolation threshold rows here from maf file
    #---------------------------
    # GET POPOOLATION DF, order it first and then remove the rows we don't want
    # genotype_file = rf'maf_{chromosome}_popoolation.genotype'
    # snp_df = pd.read_csv(rf'../final_prep/{genotype_file}', sep=" ", header=None)
    # snp_df.insert(loc = 0,
    #             column = 'temp_pos',
    #             value = post_list
    # )
    # snp_df.sort_values('temp_pos', inplace=True)
    # snp_df.drop(columns='temp_pos', inplace = True)
    # # Resetting index to not remove the wrong rows in the next step
    # snp_df.reset_index(drop=True, inplace=True)
    # print("Calculating population sums...")
    # #Switch applymap with map below (FutureWarning)
    # population_sums = snp_df.iloc[:, 3:].applymap(
    #     lambda x: sum(int(num) for num in x.split(':'))
    # )
    # rows_to_remove = (population_sums <= 54).any(axis=1)
    
    # # Remove rows from dataframes
    # df_removed_rows = snp_df[rows_to_remove]
    
    # snp_df = snp_df[~rows_to_remove]
    # Sparus_maf = Sparus_maf[~rows_to_remove]
    
    # removed_rows_count = df_removed_rows.shape[0]
    # print(f"Removed {removed_rows_count} rows from maf and popoolation files. (Threshold 54)")
    #------------------------------------
    input_baypass_path = rf'D:/final_baypass_results_resorted/Final_Sparus_{chromosome}_baypass_summary_pi_xtx.out'
    input_maf_path = rf'final_maf_files/maf_{chromosome}.csv'
    output_path = f'D:/final_baypass_results_resorted/Sparus_{chromosome}_baypass_result.csv'
    
    # input_baypass_path = rf'final_baypass_result_files_resorted/Sparus_{chromosome}_baypass_summary_pi_xtx.out'
    # input_maf_path = rf'final_maf_files/maf_{chromosome}.csv'
    # output_path = f'final_baypass_files_ready_for_manhattan_with_pos/Sparus_{chromosome}_baypass_result.csv'
    BayPass_sparus_result = pd.read_csv(input_baypass_path,  sep='\s+')
    # BayPass_sparus_result = pd.read_csv(input_baypass_path,  sep=',', header= None)
    # BayPass_sparus_result = pd.read_csv(input_baypass_path,  header=None, skiprows=1, delim_whitespace=True, engine='python')
    # Load data
    BayPass_sparus_result.rename(columns = {0:'MRK', 1:'M_P', 2:'SD_P', 3:'M_XtX', 4:'SD_XtX', 5:'XtXst', 6:'log10(1/pval)'}, inplace=True)
    
    # new addition 30/12/2024 - adding the pos column from the maf file to the baypass file
    Sparus_maf= pd.read_csv(input_maf_path)
    pos = Sparus_maf['pos']
    BayPass_sparus_result['pos'] = pos
    # BayPass_sparus_result.sort_values('temp_pos', inplace=True)
    # BayPass_sparus_result.drop(columns='temp_pos', inplace = True)
    # Resetting index to not remove the wrong rows in the next step
    # BayPass_sparus_result.reset_index(drop=True, inplace=True)
    
    # Define your column names
    column_names = [
        "chr", "pos", "index", "ref_n", 
        "A", "C", "G", "T", "A.1", "C.1", "G.1", "T.1", 
        "A.2", "C.2", "G.2", "T.2", "A.3", "C.3", "G.3", "T.3", 
        "A.4", "C.4", "G.4", "T.4", "A.5", "C.5", "G.5", "T.5", 
        "A.6", "C.6", "G.6", "T.6", "A.7", "C.7", "G.7", "T.7", 
        "A.8", "C.8", "G.8", "T.8", "A.9", "C.9", "G.9", "T.9", 
        "A.10", "C.10", "G.10", "T.10", "A.11", "C.11", "G.11", "T.11", 
        "A.12", "C.12", "G.12", "T.12", "A.13", "C.13", "G.13", "T.13", 
        "A.14", "C.14", "G.14", "T.14", "A.15", "C.15", "G.15", "T.15", 
        "A.16", "C.16", "G.16", "T.16", "A.17", "C.17", "G.17", "T.17", 
        "A.18", "C.18", "G.18", "T.18", "A.19", "C.19", "G.19", "T.19"
    ]
    
    # check on "inf" log10(1/pval) - they are 6 values for gilthead seabream (36 values for labrax in previous study)
    # emailed the author of BayPass (Mathieu Gautier)
    # answer: if you get inf, it means that p was smaller than machine precision (or rahter 1/p higher than machine precision),
    # so you can safely transform such log10(1/pval) values to the highest observed value + 1
    # BayPass_sparus_result.sort_values(by=['log10(1/pval)']).tail(20)
    
    # 16.653560 was added here for gilthead seabream
    # convert the 6 cases of "inf" in log10(1/pval) to the highest observed value + 1
    # BayPass_result.replace({float("inf"): 16.653560}, inplace=True)                                               
    # BayPass_result.sort_values(by=['log10(1/pval)']).iloc[161825:161831, :]
    BayPass_sparus_result.replace({float("inf"): 16.653560}, inplace=True)
    # BayPass_sparus_result.sort_values(by=['log10(1/pval)']).tail(20)

    # Sparus_maf.columns = column_names
    # print(Sparus_maf.head())
    #  BayPass_sparus_result['pos'] = Sparus_maf['pos']
    BayPass_sparus_result["delog_pval"] = 10 ** - BayPass_sparus_result["log10(1/pval)"]

    # multiple test correction - and repeat the analysis above - check in biostatistics
    # is_sorted = False?
    BayPass_sparus_result_BH = stats.multipletests(list(BayPass_sparus_result['delog_pval']), alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    non_zero_values = []
    for i in list(BayPass_sparus_result_BH[1]):
        if i > 0:
            non_zero_values.append(i)
    minimum_non_zero_value = min(non_zero_values)

    log_transformed_values = []
    for i in list(BayPass_sparus_result_BH[1]):
        try:
            log_transformed_values.append(-math.log10(i))
        except ValueError:
            log_transformed_values.append(-math.log10(minimum_non_zero_value))

    # append the log-transformed values as column to the df of interest
    BayPass_sparus_result["log10(1/pval)_BH"] = log_transformed_values
    BayPass_sparus_result['chromosome'] = str(chromosome)
    # TODO - ADD POS COLUMN FROM THE FINAL MAF FILES TO THE FINAL BAYPASS FILES 
    BayPass_sparus_result.to_csv(output_path, sep=",", index=False)





# merge final baypass files

# Initialize an empty list to hold dataframes
df_list = []
# filename_test = r"Sparus_LR537122_baypass_result.csv"
# temp_df = pd.read_csv(filename_test)
# temp_df.head
# temp_df.sort_values('pos')


# filename_test = r"maf_LR537121.csv"
# temp_df = pd.read_csv(filename_test, header= )
# temp_df.head

# temp_df.sort_values('1')

# one position is NaN? Check if the positions match up or not
# temp_df.to_csv("testcsv")
# temp_df.isna().sum()
# Loop over the range of file suffixes

for i in range(21, 44 + 1):
    # Format the filename
    filename = rf"D:/final_baypass_results_resorted/Sparus_LR5371{i}_baypass_result.csv"
    # filename = rf"final_baypass_files_ready_for_manhattan_with_pos/Sparus_LR5371{i}_baypass_result.csv"
    # Read the CSV file
    temp_df = pd.read_csv(filename)
    # Before merging these files we will need to sort them by their positions to get accurate manhattan plot peaks
    # Popoolation was removing any non sequential row so we will need to apply the same method there and reprocess everything
    # temp_df.sort_values()
    # Append the dataframe to the list
    df_list.append(temp_df)

# Concatenate all dataframes in the list
merged_df = pd.concat(df_list, ignore_index=True)

# Optionally, save the merged DataFrame to a new CSV file
merged_df.to_csv("D:/final_baypass_results_resorted/merged_baypass_results.csv", index=False)
# merged_df.to_csv("final_baypass_files_ready_for_manhattan_with_pos/merged_baypass_results.csv", index=False)

print("Merging completed. Merged file saved as 'merged_baypass_results.csv'.")





#------------------------------------------------------------------------------------Testing

chromosome=rf'LR537121'

input_baypass_path = rf'Sparus_{chromosome}_baypass_summary_pi_xtx.out'
input_maf_path = rf'maf_{chromosome}.csv'
output_path = f'Resorted_Results/Sparus_{chromosome}_baypass_result.csv'

# Load data
BayPass_sparus_result = pd.read_csv(input_baypass_path,  header=None, skiprows=1, delim_whitespace=True, engine='python')
BayPass_sparus_result.rename(columns = {0:'MRK', 1:'M_P', 2:'SD_P', 3:'M_XtX', 4:'SD_XtX', 5:'XtXst', 6:'log10(1/pval)'}, inplace=True)

Sparus_maf= pd.read_csv(input_maf_path)

# Define your column names
column_names = [
    "chr", "pos", "index", "ref_n", 
    "A", "C", "G", "T", "A.1", "C.1", "G.1", "T.1", 
    "A.2", "C.2", "G.2", "T.2", "A.3", "C.3", "G.3", "T.3", 
    "A.4", "C.4", "G.4", "T.4", "A.5", "C.5", "G.5", "T.5", 
    "A.6", "C.6", "G.6", "T.6", "A.7", "C.7", "G.7", "T.7", 
    "A.8", "C.8", "G.8", "T.8", "A.9", "C.9", "G.9", "T.9", 
    "A.10", "C.10", "G.10", "T.10", "A.11", "C.11", "G.11", "T.11", 
    "A.12", "C.12", "G.12", "T.12", "A.13", "C.13", "G.13", "T.13", 
    "A.14", "C.14", "G.14", "T.14", "A.15", "C.15", "G.15", "T.15", 
    "A.16", "C.16", "G.16", "T.16", "A.17", "C.17", "G.17", "T.17", 
    "A.18", "C.18", "G.18", "T.18", "A.19", "C.19", "G.19", "T.19"
]
Sparus_maf.columns = column_names
# found the missing value issue, at some point
# in the conversion from maf to sparus file, the first pos row
# gets deleted. This might be because of a header=True/False issue
print(Sparus_maf['pos'].dtype)
# Sparus_maf = Sparus_maf.dropna(subset=['pos'])
Sparus_maf
print(Sparus_maf.head())
BayPass_sparus_result['pos'] = Sparus_maf['pos'].astype(int)
BayPass_sparus_result["delog_pval"] = 10 ** - BayPass_sparus_result["log10(1/pval)"]
BayPass_sparus_result = BayPass_sparus_result.dropna(subset=['pos'])
print(BayPass_sparus_result['pos'].dtype)
BayPass_sparus_result['pos'] = BayPass_sparus_result['pos'].astype(int)

# multiple test correction - and repeat the analysis above
# is_sorted = False?
BayPass_sparus_result_BH = stats.multipletests(list(BayPass_sparus_result['delog_pval']), alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
# BayPass_sparus_result_BH.head()
non_zero_values = []
for i in list(BayPass_sparus_result_BH[1]):
    if i > 0:
        non_zero_values.append(i)
minimum_non_zero_value = min(non_zero_values)

log_transformed_values = []
for i in list(BayPass_sparus_result_BH[1]):
    try:
        log_transformed_values.append(-math.log10(i))
    except ValueError:
        log_transformed_values.append(-math.log10(minimum_non_zero_value))

# append the log-transformed values as column to the df of interest
BayPass_sparus_result["log10(1/pval)_BH"] = log_transformed_values
BayPass_sparus_result['chromosome'] = str(chromosome)
BayPass_sparus_result.head()
BayPass_sparus_result.isna().sum()

BayPass_sparus_result.to_csv(output_path, sep=",", index=False)
