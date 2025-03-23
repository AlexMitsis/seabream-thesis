# We need to remove both the rows which are removed in the initial analysis as well as the rows which are removed after the results for popoolation
# Also sort all files based on pos
# STEP 1:

from preprocessing import BayPass, PCA, PoPoolation
import pandas as pd
import numpy as np
import warnings
import numpy as np
import os
import shutil
from tqdm import tqdm

maf_files = [
    'maf_LR537121.csv', 'maf_LR537122.csv', 'maf_LR537123.csv', 'maf_LR537124.csv',
    'maf_LR537125.csv', 'maf_LR537126.csv', 'maf_LR537127.csv', 'maf_LR537128.csv',
    'maf_LR537129.csv', 'maf_LR537130.csv', 'maf_LR537131.csv', 'maf_LR537132.csv',
    'maf_LR537133.csv', 'maf_LR537134.csv', 'maf_LR537135.csv', 'maf_LR537136.csv',
    'maf_LR537137.csv', 'maf_LR537138.csv', 'maf_LR537139.csv', 'maf_LR537140.csv',
    'maf_LR537141.csv', 'maf_LR537142.csv', 'maf_LR537143.csv', 'maf_LR537144.csv'
]
genotype_files = [maf_file.replace('.csv', '_popoolation.genotype') for maf_file in maf_files]
baypass_files = [maf_file.replace('.csv', '.genotype') for maf_file in maf_files]

def detect_zero_sum_rows(input_dataframe):
    # Extract only the columns corresponding to allele counts (A, T, C, G)
    allele_columns = input_dataframe.iloc[:, 4:]  # Ignoring the first 4 columns: chr, pos, index, ref_n
    # Group the data in sets of 4 columns (A, T, C, G) per sample
    num_samples = allele_columns.shape[1] // 4
    
    # Create 3D numpy array with rows (genomic positions), columns (samples), allele counts (ATCG) as dimensions
    reshaped_data = allele_columns.to_numpy().reshape(allele_columns.shape[0], num_samples, 4)
    
    # Sum ATCG values for each sample
    allele_sums = reshaped_data.sum(axis=2)
    
    # Extract the positions and rows where any sample has an allele sum of zero
    zero_sum_positions = input_dataframe[np.any(allele_sums == 0, axis=1)]
    zero_sum_indices = zero_sum_positions.index.tolist()    
    
    # Display first few rows with zero allele sum for inspection
    # print(zero_sum_positions.head())
    
    all_values = allele_columns.to_numpy()
    # -----------------------------------------------------------
    # For logging purposes:
    # Note: This is the total count of ATCG sums equaling zero. 
    # This means that we might have many zero sum instances in the same position. 
    
    
    zero_count = np.sum(allele_sums == 0) 
    nan_count = np.sum(np.isnan(all_values))
    print(f'nan count: {nan_count}')
    
    
    # Create the log folder if it doesn't already exist
    # log_folder_path = '\zero_sum_logs'
    # os.makedirs(log_folder_path, exist_ok=True)
    # # # Define log file path
    # log_path = os.path.join(log_folder_path, f'{maf_filename}.txt')
    
    # # Write the total ATCG zero sum instances count and total nan count to the combined log file
    # with open(log_path, 'w') as log_file:
    #     log_file.write(f"Number of rows with zero allele sum: {zero_sum_positions.shape[0]}")
    #     log_file.write(f"Total ATCG zero sum instances count: {zero_count}")
    #     log_file.write(f"Total nan count: {nan_count}")
    
    print(f"Number of rows with zero allele sum: {zero_sum_positions.shape[0]}")
    print(f"Total ATCG zero sum instances count: {zero_count}")
    # print(f"Total nan count: {nan_count}")
    # TODO: Create output log files
    # -----------------------------------------------------------
    
    return zero_sum_indices

# for maf_filename in tqdm(maf_files, desc="Processing Chromosomes", unit="chromosome", colour="blue"):
#     # Read the CSV, skipping the duplicate header row
#     biallelic_fish_df = pd.read_csv(f'{maf_filename}', sep=",", header=1)
#     biallelic_fish_df.drop(index = detect_zero_sum_rows(biallelic_fish_df), inplace = True)
#     biallelic_fish_df.to_csv(f'processed_files/{maf_filename}', sep=",", header=None, index=False)

# STEP 2: -----------------------------------------------------------------------------------------------------------------------------------------------------

# warnings.filterwarnings('ignore')  # not to print any warnings
# pd.set_option('display.max_rows', 30)  # Change 10 to your desired maximum number of rows
# pd.set_option('display.max_columns', 25)



# # global variables
# INFILE = "maf_LR537121_popoolation - Copy.genotype"
# SUM_THRESHOLD = 54

def optimized_remove_rows(maf_file, genotype_file, baypass_file):
    
    # GET MAF DF
    maf_df = pd.read_csv(maf_file, sep=",", header=1)
    maf_df.drop(index = detect_zero_sum_rows(maf_df), inplace = True)
    # Resetting index to not remove the wrong rows in the next step
    maf_df.reset_index(drop=True, inplace=True)
    # Converting it to list first as directly passing the df column converted to float with some NaN values
    # (5 NaN values for chromosome 121)
    post_list = maf_df['pos'].tolist()
    maf_df.sort_values('pos', inplace=True)
    maf_df
    # maf_df.reset_index(drop=True, inplace=True)

    # GET POPOOLATION DF
    snp_df = pd.read_csv(genotype_file, sep=" ", header=None)
    snp_df.insert(loc = 0,
                column = 'temp_pos',
                value = post_list
    )
    snp_df.sort_values('temp_pos', inplace=True)
    snp_df.drop(columns='temp_pos', inplace = True)
    # Resetting index to not remove the wrong rows in the next step
    snp_df.reset_index(drop=True, inplace=True)
    
    
    # GET BAYPASS DF
    baypass_df = pd.read_csv(baypass_file, sep=" ", header=None)
    baypass_df.insert(loc = 0,
                column = 'temp_pos',
                value = post_list
    )
    baypass_df.sort_values('temp_pos', inplace=True)
    baypass_df.drop(columns='temp_pos', inplace = True)
    baypass_df.reset_index(drop=True, inplace=True)

    #---------------------------------------------------------------
    print("Calculating population sums...")
    #Switch applymap with map below (FutureWarning)
    population_sums = snp_df.iloc[:, 3:].applymap(
        lambda x: sum(int(num) for num in x.split(':'))
    )
    rows_to_remove = (population_sums <= 54).any(axis=1)
    
    # Remove rows from dataframes
    df_removed_rows = snp_df[rows_to_remove]
    
    snp_df = snp_df[~rows_to_remove]
    maf_df = maf_df[~rows_to_remove]
    baypass_df = baypass_df[~rows_to_remove]
    
    removed_rows_count = df_removed_rows.shape[0]
    print(f"Removed {removed_rows_count} rows from maf, baypass and popoolation files. (Threshold 54)")
    #---------------------------------------------------------------------------------------------
    # OUTPUT
    maf_df.to_csv(f'final_maf_files/{maf_file}', sep=",", index=False)
    print(f"Modified maf data written back to {maf_file}")
    snp_df.to_csv(f'final_popoolation_files/{genotype_file}', sep=" ", header=None, index=False)
    print(f"Modified data written back to {genotype_file}")
    baypass_df.to_csv(f'final_baypass_files/{baypass_file}', sep=" ", header=None, index=False)
    print(f"Modified data written back to {genotype_file}")
    #---------------------------------------------------------------------------------------------
    return removed_rows_count


#----------------------------------------------------------------------------------------------------------------------------

def baypass_results_optimized_remove_rows(maf_path, genotype_path, baypass_result_path):
    
    # GET ORIGINAL MAF DF
    maf_df = pd.read_csv(maf_path, sep=",", header=1)
    maf_df.drop(index = detect_zero_sum_rows(maf_df), inplace = True)
    # Resetting index to not remove the wrong rows in the next step
    maf_df.reset_index(drop=True, inplace=True)
    snp_df = pd.read_csv(genotype_path, sep=" ", header=None)
    # Converting it to list first as directly passing the df column converted to float with some NaN values
    # (5 NaN values for chromosome 121)
    # REMOVE THRESHOLD LINES FROM ORIGINAL MAF DF BASED ON POPOOLATION DATA
    #--------------------------------------------------------------------------------------------
    print("Calculating population sums...")
    #Switch applymap with map below (FutureWarning)
    population_sums = snp_df.iloc[:, 3:].applymap(
        lambda x: sum(int(num) for num in x.split(':'))
    )
    rows_to_remove = (population_sums <= 54).any(axis=1)
    
    # Remove rows from the original maf df
    df_removed_rows = snp_df[rows_to_remove]
    maf_df = maf_df[~rows_to_remove]
    
    removed_rows_count = df_removed_rows.shape[0]
    print(f"Removed {removed_rows_count} rows from original maf file. (Threshold 54)")
    #---------------------------------------------------------------------------------------------

    post_list = maf_df['pos'].tolist()
    baypass_result = pd.read_csv(baypass_result_path, header=None, skiprows=1, delim_whitespace=True, engine='python')
    baypass_result.rename(columns = {0:'MRK', 1:'M_P', 2:'SD_P', 3:'M_XtX', 
                                    4:'SD_XtX', 5:'XtXst', 6:'log10(1/pval)'}, 
                                    inplace=True)
    mrk_list = baypass_result['MRK'].tolist()
    baypass_result.drop(columns='MRK', inplace = True)
    baypass_result.insert(loc = 0,
                column = 'temp_pos',
                value = post_list
    )
    baypass_result.sort_values('temp_pos', inplace=True)
    baypass_result.drop(columns='temp_pos', inplace = True)
    baypass_result.reset_index(drop=True, inplace=True)
    baypass_result.insert(loc = 0,
            column = 'MRK',
            value = mrk_list
    )
    #---------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    # OUTPUT
    baypass_result.to_csv(f'final_baypass_result_files_resorted/{baypass_result_path}', header=0, index=False)
    print(f"Modified data written back to final_baypass_result_files_resorted/{baypass_result_path}")
    #---------------------------------------------------------------------------------------------
    return f'Resorted {baypass_result_path}!!!'

# Calling the above function
for i in range(21,44+1):
    chromosome = rf'LR5371{i}'
    # get the final maf files
    maf_path = rf'maf_{chromosome}.csv'
    genotype_path = rf'maf_{chromosome}_popoolation.genotype'
    baypass_result_path = rf'Sparus_{chromosome}_baypass_summary_pi_xtx.out'
    print(baypass_results_optimized_remove_rows(maf_path, genotype_path, baypass_result_path))
    
    
# TEST--------------

# ###############################################################
chromosome = rf'LR537121'
# get the final maf files
maf_path = rf'maf_{chromosome}.csv'
baypass_result_path = rf'Sparus_{chromosome}_baypass_summary_pi_xtx.out'
genotype_path = rf'maf_{chromosome}_popoolation.genotype'
# print(baypass_results_optimized_remove_rows(maf_path, baypass_result_path))
    
# GET MAF DF
maf_df = pd.read_csv(maf_path, sep=",", header=1)
maf_df
maf_df.drop(index = detect_zero_sum_rows(maf_df), inplace = True)
# Resetting index to not remove the wrong rows in the next step
maf_df.reset_index(drop=True, inplace=True)
snp_df = pd.read_csv(genotype_path, sep=" ", header=None)
maf_df
snp_df
print("Calculating population sums...")
#Switch applymap with map below (FutureWarning)
population_sums = snp_df.iloc[:, 3:].applymap(
    lambda x: sum(int(num) for num in x.split(':'))
)
rows_to_remove = (population_sums <= 54).any(axis=1)
rows_to_remove
# Remove rows from the original maf df
df_removed_rows = snp_df[rows_to_remove]
df_removed_rows
maf_df = maf_df[~rows_to_remove]
maf_df
removed_rows_count = df_removed_rows.shape[0]
print(f"Removed {removed_rows_count} rows from original maf file. (Threshold 54)")


# Converting it to list first as directly passing the df column converted to float with some NaN values
# (5 NaN values for chromosome 121)
post_list = maf_df['pos'].tolist()
baypass_result = pd.read_csv(baypass_result_path, header=None, skiprows=1, delim_whitespace=True, engine='python')
baypass_result
baypass_result.rename(columns = {0:'MRK', 1:'M_P', 2:'SD_P', 3:'M_XtX', 
                                4:'SD_XtX', 5:'XtXst', 6:'log10(1/pval)'}, 
                                inplace=True)
mrk_list = baypass_result['MRK'].tolist()
baypass_result.drop(columns='MRK', inplace = True)
baypass_result.insert(loc = 0,
            column = 'temp_pos',
            value = post_list
)
baypass_result.sort_values('temp_pos', inplace=True)
baypass_result.drop(columns='temp_pos', inplace = True)
baypass_result.reset_index(drop=True, inplace=True)
baypass_result.insert(loc = 0,
            column = 'MRK',
            value = mrk_list
)
baypass_result
###############################################################################################

#----------------------------------------------------------------------------------------------------------------------------
# Use a dictionary comprehension to build the data structure
# Prepare data for DataFrame

for maf_file, genotype_file, baypass_file in zip(maf_files, genotype_files, baypass_files):
    # maf_df = pd.read_csv(maf_file, sep=",", header=1)
    # maf_df.drop(index = detect_zero_sum_rows(maf_df), inplace = True)
    print(optimized_remove_rows(maf_file, genotype_file, baypass_file))






# TESTING FOR 121 CHROMOSOME 

maf_file = 'maf_LR537121.csv'
genotype_file = 'maf_LR537121_popoolation.genotype'
baypass_file = 'maf_LR537121.genotype'

# GET MAF DF
maf_df = pd.read_csv(maf_file, sep=",", header=0)
print(detect_zero_sum_rows(maf_df))
maf_df.drop(index = detect_zero_sum_rows(maf_df), inplace = True)
print(maf_df)
maf_df.reset_index(drop=True, inplace=True)
# Converting it to list first as directly passing the df column converted to float with some NaN values
# (5 NaN values for chromosome 121)
post_list = maf_df['pos'].tolist()
maf_df.sort_values('pos', inplace=True)
maf_df

# GET POPOOLATION DF
snp_df = pd.read_csv(genotype_file, sep=" ", header=None)
snp_df.insert(loc = 0,
            column = 'temp_pos',
            value = post_list
)
snp_df.sort_values('temp_pos', inplace=True)
snp_df.drop(columns='temp_pos', inplace = True)

# GET BAYPASS DF
baypass_df = pd.read_csv('maf_LR537121.genotype', sep=" ", header=None)
baypass_df.insert(loc = 0,
            column = 'temp_pos',
            value = post_list
)
baypass_df.sort_values('temp_pos', inplace=True)
baypass_df.drop(columns='temp_pos', inplace = True)

# Resetting index to not remove the wrong rows in the next stepqqqq

snp_df.reset_index(drop=True, inplace=True)
baypass_df.reset_index(drop=True, inplace=True)



baypass_df
snp_df
# OUTPUT DF TEST


maf_df.to_csv(f'test1.csv', sep=",", index=False)
print(f"Modified maf data written")
snp_df.to_csv(f'test2_pop.txt', sep=" ", header=None, index=False)
print(f"Modified data written")
baypass_df.to_csv(f'test3_BP.genotype', sep=" ", header=None, index=False)
print(f"Modified data written")
