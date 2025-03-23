import pandas as pd
import warnings
import numpy as np
import os
import shutil
warnings.filterwarnings('ignore')  # not to print any warnings
pd.set_option('display.max_rows', 30)  # Change 10 to your desired maximum number of rows
pd.set_option('display.max_columns', 25)

maf_files = ['maf_LR537121.csv']
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

# TEST CODE ---------------------------------------------------------------------------------
# TODO: Create function that counts sums for x:x:x:x:0:0 in all files
def count_pattern_sums(maf_files):
    pattern_count = {}
    for file in maf_files:
        df = pd.read_csv(file)
        # Assuming the relevant data starts from the 4th column
        # Modify the range according to your dataframe structure
        for col in df.columns[3:]:
            # Count rows matching the pattern x:x:x:x:0:0
            count = df[col].apply(lambda x: all(int(num) == 0 for num in x.split(':')[4:])).sum()
            pattern_count[col] = pattern_count.get(col, 0) + count
    return pattern_count

# Apply the function to the MAF files
pattern_counts = count_pattern_sums(maf_files)
print(pattern_counts)
#--------------------------------------------------------------------------------------

# # global variables
# INFILE = "maf_LR537121_popoolation - Copy.genotype"
# SUM_THRESHOLD = 54

def optimized_remove_rows(genotype_file, baypass_file, max_threshold):
    snp_df = pd.read_csv(genotype_file, sep=" ", header=None)
    baypass_df = pd.read_csv(baypass_file, sep=" ", header=None)
    # Helper function to calculate sum from population data
    def calculate_population_sum(population_data):
        """Calculate the sum of integers in the 'int:int' formatted string."""
        return sum(int(value) for value in population_data.split(':'))

    # Apply the helper function to each element in the dataframe starting from the 4th column
    # and calculate the sum for each population
    # population_sums = snp_df.iloc[:, 3:].applymap(calculate_population_sum)
    # Refactor to use vectorized operations for sum calculations
    # Split each string into a list of integers, convert to a NumPy array, and sum the array
    print("Calculating population sums...")
    # population_sums = snp_df.iloc[:, 3:].apply(lambda col: col.str.split(':').apply(lambda x: np.array(x, dtype=int).sum()))
    population_sums = snp_df.iloc[:, 3:].applymap(
        lambda x: sum(int(num) for num in x.split(':'))
    )
    # Dictionary to hold counts of removed rows for each threshold
    thresholds_removed = {}
    removed_rows_count_cummulative = 0
    # Apply each threshold incrementally
    for sum_threshold in range(1, max_threshold + 1):
        print(f'threshold: {sum_threshold}')
        rows_to_remove = (population_sums <= sum_threshold).any(axis=1)
        df_removed_rows = snp_df[rows_to_remove]
        snp_df = snp_df[~rows_to_remove]
        removed_rows_count = df_removed_rows.shape[0]
        removed_rows_count_cummulative += removed_rows_count
        # Store the count of removed rows for this threshold
        thresholds_removed[f"Threshold {sum_threshold}"] = removed_rows_count_cummulative
        # Update population_sums for rows not removed
        population_sums = population_sums.loc[~rows_to_remove]
        print(f"Threshold {sum_threshold}: Removed {removed_rows_count} rows cumulatively.")
    
    #---------------------------------------------------------------------------------------------
    # TODO - Do the same process below but for BayPass
    # Backup the original file before modification
    # backup_genotype_file_path = genotype_file + ".bak"
    # shutil.copy(genotype_file, backup_genotype_file_path)
    # print(f"Backup created at {backup_genotype_file_path}")
    
    snp_df.to_csv(f'{genotype_file}_V2', sep=" ", header=None, index=False)
    print(f"Modified data written back to {genotype_file}")
    baypass_df.to_csv(f'{baypass_file}_V2', sep=" ", header=None, index=False)
    print(f"Modified data written back to {genotype_file}")
    #---------------------------------------------------------------------------------------------
    return thresholds_removed
# Use a dictionary comprehension to build the data structure
# Prepare data for DataFrame
for genotype_file, baypass_file in zip(genotype_files,baypass_files):
    optimized_remove_rows(genotype_file, baypass_file, 54)
    
    
# Initial data exploration to find best threshold:
# -------------------------------------------------
data = {genotype_file: optimized_remove_rows(genotype_file, baypass_file, 54) for genotype_file, baypass_file in zip(genotype_files,baypass_files)}

# Convert dictionary to DataFrame
result_df = pd.DataFrame(data).T  # Transpose to make files as rows
print(result_df)
print(result_df.head())
result_df.to_csv('out.csv', index=True)
df_diff = result_df.diff(axis=1)
df_diff.to_csv('out2.csv', index=True)
# -------------------------------------------------

#--------------------------------------------------------------------------------
df_removed_rows = pd.DataFrame(columns=snp_df.columns)
for index, row in snp_df.iterrows():
    if index % 1000 == 0:
        print(index, end=" ")
        
    # skip first 3 columns
    for population in row[3:]:
        print(sum(int(x) for x in population.split(":")))
        pop_sum = sum(int(x) for x in population.split(":"))
        if pop_sum <= SUM_THRESHOLD:
            # Append the row to df_removed_rows
            df_removed_rows = pd.concat([df_removed_rows, pd.DataFrame([row])], ignore_index=True)
            # Remove the row from snp_df
            snp_df.drop(index, inplace=True)
            # Exit the inner loop
            break


# %%
df_removed_rows.shape

# %%
# threshold 1 -->  genotyped positions removed
# threshold 2 -->  genotyped positions removed
# threshold 3 -->  genotyped positions removed
# threshold 4 -->  genotyped positions removed
# threshold 5 -->  genotyped positions removed
# threshold 6 -->  genotyped positions removed
# threshold 7 -->  genotyped positions removed
# threshold 8 -->  genotyped positions removed
# threshold 9 -->  genotyped positions removed
# threshold 10 -->  genotyped positions removed
# threshold 11 -->  genotyped positions removed
# ...
# threshold 15 -->  genotyped positions removed
# ...
# threshold 20 -->  genotyped positions removed
# ...
# threshold 24 -->  genotyped positions removed
# threshold 25 --> 6551 genotyped positions removed

# %%
# threshold 1 --> 15 genotyped positions removed
# threshold 2 --> 19 genotyped positions removed
# threshold 3 --> 21 genotyped positions removed
# threshold 4 --> 27 genotyped positions removed
# threshold 5 --> 27 genotyped positions removed
# threshold 6 --> 32 genotyped positions removed
# threshold 7 --> 37 genotyped positions removed
# threshold 8 --> 41 genotyped positions removed
# threshold 9 --> 42 genotyped positions removed
# threshold 10 --> 43 genotyped positions removed
# threshold 11 --> 45 genotyped positions removed
# ...
# threshold 15 --> 54 genotyped positions removed
# ...
# threshold 20 --> 70 genotyped positions removed
# ...
# threshold 24 --> 96 genotyped positions removed
# threshold 25 --> 6551 genotyped positions removed

# %%
# save the desired filtered file to an outfile
snp_df.to_csv("sync_Dlabrax_LGx_popoolation_threshold_15.txt", sep=" ", header=None, index=False)


