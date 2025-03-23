import pandas as pd
import statsmodels.stats.multitest as stats
import math
from tqdm import tqdm
import numpy as np

# BayPass_sparus_result = pd.read_csv('merged_baypass_results.csv',  sep=',')
# BayPass_sparus_result

# # TEST BELOW
# chromosome_num = 29
# chromosome = rf'LR5371{chromosome_num}'

# input_baypass_path = rf'D:/final_baypass_results_Jan_2025/Final_Sparus_{chromosome}_baypass_summary_pi_xtx.out'
# input_baypass_C2_path = rf'D:/final_baypass_results_Jan_2025/Final_Sparus_{chromosome}_baypass_summary_contrast.out'
# input_maf_path = rf'C:/Users/mitsi/Desktop/seabream repo/data/processed_maf/maf_{chromosome}.csv'
# output_path = f'Sparus_{chromosome}_baypass_result.csv'

# input_baypass_C2_path
# BayPass_C2 = pd.read_csv(input_baypass_C2_path,  header=0, delim_whitespace=True)
# BayPass_C2
# BayPass_C2.describe()

# BayPass_sparus_result = pd.read_csv(input_baypass_path,  header=0, delim_whitespace=True)
# BayPass_sparus_result
# # find the max from the log10(1/pval) column
# BayPass_sparus_result["log10(1/pval)"].max()
# # BayPass_sparus_result["log10(1/pval)"].head(sortby='log10(1/pval)')
# BayPass_sparus_result.sort_values(by=['log10(1/pval)']).tail(60)

# # Filter out infinites first, then take max
# max_decimal = BayPass_sparus_result['log10(1/pval)'][~BayPass_sparus_result['log10(1/pval)'].isin([np.inf, -np.inf])].max()
# # Then replace infinites with that max + 1
# BayPass_sparus_result['log10(1/pval)'].replace([np.inf], max_decimal + 1, inplace=True)

# # BayPass_sparus_result.rename(columns = {0:'MRK', 1:'M_P', 2:'SD_P', 3:'M_XtX', 4:'SD_XtX', 5:'XtXst', 6:'log10(1/pval)'}, inplace=True)
# Sparus_maf= pd.read_csv(input_maf_path)
# Sparus_maf
# pos = Sparus_maf['pos']
# BayPass_sparus_result['pos'] = pos

def filter_and_align_chromosomes(baypass_df, popoolation_df, print_stats=True):
    """
    Filter and align chromosome data between Baypass and Popoolation datasets.
    
    Parameters:
    -----------
    baypass_df : pandas.DataFrame
        Baypass dataset with 'CHR' and 'BP' columns
    popoolation_df : pandas.DataFrame
        Popoolation dataset with 'CHR' and 'BP' columns
    print_stats : bool, optional
        Whether to print detailed statistics (default: True)
        
    Returns:
    --------
    tuple:
        - filtered_baypass_df : pandas.DataFrame
            Filtered Baypass dataset with matching positions
        - filtered_popoolation_df : pandas.DataFrame
            Filtered Popoolation dataset with matching positions
    """
    # Split dataframes by chromosome
    baypass_by_chr = dict(tuple(baypass_df.groupby('CHR')))
    popoolation_by_chr = dict(tuple(popoolation_df.groupby('CHR')))

    # Create mapping between chromosome formats (LR537121 -> 21)
    chr_mapping = {chr_id: int(chr_id.replace('LR5371', '')) for chr_id in baypass_by_chr.keys()}

    filtered_baypass = []
    filtered_popoolation = []

    if print_stats:
        print("\nDetailed Statistics per Chromosome:")
        print("-" * 80)
        print(f"{'Baypass Chr':<12} {'Pop Chr':<8} {'Baypass Rows':<12} "
            f"{'Pop Rows':<12} {'Kept':<8} {'Removed':<8} {'% Kept':<8}")
        print("-" * 80)

    total_original_baypass = 0
    total_original_pop = 0
    total_kept = 0

    # Process each chromosome
    for bay_chr, bay_data in baypass_by_chr.items():
        pop_chr = chr_mapping[bay_chr]  # Get corresponding popoolation chromosome number
        
        # Skip if chromosome not in popoolation data
        if pop_chr not in popoolation_by_chr:
            if print_stats:
                print(f"{bay_chr:<12} {pop_chr:<8} {'N/A':<12} {'N/A':<12} "
                    f"{'N/A':<8} {'N/A':<8} {'N/A':<8}")
            continue
        
        pop_data = popoolation_by_chr[pop_chr]
        
        # Find common positions
        common_positions = set(bay_data['BP']) & set(pop_data['BP'])
        
        # Filter both dataframes
        bay_filtered = bay_data[bay_data['BP'].isin(common_positions)]
        pop_filtered = pop_data[pop_data['BP'].isin(common_positions)]
        
        # Add to results
        filtered_baypass.append(bay_filtered)
        filtered_popoolation.append(pop_filtered)
        
        # Calculate statistics
        original_bay = len(bay_data)
        original_pop = len(pop_data)
        kept = len(common_positions)
        removed = original_bay - kept
        percent_kept = (kept / original_bay * 100) if original_bay > 0 else 0
        
        # Update totals
        total_original_baypass += original_bay
        total_original_pop += original_pop
        total_kept += kept
        
        if print_stats:
            print(f"{bay_chr:<12} {pop_chr:<8} {original_bay:<12} "
                f"{original_pop:<12} {kept:<8} {removed:<8} {percent_kept:>6.1f}%")

    if print_stats:
        # Print totals
        print("-" * 80)
        total_removed = total_original_baypass - total_kept
        total_percent_kept = (total_kept / total_original_baypass * 100) \
            if total_original_baypass > 0 else 0
        print(f"{'TOTAL':<12} {'':<8} {total_original_baypass:<12} "
            f"{total_original_pop:<12} {total_kept:<8} {total_removed:<8} "
            f"{total_percent_kept:>6.1f}%")

    # Combine filtered data
    filtered_baypass_df = pd.concat(filtered_baypass)
    filtered_popoolation_df = pd.concat(filtered_popoolation)

    # Reset indices
    filtered_baypass_df.reset_index(drop=True, inplace=True)
    filtered_popoolation_df.reset_index(drop=True, inplace=True)

    if print_stats:
        print(f"\nFinal filtered dataframes shape:")
        print(f"Baypass: {filtered_baypass_df.shape}")
        print(f"Popoolation: {filtered_popoolation_df.shape}")

    return filtered_baypass_df, filtered_popoolation_df






popoolation_df = pd.read_csv('C:/Users/mitsi/Desktop/seabream repo/src/4. postprocessing/remove extra rows from bp/Manhattan_plot_Sp_aurata_Popoolation_dataset.csv')
popoolation_by_chr = dict(tuple(popoolation_df.groupby('CHR')))
# Check the shape of all dataframes
for chromosome_num in tqdm(range(21,45), desc="Processing Chromosomes", unit="chromosome", colour="blue"):
    chromosome=rf'LR5371{chromosome_num}'
    # Change paths to the ones in your system
    input_baypass_path = rf'D:/final_baypass_results_Jan_2025/Final_Sparus_{chromosome}_baypass_summary_pi_xtx.out'
    input_baypass_C2_path = rf'D:/final_baypass_results_Jan_2025/Final_Sparus_{chromosome}_baypass_summary_contrast.out'
    input_maf_path = rf'C:/Users/mitsi/Desktop/seabream repo/data/processed_maf/maf_{chromosome}.csv'
    Sparus_maf= pd.read_csv(input_maf_path)
    BayPass_sparus_result = pd.read_csv(input_baypass_path,  header=0, delim_whitespace=True)
    BayPass_C2 = pd.read_csv(input_baypass_C2_path,  header=0, delim_whitespace=True)
    sync_Sparus_fst = pd.read_csv(rf'D:/final_popoolation_results_resorted/sync_Sp_aurata_{chromosome}_popoolation.fst', sep="\t", header=None)
    # sync_Sparus_fet = pd.read_csv(rf'D:/final_popoolation_results_resorted/sync_Sp_aurata_{chromosome}_popoolation.fet',  sep="\t", header=None)
    # TODO - Get positions / BP from all files and compare
    fst_pos_col = sync_Sparus_fst[1]
    # fet_pos_col = sync_Sparus_fet[1]
    
    # baypass_pos_col = Sparus_maf['pos']
    # c2_pos_col = Sparus_maf['pos']
    # maf_pos_col = Sparus_maf['pos']
    
    BayPass_sparus_result['BP'] = Sparus_maf['pos']
    # compare row count for baypass, c2, maf files and popoolation files
    print(f'Baypass: {BayPass_sparus_result.shape}')
    print(f'C2: {BayPass_C2.shape}')
    print(f'MAF: {Sparus_maf.shape}')
    print(f'Popoolation fst: {sync_Sparus_fst.shape}')
    print(f'Popoolation fet: {sync_Sparus_fst.shape}')
    
    filtered_baypass, filtered_popoolation = filter_and_align_chromosomes(
        BayPass_sparus_result, 
        sync_Sparus_fst,
        print_stats=True  # Set to False to suppress statistics output
    )
    # Drop the pos column from the Baypass dataframe
    filtered_baypass = filtered_baypass.drop(columns=['BP'])
    # filtered_baypass, filtered_popoolation = filter_and_align_chromosomes(
    # BayPass_sparus_result, 
    # popoolation_df,
    # print_stats=True  # Set to False to suppress statistics output
    # )
#---------------------
# Constants
# TODO - Read PoPoolation file and get the positions
for chromosome_num in tqdm(range(21,45), desc="Processing Chromosomes", unit="chromosome", colour="blue"):
    chromosome=rf'LR5371{chromosome_num}'
    # Change paths to the ones in your system
    input_baypass_path = rf'D:/final_baypass_results_Jan_2025/Final_Sparus_{chromosome}_baypass_summary_pi_xtx.out'
    input_baypass_C2_path = rf'D:/final_baypass_results_Jan_2025/Final_Sparus_{chromosome}_baypass_summary_contrast.out'
    input_maf_path = rf'C:/Users/mitsi/Desktop/seabream repo/data/processed_maf/maf_{chromosome}.csv'
    output_path = f'../../../data/plots/baypass/Sparus_{chromosome}_baypass_result.csv'
    C2_output_path = f'../../../data/plots/baypass/C2_Sparus_{chromosome}_baypass_result.csv'

    BayPass_sparus_result = pd.read_csv(input_baypass_path,  header=0, delim_whitespace=True)
    BayPass_C2 = pd.read_csv(input_baypass_C2_path,  header=0, delim_whitespace=True)
    
    # check on "inf" log10(1/pval) - they are 50 values for gilthead seabream (36 values for labrax in previous study)
    # emailed the author of BayPass (Mathieu Gautier)
    # answer: if you get inf, it means that p was smaller than machine precision (or rahter 1/p higher than machine precision),
    # so you can safely transform such log10(1/pval) values to the highest observed value +1 (per chromosome
    
    # Filter out infinites first, then take max
    xtx_max = BayPass_sparus_result['log10(1/pval)'][~BayPass_sparus_result['log10(1/pval)'].isin([np.inf, -np.inf])].max()
    # Then replace infinites with that max + 1
    BayPass_sparus_result['log10(1/pval)'].replace([np.inf], xtx_max + 1, inplace=True)
    print(f'max xtx value for chr {chromosome}: {xtx_max}')
    
    # Filter out infinites first, then take max
    c2_max = BayPass_C2['log10(1/pval)'][~BayPass_C2['log10(1/pval)'].isin([np.inf, -np.inf])].max()
    # Then replace infinites with that max + 1
    BayPass_C2['log10(1/pval)'].replace([np.inf], c2_max + 1, inplace=True)
    print(f'max C2 value for chr {chromosome}: {c2_max}')
    
    BayPass_sparus_result["delog_pval"] = 10 ** - BayPass_sparus_result["log10(1/pval)"]
    BayPass_C2["delog_pval"] = 10 ** - BayPass_C2["log10(1/pval)"]
    BayPass_sparus_result_BH = stats.multipletests(list(BayPass_sparus_result['delog_pval']), alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    BayPass_C2_BH = stats.multipletests(list(BayPass_C2['delog_pval']), alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    
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
            
    # do the same for the contrast file (C2)
    non_zero_values = []
    for i in list(BayPass_C2_BH[1]):
        if i > 0:
            non_zero_values.append(i)
    minimum_non_zero_value = min(non_zero_values)
    
    log_transformed_values_C2 = []
    for i in list(BayPass_C2_BH[1]):
        try:
            log_transformed_values_C2.append(-math.log10(i))
        except ValueError:
            log_transformed_values_C2.append(-math.log10(minimum_non_zero_value))
    
    # append the log-transformed values as column to the df of interest
    BayPass_sparus_result["log10(1/pval)_BH"] = log_transformed_values
    BayPass_C2["log10(1/pval)_BH"] = log_transformed_values_C2
    BayPass_sparus_result['chromosome'] = str(chromosome)
    BayPass_C2['chromosome'] = str(chromosome)
    # add the pos column from the maf file to the baypass file
    Sparus_maf= pd.read_csv(input_maf_path)
    pos = Sparus_maf['pos']
    BayPass_sparus_result['pos'] = pos
    BayPass_C2['pos'] = pos
    BayPass_sparus_result.to_csv(output_path, sep=",", index=False)
    BayPass_C2.to_csv(C2_output_path, sep=",", index=False)
    
#------------------------------------------------------------------------------------ Merging
# merge final baypass files

# Initialize an empty list to hold dataframes
df_list = []
C2_df_list = []
for chromosome_num in range(21, 44 + 1):
    chromosome=rf'LR5371{chromosome_num}'
    # Format the filename
    output_path = f'../../../data/plots/baypass/Sparus_{chromosome}_baypass_result.csv'
    C2_output_path = f'../../../data/plots/baypass/C2_Sparus_{chromosome}_baypass_result.csv'
    # Read the CSV file
    temp_df = pd.read_csv(output_path)
    temp_df_C2 = pd.read_csv(C2_output_path)
    # temp_df.sort_values()
    # Append the dataframe to the list
    df_list.append(temp_df)
    C2_df_list.append(temp_df_C2)

# Concatenate all dataframes in the list
merged_df = pd.concat(df_list, ignore_index=True)
C2_merged_df = pd.concat(C2_df_list, ignore_index=True)


popoolation_df = pd.read_csv('../../../data/plots/popoolation/Manhattan_plot_Sp_aurata_Popoolation_dataset.csv')
merged_df.to_csv("../../../data/plots/baypass/merged_baypass_results.csv", index=False)
C2_merged_df.to_csv("../../../data/plots/baypass/C2_merged_baypass_results.csv", index=False)
print("Merging completed. Merged file saved as 'merged_baypass_results.csv'.")

popoolation_df
filtered_baypass_df = pd.read_csv("../../../data/plots/baypass/merged_baypass_results.csv")
C2_filtered_baypass_df = pd.read_csv("../../../data/plots/baypass/C2_merged_baypass_results.csv")

filtered_baypass_df
C2_filtered_baypass_df
merged_df


test = pd.read_csv('../../../PoPoolation/final_prep/Manhattan_plot_Sp_aurata_Popoolation_dataset.csv')
test
#------------------------------------------------------------------------------------Testing
import pandas as pd
import statsmodels.stats.multitest as stats
import math
from tqdm import tqdm
import numpy as np

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
BayPass_sparus_result_BH.head()
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







# TESTING - 14/2/2025
# Filter the Baypass and Popoolation data so that it matches before we proceed with the rest of the analysis

