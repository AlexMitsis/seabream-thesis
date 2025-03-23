# # BayPass preparation
import pandas as pd 
import numpy as np
from tqdm import tqdm


def BayPass_genotyping_file(infile, colnames):
    print("Starting BayPass genotyping file processing")
    
    # Convert input to numpy array for faster processing
    data = infile.iloc[:, 4:].values
    
    num_rows, num_cols = data.shape

    # Ensure the number of columns is divisible by 4
    # if num_cols % 4 != 0:
    #     raise ValueError(f"Number of genotype columns ({num_cols}) is not divisible by 4.")
    
    # Reshape data to group every 4 columns
    data_reshaped = data.reshape(num_rows, -1, 4)
    num_groups = data_reshaped.shape[1]

    print(f"Input shape: {data.shape}, Reshaped: {data_reshaped.shape}")
    
    # Initialize result array
    result = np.zeros((num_rows, num_groups * 2), dtype=data.dtype)
        
    for i in tqdm(range(num_rows), desc="Processing rows", unit="row"):
        try:
            row_data = data_reshaped[i]
            non_zero_indices = np.argwhere(row_data != 0)
            first_non_zero = non_zero_indices[0]
            second_non_zero = non_zero_indices[1]
            result[i, ::2] = row_data[first_non_zero[0], first_non_zero[1]]
            result[i, 1::2] = row_data[second_non_zero[0], second_non_zero[1]]
        except Exception as e:
            print(f"Error processing row {i}: {str(e)}")
    
    return pd.DataFrame(result, columns=colnames)


def run_baypass_prep(biallelic_fish_maf, maf_filename):

    # Usage remains the same
    column_names = ['wGRE_9_1', 'wGRE_9_2', 'wGRE_13_1', 'wGRE_13_2',
                    'wSPA_5_1', 'wSPA_5_2', 'fSPA_3_1', 'fSPA_3_2',
                    'fFRA_1_1', 'fFRA_1_2', 'fGRE_10_1', 'fGRE_10_2',
                    'wTUR_14_1', 'wTUR_14_2', 'wSPA_4_1', 'wSPA_4_2',
                    'wITA_8_1', 'wITA_8_2', 'fSPA_2_1', 'fSPA_2_2',
                    'fGRE_9_1', 'fGRE_9_2', 'fGRE_8_1', 'fGRE_8_2',
                    'fCRO_5_1', 'fCRO_5_2', 'wITA_7_1', 'wITA_7_2',
                    'wGRE_12_1', 'wGRE_12_2', 'wGRE_11_1', 'wGRE_11_2',
                    'wGRE_10_1', 'wGRE_10_2', 'fITA_4_1', 'fITA_4_2',
                    'fGRE_6_1', 'fGRE_6_2', 'fGRE_7_1', 'fGRE_7_2']


    try:
        fish_maf_baypass = BayPass_genotyping_file(biallelic_fish_maf, column_names)
    except Exception as e:
        print(f"Error in main execution: {str(e)}")

    # BayPass analysis
    # save without header and index, with space separation
    output_filename = maf_filename.split('/')[-1].replace('.csv', '.genotype')
    fish_maf_baypass.to_csv(output_filename, sep=' ', index=False, header=False)
