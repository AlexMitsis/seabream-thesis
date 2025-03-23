# Importing the scripts as part of a package
from preprocessing import BayPass, PCA, PoPoolation
import pandas as pd
import numpy as np
# Library for Logging
# import loguru
# Library for better markdown console messages
# from rich.console import Console
# from rich.markdown import Markdown
# Library for cli
# import argparse
# Library alternative to pandas for quicker dataframe operations (supports multithreading as well)
# import polars
# data validation 
# import pydantic

maf_files = ['maf_LR537121.csv', 'maf_LR537122.csv', 'maf_LR537123.csv', 'maf_LR537124.csv', 
             'maf_LR537125.csv', 'maf_LR537126.csv', 'maf_LR537127.csv', 'maf_LR537128.csv', 
             'maf_LR537129.csv', 'maf_LR537130.csv', 'maf_LR537131.csv', 'maf_LR537132.csv', 
             'maf_LR537133.csv', 'maf_LR537134.csv', 'maf_LR537135.csv', 'maf_LR537136.csv', 
             'maf_LR537137.csv', 'maf_LR537138.csv', 'maf_LR537139.csv', 'maf_LR537140.csv', 
             'maf_LR537141.csv', 'maf_LR537142.csv', 'maf_LR537143.csv', 'maf_LR537144.csv', 'maf_all.csv']

# maf_files = ['maf_all.csv']

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
    zero_count = np.sum(allele_sums == 0) 

    print(f"Number of rows with zero allele sum: {zero_sum_positions.shape[0]}")
    print(f"Total ATCG zero sum instances count: {zero_count}")
    
    return zero_sum_indices

for maf_filename in maf_files:
    # Read the CSV, skipping the duplicate header row
    biallelic_fish_df = pd.read_csv(f'maf_files/{maf_filename}', sep=",", header=1)
    biallelic_fish_df.drop(index = detect_zero_sum_rows(biallelic_fish_df), inplace = True)
    # Processing with BayPass
    try:
        print(f"Processing {maf_filename} with BayPass...")
        BayPass.run_baypass_prep(biallelic_fish_df, maf_filename)
        print(f"Successfully processed {maf_filename} with BayPass.\n")
    except Exception as e:
        print(f"Error processing {maf_filename} with BayPass: {e}\n")
    
    # Processing with PCA
    try:
        print(f"Processing {maf_filename} with PCA...")
        PCA.run_pca(biallelic_fish_df, maf_filename)
        print(f"Successfully processed {maf_filename} with PCA.\n")
    except Exception as e:
        print(f"Error processing {maf_filename} with PCA: {e}\n")
    
    # Processing with PoPoolation
    try:
        print(f"Processing {maf_filename} with PoPoolation...")
        PoPoolation.run_popoolation_prep(biallelic_fish_df, maf_filename)
        print(f"Successfully processed {maf_filename} with PoPoolation.\n")
    except Exception as e:
        print(f"Error processing {maf_filename} with PoPoolation: {e}\n")
