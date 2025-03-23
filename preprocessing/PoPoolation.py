# PoPoolation preparation

import pandas as pd
import numpy as np
from tqdm import tqdm

def run_popoolation_prep(biallelic_fish_maf, maf_filename):
    # Drop the first 4 columns, leaving only the nucleotide columns in a new dataframe

    nuc_data = biallelic_fish_maf.iloc[:, 4:]

    # Create a list according to the order of populations in the biallelic maf file

    populations = [
        'wGRE_9', 'wGRE_13', 'wSPA_5', 'fSPA_3',
        'fFRA_1', 'fGRE_10', 'wTUR_14', 'wSPA_4', 
        'wITA_8', 'fSPA_2', 'fGRE_9', 'fGRE_8', 
        'fCRO_5', 'wITA_7', 'wGRE_12', 'wGRE_11', 
        'wGRE_10', 'fITA_4', 'fGRE_6', 'fGRE_7'
    ]  # list of all 20 population names

    def process_population_data(biallelic_fish_maf, populations):
        pop_data = {}
        for i, pop in tqdm(enumerate(populations), total=len(populations), desc="Processing Populations", unit='population'):
            base_col = i * 4
            A = biallelic_fish_maf.iloc[:, base_col].astype(str)
            T = biallelic_fish_maf.iloc[:, base_col+3].astype(str)
            C = biallelic_fish_maf.iloc[:, base_col+1].astype(str)
            G = biallelic_fish_maf.iloc[:, base_col+2].astype(str)
            
            pop_data[pop] = A + ':' + T + ':' + C + ':' + G + ':0:0'
        return pop_data


    # Usage
    pop_data = process_population_data(nuc_data, populations)

    # We create a dataframe based on the pop_data dictionary with the ordered_populations order

    ordered_populations = [
        'fSPA_3', 'fFRA_1', 'fGRE_10', 'fSPA_2', 
        'fGRE_9', 'fGRE_8', 'fCRO_5', 'fITA_4', 
        'fGRE_6', 'fGRE_7', 'wGRE_9', 'wGRE_13',
        'wSPA_5', 'wTUR_14', 'wSPA_4', 'wITA_8', 
        'wITA_7', 'wGRE_12', 'wGRE_11', 'wGRE_10'
    ] # all "farmed" populations first for ease of use. 

    sync_fish = pd.DataFrame.from_dict(pop_data)[ordered_populations]

    # Add the rest of the columns to the dataframe

    # Extract relevant columns from the source DataFrame
    chr_column = biallelic_fish_maf.iloc[:, 0]  # Chromosome info
    pos_column = biallelic_fish_maf.iloc[:, 1]  # Position info
    ref_n_column = biallelic_fish_maf.iloc[:, 2]  # Reference nucleotide info

    # Insert "chr", "pos", and "ref_n" columns into the synchronized file
    sync_fish.insert(0, "chr", chr_column)
    sync_fish.insert(1, "pos", pos_column)
    sync_fish.insert(2, "ref_n", ref_n_column)

    # File export
    # save the synchronized file to run the analyses
    output_filename = maf_filename.split('/')[-1].replace('.csv', '_popoolation.genotype')
    sync_fish.to_csv(output_filename, header=False, index=False, sep=" ")
