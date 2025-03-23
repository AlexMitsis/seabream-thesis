import pandas as pd
import numpy as np
import statsmodels.stats.multitest as stats
import math
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import math
from tqdm import tqdm
from pathlib import Path

def filter_and_align_chromosomes(base_dir="../../../data/plots", chromosome_range=(21, 45), print_stats=True):
    """
    Filter and align chromosome data between Baypass and Popoolation datasets.
    
    Parameters:
    -----------
    base_dir : str
        Base directory containing the data files
    chromosome_range : tuple
        Range of chromosomes to process (start, end+1)
    print_stats : bool, optional
        Whether to print detailed statistics (default: True)
        
    Returns:
    --------
    tuple:
        - filtered_baypass : pandas.DataFrame
            Merged filtered Baypass dataset
        - filtered_popoolation : pandas.DataFrame
            Merged filtered Popoolation dataset
    """
    base_dir = Path(base_dir)
    filtered_baypass_list = []
    filtered_popoolation_list = []
    
    if print_stats:
        print("\nDetailed Statistics per Chromosome:")
        print("-" * 80)
        print(f"{'Baypass Chr':<12} {'Pop Chr':<8} {'Baypass Rows':<12} "
              f"{'Pop Rows':<12} {'Kept':<8} {'Removed':<8} {'% Kept':<8}")
        print("-" * 80)

    total_original_baypass = 0
    total_original_pop = 0
    total_kept = 0
    
    try:
        # Load Popoolation data
        popoolation_df = pd.read_csv(base_dir / "popoolation/Manhattan_plot_Sp_aurata_Popoolation_dataset.csv")
        
        # Process each chromosome
        for chromosome_num in range(*chromosome_range):
            chromosome = f'LR5371{chromosome_num}'
            baypass_path = base_dir / "baypass" / f'Sparus_{chromosome}_baypass_result.csv'
            
            # Read individual chromosome results
            baypass_chr = pd.read_csv(baypass_path)
            
            # Create mapping between chromosome formats (LR537121 -> 21)
            pop_chr = int(chromosome.replace('LR5371', ''))
            
            # Get popoolation data for this chromosome
            pop_data = popoolation_df[popoolation_df['CHR'] == pop_chr]
            
            # Find common positions
            common_positions = set(baypass_chr['BP']) & set(pop_data['BP'])
            
            # Filter both dataframes
            bay_filtered = baypass_chr[baypass_chr['BP'].isin(common_positions)]
            pop_filtered = pop_data[pop_data['BP'].isin(common_positions)]
            
            # Calculate statistics
            original_bay = len(baypass_chr)
            original_pop = len(pop_data)
            kept = len(common_positions)
            removed = original_bay - kept
            percent_kept = (kept / original_bay * 100) if original_bay > 0 else 0
            
            # Update totals
            total_original_baypass += original_bay
            total_original_pop += original_pop
            total_kept += kept
            
            if print_stats:
                print(f"{chromosome:<12} {pop_chr:<8} {original_bay:<12} "
                      f"{original_pop:<12} {kept:<8} {removed:<8} {percent_kept:>6.1f}%")
            
            filtered_baypass_list.append(bay_filtered)
            filtered_popoolation_list.append(pop_filtered)
            
            # Save individual filtered results
            bay_filtered.to_csv(base_dir / "baypass" / f'filtered_Sparus_{chromosome}_baypass_result.csv', 
                              index=False)
        
        if print_stats:
            print("-" * 80)
            total_removed = total_original_baypass - total_kept
            total_percent_kept = (total_kept / total_original_baypass * 100) \
                if total_original_baypass > 0 else 0
            print(f"{'TOTAL':<12} {'':<8} {total_original_baypass:<12} "
                  f"{total_original_pop:<12} {total_kept:<8} {total_removed:<8} "
                  f"{total_percent_kept:>6.1f}%")
        
        # Merge filtered results
        filtered_baypass = pd.concat(filtered_baypass_list, ignore_index=True)
        filtered_popoolation = pd.concat(filtered_popoolation_list, ignore_index=True)
        
        # Save merged filtered results
        filtered_baypass.to_csv(base_dir / "baypass/filtered_merged_baypass_results.csv", index=False)
        filtered_popoolation.to_csv(base_dir / "popoolation/filtered_popoolation_results.csv", index=False)
        
        return filtered_baypass, filtered_popoolation
        
    except Exception as e:
        print(f"Error during filtering: {str(e)}")
        raise
    # Split dataframes by chromosome
    baypass_by_chr = dict(tuple(baypass_df.groupby('CHR')))
    popoolation_by_chr = dict(tuple(popoolation_df.groupby('CHR')))

    # Create mapping between chromosome formats (LR537121 -> 21)
    chr_mapping = {chr_id: int(chr_id.replace('LR5371', '')) 
                  for chr_id in baypass_by_chr.keys()}

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
        pop_chr = chr_mapping[bay_chr]
        
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

def _process_infinite_values(df, column):
    """Helper function to process infinite values in a dataframe column."""
    max_val = df[column][~df[column].isin([np.inf, -np.inf])].max()
    df[column].replace([np.inf], max_val + 1, inplace=True)
    return max_val

def _transform_pvalues(pvalues):
    """Helper function to transform p-values with log10 transformation."""
    non_zero_values = [x for x in pvalues if x > 0]
    minimum_non_zero = min(non_zero_values)
    
    transformed_values = []
    for pval in pvalues:
        try:
            transformed_values.append(-math.log10(pval))
        except ValueError:
            transformed_values.append(-math.log10(minimum_non_zero))
    
    return transformed_values

def process_baypass_results(chromosome_range=(21, 45), 
                          base_baypass_dir="D:/final_baypass_results_Jan_2025",
                          maf_dir="C:/Users/mitsi/Desktop/seabream repo/data/processed_maf",
                          output_dir="../../../data/plots/baypass",
                          alpha=0.05):
    """
    Process BayPass results for multiple chromosomes, handling both standard and contrast (C2) analyses.
    
    Parameters:
    -----------
    chromosome_range : tuple
        Range of chromosomes to process (start, end+1)
    base_baypass_dir : str
        Base directory containing BayPass result files
    maf_dir : str
        Directory containing MAF files
    output_dir : str
        Directory for output files
    alpha : float
        Significance level for multiple testing correction
        
    Returns:
    --------
    dict:
        Dictionary containing processed results for each chromosome
    """
    results = {}
    
    for chromosome_num in tqdm(range(*chromosome_range), desc="Processing Chromosomes", 
                             unit="chromosome", colour="blue"):
        chromosome = f'LR5371{chromosome_num}'
        
        # Define file paths
        paths = {
            'baypass': f'{base_baypass_dir}/Final_Sparus_{chromosome}_baypass_summary_pi_xtx.out',
            'baypass_c2': f'{base_baypass_dir}/Final_Sparus_{chromosome}_baypass_summary_contrast.out',
            'maf': f'{maf_dir}/maf_{chromosome}.csv',
            'output': f'{output_dir}/Sparus_{chromosome}_baypass_result.csv',
            'output_c2': f'{output_dir}/C2_Sparus_{chromosome}_baypass_result.csv'
        }
        
        # Read input files
        baypass_result = pd.read_csv(paths['baypass'], header=0, delim_whitespace=True)
        baypass_c2 = pd.read_csv(paths['baypass_c2'], header=0, delim_whitespace=True)
        
        # Process XTX values
        xtx_max = _process_infinite_values(baypass_result, 'log10(1/pval)')
        c2_max = _process_infinite_values(baypass_c2, 'log10(1/pval)')
        print(f'max xtx value for chr {chromosome}: {xtx_max}')
        print(f'max C2 value for chr {chromosome}: {c2_max}')
        
        # Calculate delogged p-values
        for df in [baypass_result, baypass_c2]:
            df["delog_pval"] = 10 ** -df["log10(1/pval)"]
        
        # Perform multiple testing correction
        baypass_bh = stats.multipletests(baypass_result['delog_pval'], 
                                       alpha=alpha, method='fdr_bh', 
                                       is_sorted=False, returnsorted=False)
        baypass_c2_bh = stats.multipletests(baypass_c2['delog_pval'], 
                                          alpha=alpha, method='fdr_bh', 
                                          is_sorted=False, returnsorted=False)
        
        # Transform BH-corrected p-values
        baypass_result["log10(1/pval)_BH"] = _transform_pvalues(baypass_bh[1])
        baypass_c2["log10(1/pval)_BH"] = _transform_pvalues(baypass_c2_bh[1])
        
        # Add chromosome and position information
        maf_data = pd.read_csv(paths['maf'])
        for df in [baypass_result, baypass_c2]:
            df['chromosome'] = str(chromosome)
            df['pos'] = maf_data['pos']
        
        # Save results
        baypass_result.to_csv(paths['output'], sep=",", index=False)
        baypass_c2.to_csv(paths['output_c2'], sep=",", index=False)
        
        # Store results in dictionary
        results[chromosome] = {
            'baypass': baypass_result,
            'baypass_c2': baypass_c2
        }
    
    return results

def merge_baypass_results(base_dir="../../../data/plots/baypass", 
                         chromosome_range=(21, 45)):
    """
    Merge BayPass results from multiple chromosomes into single files.
    
    Parameters:
    -----------
    base_dir : str
        Base directory containing individual chromosome results
    chromosome_range : tuple
        Range of chromosomes to process (start, end+1)
        
    Returns:
    --------
    tuple:
        merged_df, C2_merged_df : DataFrames containing merged results
    """
    df_list = []
    C2_df_list = []
    
    for chromosome_num in range(*chromosome_range):
        chromosome = f'LR5371{chromosome_num}'
        
        # Define file paths
        output_path = f'{base_dir}/Sparus_{chromosome}_baypass_result.csv'
        C2_output_path = f'{base_dir}/C2_Sparus_{chromosome}_baypass_result.csv'
        
        # Read and append files
        df_list.append(pd.read_csv(output_path))
        C2_df_list.append(pd.read_csv(C2_output_path))
    
    # Merge all results
    merged_df = pd.concat(df_list, ignore_index=True)
    C2_merged_df = pd.concat(C2_df_list, ignore_index=True)
    
    # Save merged results
    merged_df.to_csv(f"{base_dir}/merged_baypass_results.csv", index=False)
    C2_merged_df.to_csv(f"{base_dir}/C2_merged_baypass_results.csv", index=False)
    
    print("Merging completed. Merged files saved successfully.")
    return merged_df, C2_merged_df

def main():
    """
    Main function to orchestrate the entire BayPass analysis workflow.
    """
    # Define common parameters
    base_dir = Path("../../../data/plots")
    baypass_dir = "D:/final_baypass_results_Jan_2025"
    maf_dir = "C:/Users/mitsi/Desktop/seabream repo/data/processed_maf"
    
    try:
        # First step: Filter and align chromosomes
        filtered_baypass, filtered_popoolation = filter_and_align_chromosomes(
            base_dir=str(base_dir),
            chromosome_range=(21, 45),
            print_stats=True
        )
        
        # Second step: Process BayPass results
        results = process_baypass_results(
            chromosome_range=(21, 45),
            base_baypass_dir=baypass_dir,
            maf_dir=maf_dir,
            output_dir=str(base_dir / "baypass"),
            alpha=0.05
        )
        
        print("Analysis workflow completed successfully.")
        return results, filtered_baypass, filtered_popoolation
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        raise
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        raise

if __name__ == "__main__":
    results, merged_df, C2_merged_df, popoolation_df = main()
