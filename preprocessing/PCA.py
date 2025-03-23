# # PCA

# We work with the biallelic SNP file
# 
# step 1: note the order of the populations as the order will prove useful to know in the final figure (to know which population is which)
# 
# step 2: for the PCA we will need the frequency of one of the two alleles,
# otherwise the dataset becomes orthogonal and it doesn't work with the PCA.

import pandas as pd 
import numpy as np
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# import os
# ## Import the maf csv file

# The function "major_allele_freq_calc" below
# estimates the major allele frequency of the FIRST population rounded to third decimal for every pop - position
# then it estimates the allele frequency of the SAME allele in the rest of the populations
# pca will be based on the allele frequency of that same allele across populations
def major_allele_freq_calc(input_dataframe):
    # Convert DataFrame to numpy array for faster operations (around 300x speed benefit, 00:03 vs 09:14)
    # We convert the dataframe but without the first 4 columns
    genetic_data = input_dataframe.iloc[:, 4:].to_numpy() 
    
    # Calculate the number of samples (assuming 4 columns per sample)
    num_samples = genetic_data.shape[1] // 4
        
    # Group alleles for each sample. Reorganize the data into a 3D array: (number of SNPs, number of samples, 4 alleles per sample)
    reshaped_data = genetic_data.reshape(genetic_data.shape[0], num_samples, 4)
    
    # Find the index of the maximum value in the first sample. Identify the major allele based on the first sample
    max_indices = reshaped_data[:, 0, :].argmax(axis=1)
    
    # Calculate major allele frequencies
    major_allele_freqs = np.zeros((reshaped_data.shape[0], num_samples))
    for i in range(num_samples):
        sample_data = reshaped_data[:, i, :]
        major_alleles = sample_data[np.arange(sample_data.shape[0]), max_indices]
        allele_sums = sample_data.sum(axis=1)
        major_allele_freqs[:, i] = np.round(major_alleles / allele_sums, 3) 
    
    
    result_dict = {input_dataframe.iloc[i, 1]: list(major_allele_freqs[i]) for i in tqdm(range(genetic_data.shape[0]), 
                                                                                        desc="Creating dictionary", 
                                                                                        unit="row")
    }
    
    return result_dict

def run_pca(biallelic_fish_df, maf_filename):
    
    result_dict = major_allele_freq_calc(biallelic_fish_df)
    
    major_alllele_freq_df = pd.DataFrame(result_dict)
    # conduct the PCA analysis using the PCA function of the sklearn package
    # 
    # standard scaler standardizes the variance of each row to unit and scales the data to zero
    # 
    # this way makes possible for PCA to work as unequal variances ineterfere with the PCA results
    scaler = StandardScaler()
    scaler.fit(major_alllele_freq_df)
    scaled_data=pd.DataFrame(scaler.transform(major_alllele_freq_df))
    # Replace 'nan' and 'inf' values with 0 (check if this is needed)
    scaled_data = scaled_data.replace([np.nan, -np.inf], 0)
    # perform the PCA analysis
    pca = PCA(n_components=2)
    pca.fit(scaled_data)
    pca_fin = pca.transform(scaled_data)
    # get the explained variance
    print(pca.explained_variance_ratio_)
    populations = [
        'wGRE_9', 'wGRE_13', 'wSPA_5', 'fSPA_3',
        'fFRA_1', 'fGRE_10', 'wTUR_14', 'wSPA_4', 
        'wITA_8', 'fSPA_2', 'fGRE_9', 'fGRE_8', 
        'fCRO_5', 'wITA_7', 'wGRE_12', 'wGRE_11', 
        'wGRE_10', 'fITA_4', 'fGRE_6', 'fGRE_7'
    ]  # list of all 20 population names, they will need to match the input maf file

    scaled_data=pd.DataFrame(scaled_data)
    scaled_data.insert(0, "pops", populations)
    pca_fin=pd.DataFrame(pca_fin)
    finalDf = pd.concat([pca_fin, scaled_data[['pops']]], axis =1)
    
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)

    populations = [
        'wGRE_9', 'wGRE_13', 'wSPA_5', 'fSPA_3',
        'fFRA_1', 'fGRE_10', 'wTUR_14', 'wSPA_4', 
        'wITA_8', 'fSPA_2', 'fGRE_9', 'fGRE_8', 
        'fCRO_5', 'wITA_7', 'wGRE_12', 'wGRE_11', 
        'wGRE_10', 'fITA_4', 'fGRE_6', 'fGRE_7'
    ]  # list of all 20 population names

    # Color assignments for 'populations' list
    colors = [
        'blue', 'blue', 'blue', 'red',  
        'red',  'red',  'blue', 'blue',
        'blue', 'red',  'red',  'red',  
        'red', 'blue', 'blue', 'blue', 
        'blue', 'red',  'red',  'red'
    ]

    for labels, color in zip(populations, colors):
        indicesToKeep = finalDf['pops'] == populations
        ax.scatter(finalDf.iloc[:,0]
                , finalDf.iloc[:,1]
                , c = colors
                , s = 50)
    # ax.legend(populations)
    ax.grid()

    plt.figure(figsize=(12,8))
    plt.scatter(finalDf.iloc[:,0],finalDf.iloc[:,1],c=colors, s=200)
    plt.xlabel('First principal component (12,1% explained variation)', fontsize=22, labelpad=15)
    plt.ylabel('Second Principal Component (7,3% explained variation)', fontsize=14)

    for i, txt in enumerate(finalDf.iloc[:,2]):
        plt.annotate(txt, (finalDf.iloc[i, 0], finalDf.iloc[i, 1]), fontsize=8)

    # Save the plot with a name derived from the input file
    output_filename = maf_filename.split('/')[-1].replace('.csv', '_pca_plot.png')
    plt.savefig(output_filename)
    plt.close()  # Close the plot to free memory
    print(f"{output_filename} created\n")


