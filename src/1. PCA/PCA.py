'''

PCA
We work with the biallelic SNP file
step 1: note the order of the populations as the order will prove useful to know in the final figure (to know which population is which)
step 2: for the PCA we will need the frequency of one of the two alleles,
otherwise the dataset becomes orthogonal and it doesn't work with the PCA.

'''
import pandas as pd 
import numpy as np
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# import os
# Import the maf csv file

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
                                                                                        desc='Creating dictionary', 
                                                                                        unit='row')
    }
    
    return result_dict




import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# def run_pca(biallelic_fish_df, maf_filename):
#     result_dict = major_allele_freq_calc(biallelic_fish_df)
#     major_alllele_freq_df = pd.DataFrame(result_dict)
    
#     scaler = StandardScaler()
#     scaled_data = scaler.fit_transform(major_alllele_freq_df)
#     scaled_data = pd.DataFrame(scaled_data).replace([np.nan, -np.inf], 0)
    
#     pca = PCA(n_components=2)
#     pca_fin = pca.fit_transform(scaled_data)
    
#     populations = [
#         'wGRE_9', 'wGRE_13', 'wSPA_5', 'fSPA_3', 'fFRA_1', 'fGRE_10', 'wTUR_14', 'wSPA_4', 
#         'wITA_8', 'fSPA_2', 'fGRE_9', 'fGRE_8', 'fCRO_5', 'wITA_7', 'wGRE_12', 'wGRE_11', 
#         'wGRE_10', 'fITA_4', 'fGRE_6', 'fGRE_7'
#     ]

#     finalDf = pd.DataFrame(pca_fin, columns=['PC1', 'PC2'])
#     finalDf.insert(0, 'pops', populations)

#     colors = ['blue' if 'w' in pop else 'red' for pop in populations]

#     fig, ax = plt.subplots(figsize=(12, 8))
#     for label, color in zip(set(finalDf['pops']), set(colors)):
#         indicesToKeep = finalDf['pops'] == label
#         ax.scatter(finalDf.loc[indicesToKeep, 'PC1'], finalDf.loc[indicesToKeep, 'PC2'], c=color, s=50, label=label)
    
#     ax.set_xlabel(f'Principal Component 1', fontsize=15)
#     ax.set_ylabel(f'Principal Component 2', fontsize=15)
#     ax.set_title('Scatter Plot of PCA data')
#     ax.legend(title='Population', frameon=True, facecolor='white', edgecolor='dimgray')
#     ax.grid()

#     output_filename = maf_filename.split('/')[-1].replace('.csv', '_pca_plot.png')
#     fig.savefig(f'v1/{output_filename}')  # Save the figure using fig
#     plt.close(fig)  # Close the figure explicitly
#     print(f'{output_filename} created in v1 folder\n')


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
    scaled_data.replace([np.nan, -np.inf], 0, inplace=True)
    # perform the PCA analysis
    pca = PCA(n_components=2)
    pca.fit(scaled_data)
    pca_fin = pca.transform(scaled_data)
    # get the explained variance
    # print(pca.explained_variance_ratio_)
    populations = [
        'wGRE_9', 'wGRE_13', 'wSPA_5', 'fSPA_3',
        'fFRA_1', 'fGRE_10', 'wTUR_14', 'wSPA_4', 
        'wITA_8', 'fSPA_2', 'fGRE_9', 'fGRE_8', 
        'fCRO_5', 'wITA_7', 'wGRE_12', 'wGRE_11', 
        'wGRE_10', 'fITA_4', 'fGRE_6', 'fGRE_7'
    ]  # list of all 20 population names, they will need to match the input maf file

    wild_color = '#FFFF00'  # Yellow
    farmed_color = '#008000'  # Green
    colors = [
        wild_color, wild_color, wild_color, farmed_color,
        farmed_color, farmed_color, wild_color, wild_color,
        wild_color, farmed_color, farmed_color, farmed_color,
        farmed_color, wild_color, wild_color, wild_color,
        wild_color, farmed_color, farmed_color, farmed_color
    ]
    # colors = [
    #     '#003f5c', '#003f5c', '#003f5c', '#aacc00',  
    #     '#aacc00',  '#aacc00',  '#003f5c', '#003f5c',
    #     '#003f5c', '#aacc00',  '#aacc00',  '#aacc00',  
    #     '#aacc00', '#003f5c', '#003f5c', '#003f5c', 
    #     '#003f5c', '#aacc00',  '#aacc00',  '#aacc00'
    # ] # Color assignments for 'populations' list
    
    scaled_data=pd.DataFrame(scaled_data)
    scaled_data.insert(0, 'pops', populations)
    pca_fin=pd.DataFrame(pca_fin)
    finalDf = pd.concat([pca_fin, scaled_data[['pops']]], axis =1)
    print(finalDf.head)
    #------------------------------------------
    # PCA Plotting
    # Set global style and background color
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = 'white'
    # fig = plt.figure(figsize = (8,8))
    # ax = fig.add_subplot(1,1,1) 
    # ax.set_xlabel('Principal Component 1', fontsize = 15)
    # ax.set_ylabel('Principal Component 2', fontsize = 15)
    # ax.set_title('2 component PCA', fontsize = 20)

    # for labels, color in zip(populations, colors):
    #     indicesToKeep = finalDf['pops'] == populations
    #     ax.scatter(finalDf.iloc[:,0]
    #             , finalDf.iloc[:,1]
    #             , c = colors
    #             , s = 50)
    
    # ax.legend(['Farmed', 'Wild'], frameon=True, facecolor='white', edgecolor='dimgray')
    # ax.grid()
    #------------------------------------------
    # plt.figure(figsize=(12,8))
    # plt.scatter(finalDf.iloc[:,0],finalDf.iloc[:,1],c=colors, s=120)
    # # below previously had  (12,1% explained variation) added
    # plt.xlabel('Principal Component 1', fontsize=20, labelpad=15)
    # # below previously had  (7,3% explained variation) added, might need to calculate this value dynamically
    # plt.ylabel('Principal Component 2', fontsize=20)
    # plt.title('Scatter Plot of PCA data')
    # # plt.legend(['Farmed', 'Wild'], frameon=True, facecolor='white', edgecolor='dimgray')
    # # (Optional) - Add labels to each point
    # for i, txt in enumerate(finalDf.iloc[:,2]):
    #     plt.annotate(txt, (finalDf.iloc[i, 0], finalDf.iloc[i, 1]), fontsize=8)

    # # Save the plot with a name derived from the input file
    # output_filename = maf_filename.split('/')[-1].replace('.csv', '_pca_plot.png')
    # plt.savefig(f'{output_filename}')
    # plt.close()  # Close the plot to free memory
    # print(f'{output_filename} created in folder\n')

    #---------------------------------------

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(10, 10))

    # Hide grid lines
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Plotting
    ax.set_title('Scatter Plot of PCA Data', fontsize=24)
    ax.set_xlabel('Principal Component 1', fontsize=18)
    ax.set_ylabel('Principal Component 2', fontsize=18)
    
    scatter = ax.scatter(finalDf.iloc[:, 0], finalDf.iloc[:, 1], c=colors, edgecolors='dimgray', s=150)
    
    legend_handles = [Patch(facecolor=farmed_color, edgecolor=farmed_color, label='Farmed'),
                        Patch(facecolor=wild_color, edgecolor=wild_color, label='Wild')]
    # legend_handles = [Patch(facecolor='#aacc00', edgecolor='#aacc00', label='Farmed'),
    #                     Patch(facecolor='#003f5c', edgecolor='#003f5c', label='Wild')]
    ax.legend(handles=legend_handles, frameon=False, facecolor='white', edgecolor='dimgray', fontsize=14)
    # add padding to the title and axes labels
    # ax.title.set_pad(20)
    # ax.xaxis.label.set_label_coords(0.5, -0.1)
    # ax.yaxis.label.set_label_coords(-0.1, 0.5)
    
    plt.show()

    
    # ax.legend(['Farmed', 'Wild'], frameon=True, facecolor='white', edgecolor='dimgray')
    # ax.legend(handles=scatter.legend_elements()[0], labels=['Farmed', 'Wild'])



    # Optional - Add labels to each point
    # for i, txt in enumerate(finalDf.iloc[:,2]):
    #     ax.annotate(txt, (finalDf.iloc[i, 0], finalDf.iloc[i, 1]), fontsize=8)

    # Save the plot
    output_filename = maf_filename.split('/')[-1].replace('.csv', '_pca_plot.png')
    fig.savefig(output_filename)
    plt.close(fig)  # Close the figure to free memory
    print(f'{output_filename} created in folder')
    # todo - create sample data and configure plot to them to finalize it





# Call PCA here or in the main function
chr_list = range(21, 45)

for CHR in chr_list:
    maf_filename = rf'../../data/processed_maf/maf_LR5371{CHR}.csv'
    biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
    try:
        print(f"Processing {maf_filename} with PCA...")
        run_pca(biallelic_fish_df, maf_filename)
        print(f"Successfully processed {maf_filename} with PCA.\n")
    except Exception as e:
        print(f"Error processing {maf_filename} with PCA: {e}\n")
        
        
    
    run_pca(biallelic_fish_df, maf_filename)

# merge all maf files:
merged_maf = pd.DataFrame()
for CHR in chr_list:
    maf_filename = rf'../../data/processed_maf/maf_LR5371{CHR}.csv'
    biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
    merged_maf = pd.concat([merged_maf, biallelic_fish_df], ignore_index=True)
# Save the merged DataFrame to a new CSV file
merged_maf.to_csv('../../data/processed_maf/maf_merged.csv', index=False)




# Now run PCA on the merged file
maf_filename = rf'../../data/processed_maf/maf_all.csv'
biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
try:
    print(f"Processing {maf_filename} with PCA...")
    run_pca(biallelic_fish_df, maf_filename)
    print(f"Successfully processed {maf_filename} with PCA.\n")
except Exception as e:
    print(f"Error processing {maf_filename} with PCA: {e}\n")
    
    

run_pca(biallelic_fish_df, maf_filename)




# Now run PCA on the merged file
maf_filename = rf'../../data/raw_maf/maf_all.csv'
biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
try:
    print(f"Processing {maf_filename} with PCA...")
    run_pca(biallelic_fish_df, maf_filename)
    print(f"Successfully processed {maf_filename} with PCA.\n")
except Exception as e:
    print(f"Error processing {maf_filename} with PCA: {e}\n")
    
    

run_pca(biallelic_fish_df, maf_filename)
















#---------------------------------------------------------------------
# PLOTLY VISUALIZATION
#---------------------------------------------------------------------


import pandas as pd 
import numpy as np
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.express as px
import plotly.graph_objects as go


# The function estimates the major allele frequency of the FIRST population 
# then estimates the allele frequency of the SAME allele in the rest of the populations
def major_allele_freq_calc(input_dataframe):
    # Convert DataFrame to numpy array for faster operations
    genetic_data = input_dataframe.iloc[:, 4:].to_numpy() 
    
    # Calculate the number of samples (assuming 4 columns per sample)
    num_samples = genetic_data.shape[1] // 4
        
    # Group alleles for each sample
    reshaped_data = genetic_data.reshape(genetic_data.shape[0], num_samples, 4)
    
    # Find the index of the maximum value in the first sample
    max_indices = reshaped_data[:, 0, :].argmax(axis=1)
    
    # Calculate major allele frequencies
    major_allele_freqs = np.zeros((reshaped_data.shape[0], num_samples))
    for i in range(num_samples):
        sample_data = reshaped_data[:, i, :]
        major_alleles = sample_data[np.arange(sample_data.shape[0]), max_indices]
        allele_sums = sample_data.sum(axis=1)
        major_allele_freqs[:, i] = np.round(major_alleles / allele_sums, 3) 
    
    result_dict = {input_dataframe.iloc[i, 1]: list(major_allele_freqs[i]) for i in tqdm(range(genetic_data.shape[0]), 
                                                                                        desc='Creating dictionary', 
                                                                                        unit='row')
    }
    
    return result_dict


def run_pca(biallelic_fish_df, maf_filename):
    # Calculate allele frequencies
    result_dict = major_allele_freq_calc(biallelic_fish_df)
    major_allele_freq_df = pd.DataFrame(result_dict)
    
    # Standardize the data
    scaler = StandardScaler()
    scaler.fit(major_allele_freq_df)
    scaled_data = pd.DataFrame(scaler.transform(major_allele_freq_df))
    
    # Replace 'nan' and 'inf' values with 0
    scaled_data.replace([np.nan, -np.inf], 0, inplace=True)
    
    # Perform PCA analysis
    pca = PCA(n_components=2)
    pca.fit(scaled_data)
    pca_result = pca.transform(scaled_data)
    
    # Get explained variance for axis labels
    explained_variance = pca.explained_variance_ratio_ * 100
    
    # Population names
    populations = [
        'wGRE_9', 'wGRE_13', 'wSPA_5', 'fSPA_3',
        'fFRA_1', 'fGRE_10', 'wTUR_14', 'wSPA_4', 
        'wITA_8', 'fSPA_2', 'fGRE_9', 'fGRE_8', 
        'fCRO_5', 'wITA_7', 'wGRE_12', 'wGRE_11', 
        'wGRE_10', 'fITA_4', 'fGRE_6', 'fGRE_7'
    ]
    
    # Create category labels (wild vs farmed)
    categories = ['Wild' if pop.startswith('w') else 'Farmed' for pop in populations]
    
    # Create DataFrame for plotting
    pca_df = pd.DataFrame({
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1],
        'Population': populations,
        'Category': categories
    })
    
    # Create Plotly figure
    fig = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Category',
        hover_name='Population',
        color_discrete_map={'Wild': '#FFFF00', 'Farmed': '#008000'},
        labels={
            'PC1': f'Principal Component 1 ({explained_variance[0]:.1f}%)',
            'PC2': f'Principal Component 2 ({explained_variance[1]:.1f}%)'
        },
        title='PCA Analysis of Fish Populations',
        size_max=15,
        template='plotly_white'
    )
    
    # Update marker size and add borders
    fig.update_traces(
        marker=dict(
            size=15,
            line=dict(width=1, color='DarkSlateGrey')
        )
    )
    
    # Customize layout
    fig.update_layout(
        title_font_size=24,
        legend_title_font_size=14,
        legend_font_size=12,
        xaxis_title_font_size=16,
        yaxis_title_font_size=16,
        width=900,
        height=700
    )
    
    # Save as HTML and PNG
    output_filename_base = maf_filename.split('/')[-1].replace('.csv', '_pca')
    fig.write_html(f'{output_filename_base}.html')
    fig.write_image(f'{output_filename_base}plotly.png')
    
    # Display the plot
    fig.show()
    
    print(f"PCA plots saved as {output_filename_base}.html and {output_filename_base}.png")
    
    return pca_df


# Process individual chromosomes
def process_chromosomes(chr_list):
    for CHR in chr_list:
        maf_filename = f'../../data/processed_maf/maf_LR5371{CHR}.csv'
        try:
            print(f"Processing {maf_filename} with PCA...")
            biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
            run_pca(biallelic_fish_df, maf_filename)
            print(f"Successfully processed {maf_filename} with PCA.\n")
        except Exception as e:
            print(f"Error processing {maf_filename} with PCA: {e}\n")


# Process merged data
def process_merged_data():
    # Merge all chromosome files
    print("Creating merged file from all chromosomes...")
    chr_list = range(21, 45)
    merged_maf = pd.DataFrame()
    
    for CHR in chr_list:
        maf_filename = f'../../data/processed_maf/maf_LR5371{CHR}.csv'
        try:
            biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
            merged_maf = pd.concat([merged_maf, biallelic_fish_df], ignore_index=True)
        except Exception as e:
            print(f"Error reading {maf_filename}: {e}")
    
    # Save the merged DataFrame
    merged_output = '../../data/processed_maf/maf_merged.csv'
    merged_maf.to_csv(merged_output, index=False)
    print(f"Merged file saved as {merged_output}")
    
    # Run PCA on merged file
    try:
        print("Processing merged file with PCA...")
        run_pca(merged_maf, merged_output)
        print("Successfully processed merged file with PCA.\n")
    except Exception as e:
        print(f"Error processing merged file with PCA: {e}\n")


# Process additional files
def process_additional_files():
    additional_files = [
        '../../data/processed_maf/maf_all.csv',
        '../../data/raw_maf/maf_all.csv'
    ]
    
    for filename in additional_files:
        try:
            print(f"Processing {filename} with PCA...")
            biallelic_fish_df = pd.read_csv(filename, sep=",", header=1)
            run_pca(biallelic_fish_df, filename)
            print(f"Successfully processed {filename} with PCA.\n")
        except Exception as e:
            print(f"Error processing {filename} with PCA: {e}\n")


# Main execution
if __name__ == "__main__":
    # Process chromosomes 21-44
    chr_list = range(21, 45)
    process_chromosomes(chr_list)
    
    # Process merged data
    process_merged_data()
    
    # Process additional files
    process_additional_files()









#----------------------------------------------------------------------------------------------------



import pandas as pd 
import numpy as np
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# The function estimates the major allele frequency of the FIRST population 
# then estimates the allele frequency of the SAME allele in the rest of the populations
def major_allele_freq_calc(input_dataframe):
    # Convert DataFrame to numpy array for faster operations
    genetic_data = input_dataframe.iloc[:, 4:].to_numpy() 
    
    # Calculate the number of samples (assuming 4 columns per sample)
    num_samples = genetic_data.shape[1] // 4
        
    # Group alleles for each sample
    reshaped_data = genetic_data.reshape(genetic_data.shape[0], num_samples, 4)
    
    # Find the index of the maximum value in the first sample
    max_indices = reshaped_data[:, 0, :].argmax(axis=1)
    
    # Calculate major allele frequencies
    major_allele_freqs = np.zeros((reshaped_data.shape[0], num_samples))
    for i in range(num_samples):
        sample_data = reshaped_data[:, i, :]
        major_alleles = sample_data[np.arange(sample_data.shape[0]), max_indices]
        allele_sums = sample_data.sum(axis=1)
        major_allele_freqs[:, i] = np.round(major_alleles / allele_sums, 3) 
    
    result_dict = {input_dataframe.iloc[i, 1]: list(major_allele_freqs[i]) for i in tqdm(range(genetic_data.shape[0]), 
                                                                                        desc='Creating dictionary', 
                                                                                        unit='row')
    }
    
    return result_dict


def run_pca(biallelic_fish_df, maf_filename):
    # Calculate allele frequencies
    result_dict = major_allele_freq_calc(biallelic_fish_df)
    major_allele_freq_df = pd.DataFrame(result_dict)
    
    # Standardize the data
    scaler = StandardScaler()
    scaler.fit(major_allele_freq_df)
    scaled_data = pd.DataFrame(scaler.transform(major_allele_freq_df))
    
    # Replace 'nan' and 'inf' values with 0
    scaled_data.replace([np.nan, -np.inf], 0, inplace=True)
    
    # Perform PCA analysis - get 3 components for 3D visualization
    pca = PCA(n_components=3)
    pca.fit(scaled_data)
    pca_result = pca.transform(scaled_data)
    
    # Get explained variance for axis labels
    explained_variance = pca.explained_variance_ratio_ * 100
    
    # Population names
    populations = [
        'wGRE_9', 'wGRE_13', 'wSPA_5', 'fSPA_3',
        'fFRA_1', 'fGRE_10', 'wTUR_14', 'wSPA_4', 
        'wITA_8', 'fSPA_2', 'fGRE_9', 'fGRE_8', 
        'fCRO_5', 'wITA_7', 'wGRE_12', 'wGRE_11', 
        'wGRE_10', 'fITA_4', 'fGRE_6', 'fGRE_7'
    ]
    
    # Extract country codes for additional grouping
    countries = [pop.split('_')[0][1:] for pop in populations]  # Remove first letter (w/f) and get country code
    
    # Create category labels (wild vs farmed)
    categories = ['Wild' if pop.startswith('w') else 'Farmed' for pop in populations]
    
    # Create DataFrame for plotting with all metadata
    pca_df = pd.DataFrame({
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1],
        'PC3': pca_result[:, 2],
        'Population': populations,
        'Category': categories,
        'Country': countries,
        'ID': [pop.split('_')[1] for pop in populations]
    })
    
    # Create 2D and 3D visualizations in a single HTML file
    create_interactive_plots(pca_df, explained_variance, maf_filename)
    
    return pca_df


def create_interactive_plots(pca_df, explained_variance, maf_filename):
    # Define color schemes
    color_map = {'Wild': '#FFFF00', 'Farmed': '#008000'}
    
    # Create figures
    fig_2d = create_2d_plot(pca_df, explained_variance, color_map)
    fig_3d = create_3d_plot(pca_df, explained_variance, color_map)
    fig_country = create_country_plot(pca_df, explained_variance)
    fig_scree = create_scree_plot(explained_variance)
    
    # Create dashboard with tabs
    output_filename_base = maf_filename.split('/')[-1].replace('.csv', '_pca')
    
    # Create a dashboard HTML with all plots
    dashboard = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "scatter"}, {"type": "scatter3d"}],
               [{"type": "scatter"}, {"type": "scatter"}]],
        subplot_titles=("2D PCA by Category", "3D PCA Visualization", 
                        "PCA by Country", "Scree Plot"),
        vertical_spacing=0.1,
        horizontal_spacing=0.05
    )
    
    # Add 2D plot traces
    for trace in fig_2d.data:
        dashboard.add_trace(trace, row=1, col=1)
    
    # Add 3D plot traces
    for trace in fig_3d.data:
        dashboard.add_trace(trace, row=1, col=2)
    
    # Add country plot traces
    for trace in fig_country.data:
        dashboard.add_trace(trace, row=2, col=1)
    
    # Add scree plot trace
    dashboard.add_trace(fig_scree.data[0], row=2, col=2)
    
    # Update layout
    dashboard.update_layout(
        title=f"Interactive PCA Analysis Dashboard - {output_filename_base}",
        height=1000,
        width=1200,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    # Save individual plots and dashboard
    fig_2d.write_html(f'{output_filename_base}_2d.html')
    fig_3d.write_html(f'{output_filename_base}_3d.html')
    fig_country.write_html(f'{output_filename_base}_country.html')
    dashboard.write_html(f'{output_filename_base}_dashboard.html')
    
    # Save as static image for reports
    fig_2d.write_image(f'{output_filename_base}_2d.png')
    
    print(f"Interactive plots saved as {output_filename_base}_dashboard.html")
    print(f"Individual plots also available as separate HTML files")


def create_2d_plot(pca_df, explained_variance, color_map):
    """Create the 2D PCA plot with additional interactivity"""
    fig = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Category',
        symbol='Category',  # Use different symbols for wild vs farmed
        hover_name='Population',
        hover_data=['Country', 'ID'],  # Show these fields on hover
        color_discrete_map=color_map,
        labels={
            'PC1': f'Principal Component 1 ({explained_variance[0]:.1f}%)',
            'PC2': f'Principal Component 2 ({explained_variance[1]:.1f}%)'
        },
        title='PCA Analysis of Fish Populations',
        size_max=15,
        template='plotly_white'
    )
    
    # Update marker size and add borders
    fig.update_traces(
        marker=dict(
            size=15,
            line=dict(width=1, color='DarkSlateGrey')
        )
    )
    
    # Add annotations for each point
    for i, row in pca_df.iterrows():
        fig.add_annotation(
            x=row['PC1'],
            y=row['PC2'],
            text=row['Population'],
            showarrow=False,
            font=dict(size=8),
            xshift=10,
            yshift=10,
            visible=False,  # Hidden by default, will be shown via button
            hovertext=f"Country: {row['Country']}<br>ID: {row['ID']}"
        )
    
    # Add buttons for interactive features
    fig.update_layout(
        updatemenus=[
            # Button to show/hide labels
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[{"annotations[0:20].visible": True}],
                        label="Show Labels",
                        method="relayout"
                    ),
                    dict(
                        args=[{"annotations[0:20].visible": False}],
                        label="Hide Labels",
                        method="relayout"
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=1.1,
                yanchor="top"
            ),
            # Button to toggle marker size by explained variance
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[{"marker.size": 15}],
                        label="Reset Size",
                        method="update"
                    ),
                    dict(
                        args=[{"marker.size": [10 + 40 * abs(val) for val in pca_df['PC1']]}],
                        label="Size by PC1",
                        method="update"
                    ),
                    dict(
                        args=[{"marker.size": [10 + 40 * abs(val) for val in pca_df['PC2']]}],
                        label="Size by PC2",
                        method="update"
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.5,
                xanchor="left",
                y=1.1,
                yanchor="top"
            )
        ]
    )
    
    # Add shape outlines for groups (wild/farmed)
    for category in pca_df['Category'].unique():
        category_data = pca_df[pca_df['Category'] == category]
        
        # Calculate convex hull (simplified approach - just use min/max)
        x_min, x_max = category_data['PC1'].min() - 0.1, category_data['PC1'].max() + 0.1
        y_min, y_max = category_data['PC2'].min() - 0.1, category_data['PC2'].max() + 0.1
        
        color = color_map[category]
        
        # Add a semi-transparent rectangle
        fig.add_shape(
            type="rect",
            x0=x_min, y0=y_min,
            x1=x_max, y1=y_max,
            line=dict(color=color, width=1),
            fillcolor=color,
            opacity=0.1,
            layer="below",
            visible="legendonly"  # Toggle with legend
        )
    
    # Customize layout
    fig.update_layout(
        title_font_size=24,
        legend_title_font_size=14,
        legend_font_size=12,
        xaxis_title_font_size=16,
        yaxis_title_font_size=16,
        width=900,
        height=700
    )
    
    return fig


def create_3d_plot(pca_df, explained_variance, color_map):
    """Create a 3D PCA plot"""
    fig = px.scatter_3d(
        pca_df, 
        x='PC1', 
        y='PC2', 
        z='PC3',
        color='Category',
        hover_name='Population',
        hover_data=['Country', 'ID'],
        color_discrete_map=color_map,
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.1f}%)',
            'PC2': f'PC2 ({explained_variance[1]:.1f}%)',
            'PC3': f'PC3 ({explained_variance[2]:.1f}%)'
        },
        title='3D PCA Visualization'
    )
    
    # Update marker size
    fig.update_traces(
        marker=dict(
            size=8,
            line=dict(width=0.5, color='DarkSlateGrey')
        )
    )
    
    # Update layout
    fig.update_layout(
        scene=dict(
            xaxis_title=f'PC1 ({explained_variance[0]:.1f}%)',
            yaxis_title=f'PC2 ({explained_variance[1]:.1f}%)',
            zaxis_title=f'PC3 ({explained_variance[2]:.1f}%)',
        ),
        width=800,
        height=700,
        margin=dict(l=0, r=0, b=0, t=30)
    )
    
    return fig


def create_country_plot(pca_df, explained_variance):
    """Create a PCA plot colored by country"""
    # Define a color map for countries (add more as needed)
    country_colors = {
        'GRE': '#1f77b4',  # Blue
        'SPA': '#ff7f0e',  # Orange
        'ITA': '#2ca02c',  # Green
        'FRA': '#d62728',  # Red
        'TUR': '#9467bd',  # Purple
        'CRO': '#8c564b',  # Brown
    }
    
    fig = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Country',
        symbol='Category',  # Use different symbols for wild vs farmed
        hover_name='Population',
        color_discrete_map=country_colors,
        labels={
            'PC1': f'Principal Component 1 ({explained_variance[0]:.1f}%)',
            'PC2': f'Principal Component 2 ({explained_variance[1]:.1f}%)'
        },
        title='PCA by Country of Origin',
        template='plotly_white'
    )
    
    # Update marker properties
    fig.update_traces(
        marker=dict(
            size=15,
            line=dict(width=1, color='DarkSlateGrey')
        )
    )
    
    # Add categorical markers for wild vs farmed
    for category in ['Wild', 'Farmed']:
        category_data = pca_df[pca_df['Category'] == category]
        fig.add_traces(
            go.Scatter(
                x=category_data['PC1'],
                y=category_data['PC2'],
                mode='markers',
                marker=dict(
                    symbol='circle' if category == 'Wild' else 'square',
                    size=30,
                    opacity=0.3,
                    line=dict(width=2, color='black')
                ),
                name=f"{category} Outline",
                showlegend=True,
                hoverinfo='skip'
            )
        )
    
    return fig


def create_scree_plot(explained_variance):
    """Create a scree plot showing explained variance"""
    # Add data for all possible components
    components = list(range(1, len(explained_variance) + 1))
    
    # Create cumulative variance
    cumulative_variance = np.cumsum(explained_variance)
    
    fig = go.Figure()
    
    # Add bar chart for individual explained variance
    fig.add_trace(
        go.Bar(
            x=components,
            y=explained_variance,
            name='Explained Variance',
            marker_color='royalblue'
        )
    )
    
    # Add line chart for cumulative explained variance
    fig.add_trace(
        go.Scatter(
            x=components,
            y=cumulative_variance,
            name='Cumulative Variance',
            marker_color='red',
            yaxis='y2'
        )
    )
    
    # Update layout with second y-axis
    fig.update_layout(
        title='Scree Plot - Explained Variance by Principal Component',
        xaxis_title='Principal Component',
        yaxis_title='Explained Variance (%)',
        yaxis2=dict(
            title='Cumulative Variance (%)',
            overlaying='y',
            side='right',
            range=[0, 100]
        ),
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='right',
            x=1
        ),
        template='plotly_white'
    )
    
    return fig


# Process individual chromosomes
def process_chromosomes(chr_list):
    for CHR in chr_list:
        maf_filename = f'../../data/processed_maf/maf_LR5371{CHR}.csv'
        try:
            print(f"Processing {maf_filename} with PCA...")
            biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
            run_pca(biallelic_fish_df, maf_filename)
            print(f"Successfully processed {maf_filename} with PCA.\n")
        except Exception as e:
            print(f"Error processing {maf_filename} with PCA: {e}\n")


# Process merged data
def process_merged_data():
    # Merge all chromosome files
    print("Creating merged file from all chromosomes...")
    chr_list = range(21, 45)
    merged_maf = pd.DataFrame()
    
    for CHR in chr_list:
        maf_filename = f'../../data/processed_maf/maf_LR5371{CHR}.csv'
        try:
            biallelic_fish_df = pd.read_csv(maf_filename, sep=",", header=1)
            merged_maf = pd.concat([merged_maf, biallelic_fish_df], ignore_index=True)
        except Exception as e:
            print(f"Error reading {maf_filename}: {e}")
    
    # Save the merged DataFrame
    merged_output = '../../data/processed_maf/maf_merged.csv'
    merged_maf.to_csv(merged_output, index=False)
    print(f"Merged file saved as {merged_output}")
    
    # Run PCA on merged file
    try:
        print("Processing merged file with PCA...")
        run_pca(merged_maf, merged_output)
        print("Successfully processed merged file with PCA.\n")
    except Exception as e:
        print(f"Error processing merged file with PCA: {e}\n")


# Process additional files
def process_additional_files():
    additional_files = [
        # '../../data/processed_maf/maf_all.csv',
        '../../data/raw_maf/maf_all.csv'
    ]
    
    for filename in additional_files:
        try:
            print(f"Processing {filename} with PCA...")
            biallelic_fish_df = pd.read_csv(filename, sep=",", header=1)
            run_pca(biallelic_fish_df, filename)
            print(f"Successfully processed {filename} with PCA.\n")
        except Exception as e:
            print(f"Error processing {filename} with PCA: {e}\n")


# Main execution
if __name__ == "__main__":
    # Process chromosomes 21-44
    # chr_list = [21]
    # process_chromosomes(chr_list)
    
    # Process merged data
    # process_merged_data()
    
    # Process additional files
    process_additional_files()
