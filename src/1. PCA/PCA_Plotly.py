
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
