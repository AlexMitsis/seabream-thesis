# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline
# import allel
import statistics as stat
import statsmodels.stats.multitest as stats
import math

# %%
new_address = r'D:\final_popoolation_results_resorted'

for CHR in range(24, 25):
    sync_Sparus_fst = pd.read_csv(rf'{new_address}/sync_Sp_aurata_LR5371{CHR}_popoolation.fst', sep="\t", header=None)
    sync_Sparus_fet = pd.read_csv(rf'{new_address}/sync_Sp_aurata_LR5371{CHR}_popoolation.fet',  sep="\t", header=None)

    # %%
    # focus on fst differences between the farmed and the wild populations
    # pops1-10 are farmed; pops11-20 are wild

    #populations 5 (80-89) and 8 (119-128) were missing 
    indices = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
                65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
                80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
                94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
                130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
                140, 141, 142, 143, 144, 145, 146, 147, 148, 149
                ]


    sync_Sparus_fst_filtered = sync_Sparus_fst.iloc[:, indices]
    sync_Sparus_fst_filtered

    # %%
    import pandas as pd
    import numpy as np

    # Extract FST values directly from DataFrame without looping through each row manually
    # Split each string at '=' and expand to new DataFrame, then select the second column containing the values
    fst_values = sync_Sparus_fst_filtered.map(lambda x: float(x.split('=')[1]))

    # Calculate mean across columns for each row, this provides the average FST per SNP per comparison
    avg_FST_perSNP_Sparus = fst_values.mean(axis=1)

    # Optional: Print indices at intervals
    for index in range(0, len(avg_FST_perSNP_Sparus), 10000):
        print(index, end=" ")

    # %%
    dfex = pd.DataFrame(avg_FST_perSNP_Sparus)

    # %%
    sync_Sparus_fst_df = pd.DataFrame(avg_FST_perSNP_Sparus, columns=["Fst"])
    sync_Sparus_fst_df['index1'] = sync_Sparus_fst_df.index
    pos = sync_Sparus_fst[1]
    sync_Sparus_fst_df[1] = pos
    sync_Sparus_fst_df[0] = CHR
    sync_Sparus_fst_df

    # %%
    # focus on fet differences between the farmed and the wild populations
    # pops1-10 are farmed; pops11-20 are wild

    #populations 5 (80-89) and 8 (119-128) were missing 
    indices = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
                65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
                80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
                94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
                130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
                140, 141, 142, 143, 144, 145, 146, 147, 148, 149
                ]



    sync_Sparus_fet_filtered = sync_Sparus_fet.iloc[:, indices]
    # %%
    # CORRECT FOR MULTIPLE TESTING WITH THE BH METHOD

    # the "for" loop of above, which estimates the AVERAGE -LOG(P-val) per row
    # modified to correct for multiple test

    AVERAGE_FET_per_SNP_Sparus_BHcor = []
    dropped_indices = []

    for row in range(len(sync_Sparus_fet_filtered)):
        FET_values_per_SNP = []
        if row%10000 == 0:
            print(row, end=" ")
        try:
            for element in range(len(list(sync_Sparus_fet_filtered.iloc[row,:]))):
                num_ = round(float(sync_Sparus_fet_filtered.iloc[row,element].split("=")[1]), 3)
                num_ = 10**(-num_)  # convert num_ from log10(p-value) to p-value
                FET_values_per_SNP.append(num_)
            # run the multiple correction
            FET_values_per_SNP_cor = stats.multipletests(FET_values_per_SNP, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
            # convert the corrected values back to log10
            FET_values_per_SNP_cor = [-math.log10(x) for x in list(FET_values_per_SNP_cor)]
            row_mean_cor = stat.mean(FET_values_per_SNP_cor)
            AVERAGE_FET_per_SNP_Sparus_BHcor.append(row_mean_cor)
        except ValueError:
            dropped_indices.append(row)
            pass

    # %%
    len(dropped_indices)

    # %%
    dropped_indices

    # %%
    sync_Sparus_fet_filtered = sync_Sparus_fet_filtered.drop(index = dropped_indices)
    sync_Sparus_fst_df = sync_Sparus_fst_df.drop(index = dropped_indices)

    # %%
    # make certain to round the resulted averages to 3 decimals
    AVERAGE_FET_per_SNP_Sparus_BHcor = [round(x, 3) for x in AVERAGE_FET_per_SNP_Sparus_BHcor]

    # add averages as new column to the "sync_Sparus_fet_filtered" df
    sync_Sparus_fet_filtered["AVERAGE_SIGN_BHcor"] = AVERAGE_FET_per_SNP_Sparus_BHcor

    # %%
    sync_Sparus_fst_df ['FET_AVG_BHcor'] = sync_Sparus_fet_filtered["AVERAGE_SIGN_BHcor"]

    # %%
    data1 = sync_Sparus_fst_df["Fst"]
    data2 = sync_Sparus_fst_df["FET_AVG_BHcor"]

    fig, ax1 = plt.subplots()

    color = 'tab:blue'
    ax1.set_xlabel(str(CHR) )
    ax1.set_ylabel('Fst', color=color)
    ax1.plot(data1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    plt.ylim(0.05 , 0.2)

    ax2 = ax1.twinx()

    color = 'tab:orange'
    ax2.set_ylabel('AVERAGE_SIGN_BHcor', color=color) 
    ax2.plot(data2, color=color)
    plt.axhline(3, color = 'red',linestyle='dashed' )
    ax2.tick_params(axis='y', labelcolor=color)
    plt.ylim(1.5, 100)

    plt.figure(figsize=(20, 20))
    fig.tight_layout() 
    plt.savefig(rf'Sparus_popoolation_plot_LR5371{CHR}')
    plt.show()

    # %%
    sync_Sparus_fst_df.to_csv(rf"LR5371{CHR}_Sparus_popoolation_result.csv", header=True, index=False, sep=",")
















#TEST
# %%

df = pd.read_csv(rf"LR537121_Sparus_popoolation_result.csv", sep=",")
df
frames = []
CHR = 21
for CHR in range(21,44+1):
    df = pd.read_csv(rf"LR5371{CHR}_Sparus_popoolation_result.csv", sep=",")
    frames.append(df)
# Concatenate all DataFrames into one
merged_df = pd.concat(frames, ignore_index=True)

# Optionally, save the merged DataFrame to a new CSV file
merged_df.to_csv("merged_Sparus_popoolation_results.csv", index=False)


# %% [markdown]
# # -------------------------------- OLD RUN -------------------------------------------------

# %%


# %%
CHR = '21'
sync_Sparus_fst = pd.read_csv(rf'sync_Sp_aurata_LR5371{CHR}_popoolation.fst', sep="\t", header=None)
sync_Sparus_fet = pd.read_csv(rf'sync_Sp_aurata_LR5371{CHR}_popoolation.fet',  sep="\t", header=None)

# sync_Sparus_fst = pd.read_csv(r"C:\Users\arist\Downloads\sync_Sp_aurata_LR5371" + str(CHR) + "_popoolation.fst", sep="\t", header=None)
# sync_Sparus_fet = pd.read_csv(r"C:\Users\arist\Downloads\sync_Sp_aurata_LR5371" + str(CHR) + "_popoolation.fet", sep="\t", header=None, low_memory=False)

# %%
sync_Sparus_fet

# %%
sync_Sparus_fst

# %%
sync_Sparus_fst.iloc[:, 128]



# %%
# focus on fst differences between the farmed and the wild populations
# pops1-10 are farmed; pops11-20 are wild

#populations 5 (80-89) and 8 (119-128) were missing 
indices = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
            49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
            65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
            80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
            94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
            107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
            119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
            130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
            140, 141, 142, 143, 144, 145, 146, 147, 148, 149
            ]


sync_Sparus_fst_filtered = sync_Sparus_fst.iloc[:, indices]
sync_Sparus_fst_filtered

# %%
import pandas as pd
import numpy as np

# Extract FST values directly from DataFrame without looping through each row manually
# Split each string at '=' and expand to new DataFrame, then select the second column containing the values
fst_values = sync_Sparus_fst_filtered.map(lambda x: float(x.split('=')[1]))

# Calculate mean across columns for each row, this provides the average FST per SNP per comparison
avg_FST_perSNP_Sparus = fst_values.mean(axis=1)

# Optional: Print indices at intervals
for index in range(0, len(avg_FST_perSNP_Sparus), 10000):
    print(index, end=" ")

# %%
# Old method - replaced with the above cell
avg_FST_perSNP_Sparus = []

for row in range(len(sync_Sparus_fst_filtered)):
    if row%10000 == 0:
        print(row, end=" ")
    avg_FST_perSNP_perComparison = []
    for element in range(len(list(sync_Sparus_fst_filtered.iloc[row,:]))):
        num_ = float(list(sync_Sparus_fst_filtered.iloc[row,:])[element].split("=")[1])
        avg_FST_perSNP_perComparison.append(num_)
    avg_FST_perSNP_Sparus.append(stat.mean(avg_FST_perSNP_perComparison))

# %%
dfex = pd.DataFrame(avg_FST_perSNP_Sparus)

# %%
sync_Sparus_fst_df = pd.DataFrame(avg_FST_perSNP_Sparus, columns=["Fst"])
sync_Sparus_fst_df['index1'] = sync_Sparus_fst_df.index
pos = sync_Sparus_fst[1]
sync_Sparus_fst_df[1] = pos
sync_Sparus_fst_df[0] = CHR
sync_Sparus_fst_df

# %%
# focus on fet differences between the farmed and the wild populations
# pops1-10 are farmed; pops11-20 are wild

#populations 5 (80-89) and 8 (119-128) were missing 
indices = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
            49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
            65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
            80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
            94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
            107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
            119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
            130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
            140, 141, 142, 143, 144, 145, 146, 147, 148, 149
            ]



sync_Sparus_fet_filtered = sync_Sparus_fet.iloc[:, indices]
sync_Sparus_fet_filtered

# %%
# CORRECT FOR MULTIPLE TESTING WITH THE BH METHOD

# the "for" loop of above, which estimates the AVERAGE -LOG(P-val) per row
# modified to correct for multiple test

AVERAGE_FET_per_SNP_Sparus_BHcor = []
dropped_indices = []

for row in range(len(sync_Sparus_fet_filtered)):
    FET_values_per_SNP = []
    if row%10000 == 0:
        print(row, end=" ")
    try:
        for element in range(len(list(sync_Sparus_fet_filtered.iloc[row,:]))):
            num_ = round(float(sync_Sparus_fet_filtered.iloc[row,element].split("=")[1]), 3)
            num_ = 10**(-num_)  # convert num_ from log10(p-value) to p-value
            FET_values_per_SNP.append(num_)
        # run the multiple correction
        FET_values_per_SNP_cor = stats.multipletests(FET_values_per_SNP, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
        # convert the corrected values back to log10
        FET_values_per_SNP_cor = [-math.log10(x) for x in list(FET_values_per_SNP_cor)]
        row_mean_cor = stat.mean(FET_values_per_SNP_cor)
        AVERAGE_FET_per_SNP_Sparus_BHcor.append(row_mean_cor)
    except ValueError:
        dropped_indices.append(row)
        pass

# %%
len(dropped_indices)

# %%
dropped_indices

# %%
sync_Sparus_fet_filtered = sync_Sparus_fet_filtered.drop(index = dropped_indices)
sync_Sparus_fst_df = sync_Sparus_fst_df.drop(index = dropped_indices)

# %%
# make certain to round the resulted averages to 3 decimals
AVERAGE_FET_per_SNP_Sparus_BHcor = [round(x, 3) for x in AVERAGE_FET_per_SNP_Sparus_BHcor]

# add averages as new column to the "sync_Sparus_fet_filtered" df
sync_Sparus_fet_filtered["AVERAGE_SIGN_BHcor"] = AVERAGE_FET_per_SNP_Sparus_BHcor

# %%
sync_Sparus_fet_filtered

# %%
sync_Sparus_fst_df ['FET_AVG_BHcor'] = sync_Sparus_fet_filtered["AVERAGE_SIGN_BHcor"]
sync_Sparus_fst_df

# %%
data1 = sync_Sparus_fst_df["Fst"]
data2 = sync_Sparus_fst_df["FET_AVG_BHcor"]

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel(str(CHR) )
ax1.set_ylabel('Fst', color=color)
ax1.plot(data1, color=color)
ax1.tick_params(axis='y', labelcolor=color)
plt.ylim(0.05 , 0.2)

ax2 = ax1.twinx()

color = 'tab:orange'
ax2.set_ylabel('AVERAGE_SIGN_BHcor', color=color) 
ax2.plot(data2, color=color)
plt.axhline(3, color = 'red',linestyle='dashed' )
ax2.tick_params(axis='y', labelcolor=color)
plt.ylim(1.5, 100)

plt.figure(figsize=(20, 20))
fig.tight_layout() 
plt.savefig(rf'Sparus_popoolation_plot_LR5371{CHR}')
plt.show()

# %%
sync_Sparus_fst_df.to_csv(rf"LR5371{CHR}_Sparus_popoolation_result.csv", header=True, index=False, sep=",")


