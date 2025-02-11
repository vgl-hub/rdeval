import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import matplotlib
plt.rcParams['svg.fonttype'] = 'none'
matplotlib.use('Agg')

output_dir = "info_instruments"
os.makedirs(output_dir, exist_ok=True)

# Step 1: Import data
vgp_assembly_dataset_nxSeqRun = pd.read_csv("new_rdeval_30.01.filtered.tsv", header=None, sep="\0")

# Step 2: Split and organize columns
vgp_assembly_dataset_nxSeqRun = vgp_assembly_dataset_nxSeqRun[0].str.split('\t',n=1, expand=True)
vgp_assembly_dataset_nxSeqRun.columns = ['Accession', 'Data']

values = vgp_assembly_dataset_nxSeqRun['Data'].str.split(';', expand=True)

vgp_assembly_dataset_nxSeqRun.drop('Data', axis=1, inplace=True)

vgp_assembly_dataset_nxSeqRun = pd.concat([vgp_assembly_dataset_nxSeqRun, values], axis=1)

vgp_assembly_dataset_nxSeqRun = vgp_assembly_dataset_nxSeqRun.melt(
        id_vars=['Accession'], var_name="Data", value_vars=range(values.shape[1])
    ).dropna()

vgp_assembly_dataset_nxSeqRun[['Length', 'Nr','Partial','Size']] = vgp_assembly_dataset_nxSeqRun['value'].str.split("\t", expand=True).apply(pd.to_numeric)


vgp_assembly_dataset_nxSeqRun.drop('value', axis=1, inplace=True)

sequencing_platform = pd.read_csv("SRS_platform.tsv", sep ='\t', names=['Accession','Instrument'])

vgp_assembly_dataset_nxSeqRun= vgp_assembly_dataset_nxSeqRun.merge(sequencing_platform, on ="Accession")

custom_colors = {
    "HiFi_Revio": "#40a408",
    "CLR_Sequel": "#e54ecc",
    "CLR_Sequel II": "#afd18d",
    "HiFi_Sequel II": "#8B0000"
}

def process_group(grp):
    # Sort the group by 'Length' in ascending order
    grp = grp.sort_values(by="Length", ascending=False)
    total_size = grp["Partial"].sum()
    grp["Nx"] =  grp["Size"]/total_size
    return grp

accession_agg = vgp_assembly_dataset_nxSeqRun.groupby("Accession").apply(process_group).reset_index(drop=True)
platform_results = []
num_points = 5000
x_values = np.linspace(0, 1, num_points)  # Generate 5000 equally spaced x-values

for platform, platform_grp in accession_agg.groupby('Instrument'):
    aggregated_data = []

    for x in x_values:
        # Select data points near the current x-value
        nearby = platform_grp[np.abs(platform_grp['Nx'] - x) <= (1 / num_points)]
        if not nearby.empty:
            # Calculate statistics for the current range
            mean = nearby['Length'].mean()
            max_val = nearby['Length'].max()
            min_val = nearby['Length'].min()
            aggregated_data.append((x, mean, max_val, min_val))

    # Store results for each platform
    platform_results.append((platform, pd.DataFrame(aggregated_data, columns=["Nx", "Mean", "Max", "Min"])))

# Step 3: Plot Results
fig, ax = plt.subplots()
ax.set_yscale('log')

for platform, agg_df in platform_results:
    color = custom_colors.get(platform, "#000000")  # Default to black if color not specified

    # Plot median as a solid line
    ax.plot(agg_df["Nx"]*100, agg_df["Mean"], label=f'{platform} Mean', color=color)
    # Save agg_df to a CSV inside platform-specific folder
    platform_dir = os.path.join(output_dir, platform)
    os.makedirs(platform_dir, exist_ok=True)  # Create platform-specific directory if it doesn't exist

    file_path = os.path.join(platform_dir, f"{platform}_agg_df.csv")
    agg_df.to_csv(file_path, index=False)  # Save as CSV
    print(f"Saved: {file_path}")

# Finalize plot
def bp_kbp_formatter(x, pos):
    if x >= 1000:  # Convert to Kbp
        return f"{int(x/1000)} Kbp"
    else:  # Show values below 100 as bp
        return f"{int(x)} bp"

ax.set_xlabel('Nx', fontsize=26)
ax.set_ylabel('Length', fontsize=26)
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_yticks([ 100, 1000, 10000, 20000,100000])

ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(bp_kbp_formatter))
ax.legend(loc='lower left')


plt.tight_layout()
plt.savefig('Figure_2a.svg',dpi=300, format='svg')

