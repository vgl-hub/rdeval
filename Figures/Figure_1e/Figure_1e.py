import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import linregress
from matplotlib.lines import Line2D

matplotlib.use('Agg')  # Use a non-GUI backend

plt.rcParams['svg.fonttype'] = 'none'

#import empirical .rd size along downsampling process
emp_df = pd.read_csv("emp_rd_size.tsv", sep="\t", names=["Sample","Downsampling","Size"])

#import theoretical projection of .rd size along downsampling process
the_df = pd.read_csv("projection.tsv", sep="\t", names=["Sample","Downsampling","Size"])

#manipulation on dataframes
emp_df["Downsampling"] = emp_df["Downsampling"].astype(float)
grouped = emp_df.groupby("Sample")

the_df["Downsampling"] = the_df["Downsampling"].astype(float)
the_df['Sample'] = the_df["Sample"].str.replace("_p", "", regex=False)
merged_grouped = the_df.groupby("Sample")

# Example: Replace 'grouped' and 'merged_grouped' with your actual grouped DataFrames
samples = set([sample for sample, _ in grouped] + [sample for sample, _ in merged_grouped])

# Define your custom colors for each sample
custom_colors = {'CHM13.illumina': "#154CAC", 'CHM13.Hifi': "#EF476F", 'CHM13.ONT': "#118AB2", "ASHG-C063_R1.AVITI": "#06D6A0" , "CHM13.10x":"#FFD166"}

# Initialize a figure with a specific size
plt.figure(figsize=(10, 10))

# Scatter plot for the first dataset with consistent coloring
for sample, data in grouped:
    plt.scatter(
        data["Downsampling"],            # X-axis: Downsampling ratio
        data["Size"] / 1048576,          # Y-axis: File size converted to Mb
        label=sample,                    # Label for the sample
        alpha=0.9,                       # Transparency for points
        s=80,                            # Marker size
        color=custom_colors[sample]      # Color based on the custom color dictionary
    )

# Scatter plot and trendline for the second dataset (merged_grouped)
for sample, data in merged_grouped:
    x = data["Downsampling"].sort_values()          # X-axis: Downsampling ratio
    y = data["Size"].sort_values() / 1048576        # Y-axis: File size converted to Mb

    # Add a trendline using linear regression
    slope, intercept, _, _, _ = linregress(x, y)    # Calculate slope and intercept of the regression line
    trendline = slope * x + intercept               # Compute trendline values
    plt.plot(
        x, 
        trendline, 
        color=custom_colors[sample],               # Match the color with scatter points
        linestyle="--",                            # Dashed line for projection
        label=f"{sample} Projection"              # Label for the trendline
    )

# Set x-axis and y-axis labels
plt.xlabel("Downsampling", fontsize=26)
plt.ylabel(".rd file size (Mb)", fontsize=26)

# Customize tick label font sizes
plt.xticks(fontsize=26)
plt.yticks(fontsize=26)

plt.grid(False)

plt.tight_layout()
plt.legend(loc="upper left")

plt.savefig("Figure_1e.svg", format="svg", dpi=300)
