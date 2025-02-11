import pandas as pd
import os
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import matplotlib.patches as mpatches
import matplotlib
matplotlib.use('Agg')


sequencing_platform=pd.read_csv("SRS_platform.tsv",sep='\t',names=['Accession','Type'])
stats=pd.read_csv("filtered_rdeval.tsv",sep='\t',names=['Accession','Average length','Average quality'])

merged_df=stats.merge(sequencing_platform, on = "Accession")

percentiles = [5,50,95]
merged_df["Average length"] = pd.to_numeric(merged_df["Average length"], errors="coerce")
merged_df["Average quality"] = pd.to_numeric(merged_df["Average quality"], errors="coerce")

length_percentiles = {p: merged_df["Average length"].quantile(p / 100) for p in percentiles}
custom_colors = {
    "HiFi_Revio": "#40a408",
    "CLR_Sequel": "#e54ecc",
    "CLR_Sequel II": "#afd18d",
    "HiFi_Sequel II": "#8B0000"
}

plt.figure(figsize=(10, 10))
plt.scatter(
    merged_df["Average length"],
    merged_df["Average quality"],
    alpha=0.9,
    c=merged_df["Type"].map(custom_colors),  # Use the color mapping
    edgecolors="black",
    s=100
)
for p in percentiles:
    plt.axvline(length_percentiles[p], color="black", linestyle="--", label=f"{p}th Length: {length_percentiles[p]:.2f}")
plt.text(length_percentiles[5]  + 100 , 44 , "5th" , fontsize=16, color="black")
plt.text(length_percentiles[50] + 100 , 44 , "50th", fontsize=16, color="black")
plt.text(length_percentiles[95] + 100 , 44 , "95th", fontsize=16, color="black")

# Add labels and title
plt.xlabel("Average read length (bp)", fontsize=26)
plt.ylabel("Average read quality (mean BAM QUAL)", fontsize=26)
plt.ylim(-1,45)
plt.xticks(fontsize=26)
yticks = [0,10,20,30,40]
ytick_labels = ["NA" if tick == 0 else str(int(tick)) for tick in yticks]
plt.yticks(yticks, ytick_labels,fontsize=26)

# Create legend patches for each category
legend_patches = [
    mpatches.Patch(color=color, label=type)
    for type, color in custom_colors.items()
]
# Add the legend to the plot
plt.legend(handles=legend_patches, title="Sequencing Type")
plt.grid(False)
plt.tight_layout()
plt.savefig("Figure_2b.svg",format='svg', dpi=300)
