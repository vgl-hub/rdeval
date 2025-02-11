import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')
plt.rcParams['svg.fonttype'] = 'none'

total_size=pd.read_csv("original_dataset_size.tsv",sep='\t', names=['Accession','Size'])
df1_illumina=pd.read_csv("final_homopolymer_compression_illumina.tsv",sep=' ', names=['Accession','HP-Size'])
df2_pacbio=pd.read_csv("final_homopolymer_compression_pacbio.tsv",sep=' ', names=['Accession','HP-Size'])
df1_illumina["Platform"] = "Illumina"
df2_pacbio["Platform"] = "PacBio"

concatenated_df = pd.concat([df1_illumina, df2_pacbio], ignore_index=True)
merged_df=concatenated_df.merge(total_size, on="Accession")

custom_colors = {'Illumina': "#154CAC", 'PacBio': "#EF476F"}
# Map the 'Platform' column to colors
merged_df['Color'] = merged_df['Platform'].map(custom_colors)

# Create the scatter plot
fig, ax = plt.subplots(figsize=(10, 10))
scatter = ax.scatter(
    merged_df['Size'] / 1000000000, 
    merged_df['HP-Size'] / 1000000000, 
    c=merged_df['Color'],  # Use the mapped colors
    alpha=0.7, s=100
)

# Add labels and title
ax.set_xlabel("Original size (Gbp)", fontsize=26)
ax.set_ylabel("HC size (Gbp)", fontsize=26)
ticks = [0, 25, 50, 75, 100, 125, 150]
#ax.tick_params(ticks,axis='x', labelsize=26)  # Sets x-tick label font size to 26
#ax.tick_params(ticks,axis='y', labelsize=26)  # Sets y-tick label font size to 26
plt.xticks(ticks, fontsize=26)
plt.yticks(ticks, fontsize=26)

for platform, color in custom_colors.items():
    ax.scatter([], [], color=color, label=platform, s=100)
ax.legend(title="Platform", fontsize=20, title_fontsize=22, loc="upper left")
plt.xlim(0,150)
plt.ylim(0,150)
# Add a dashed diagonal line from the bottom-left to the top-right
min_val = min(plt.xlim()[0], plt.ylim()[0])  # Get the minimum value of axes
max_val = max(plt.xlim()[1], plt.ylim()[1])  # Get the maximum value of axes
plt.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="black", linewidth=1, label="Diagonal")

# Ensure the axes are equal to make the diagonal accurate
plt.axis("equal")
# Adjust layout and save the plot
fig.tight_layout()
plt.savefig("Figure_Suppl_figure.svg",format="svg", dpi=300)
