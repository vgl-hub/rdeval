import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
matplotlib.use('Agg')  # Use a non-GUI backend
plt.rcParams['svg.fonttype'] = 'none'

combined_df=pd.read_csv("combined_filetype_size.tsv",sep='\t',names=['Sample Name','Size', 'File Type','Sequencing Technology'])

# Filter for Illumina data
illumina_data = combined_df[combined_df["Sequencing Technology"].str.lower() == "illumina"]

# Filter for PacBio data
pacbio_data = combined_df[combined_df["Sequencing Technology"].str.lower() == "pacbio"]

illumina_order = (
    combined_df[combined_df["Sequencing Technology"] == "illumina"]
    .groupby("File Type")["Size"]
    .median()
    .sort_values(ascending=False)
    .index
)

pacbio_order = (
    combined_df[combined_df["Sequencing Technology"] == "pacbio"]
    .groupby("File Type")["Size"]
    .median()
    .sort_values(ascending=False)
    .index
)

fig, axes = plt.subplots(1, 2, figsize=(10, 10), sharey=True)

# Plot Illumina data
sns.boxplot(
    data=illumina_data,
    x="File Type",
    y="Size",
    order=illumina_order,
    palette=["#154CAC"] * len(illumina_order),  
    ax=axes[0],
    showcaps=True,
    whiskerprops={'color': 'black', 'linewidth': 1.5},
    boxprops={'edgecolor': 'black', 'linewidth': 1.5},
    medianprops={'color': 'black', 'linewidth': 1.5},
    showfliers=False
)
axes[0].set_title("Illumina", fontsize=26)
axes[0].set_ylabel("File size log$_{10}$(bytes)", fontsize=26)
axes[0].set_yscale("log")
axes[0].tick_params(axis='x', rotation=45,labelsize=26)
axes[0].tick_params(axis='y', labelsize=26)


# Plot PacBio data
sns.boxplot(
    data=pacbio_data,
    x="File Type",
    y="Size",
    order=pacbio_order,
    palette=["#EF476F"] * len(pacbio_order),
    ax=axes[1],
    showcaps=True,
    whiskerprops={'color': 'black', 'linewidth': 1.5},
    boxprops={'edgecolor': 'black', 'linewidth': 1.5},
    medianprops={'color': 'black', 'linewidth': 1.5},
    showfliers=False
)
axes[1].set_title("PacBio", fontsize=26)
axes[1].set_ylabel("File size log$_{10}$(bytes)", fontsize=26)
axes[1].tick_params(axis='x', rotation=45,labelsize=26)

# Adjust layout for better spacing
plt.tight_layout(rect=[0, 0, 1, 0.95])

plt.savefig("Figure_1b.svg", format="svg",dpi=300)
