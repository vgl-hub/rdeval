import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import pandas as pd
matplotlib.use('Agg')  # Use a non-GUI backend

# Preprocess the data: Convert Size to GB and Runtime to minutes
merged_df=pd.read_csv("runtime_size.tsv", sep='\t', names=['Sample Name','Runtime','File Type','Sequencing Technology','Size'])
merged_df['Size'] = merged_df['Size'] / 1000000000  # Convert bytes to GB
merged_df['Runtime'] = merged_df['Runtime'] / 60       # Convert seconds to minutes

custom_colors = {'Illumina': "#154CAC", 'PacBio': "#EF476F"}

plt.figure(figsize=(10, 10))

# Initialize a list to track annotation positions
annotation_positions = []

# Loop through each group of 'Sequencing Technology' and 'File Type'
unique_groups = merged_df.groupby(['Sequencing Technology', 'File Type'])

for (tech, file_type), group_data in unique_groups:
    color = custom_colors[tech] if tech in custom_colors else ('#1f77b4' if color_switch else '#2ca02c')

    sns.regplot(
        data=group_data,
        x='Size',
        y='Runtime',
        scatter=False,  # Enable scatter plot
        label=f"{tech} - {file_type}",  # Add a unique label for each trendline
        line_kws={'linewidth': 2, 'linestyle': '-', 'color': color},  # Customize line style and color
        scatter_kws={'s': 40, 'color': color}#, 'marker': marker_styles}  # Customize scatter marker style
    )

    # Annotate maximum points manually to better visualization of labels
    max_point = group_data.loc[group_data['Runtime'].idxmax()]
    if file_type == 'fastq.gz' and tech =='PacBio':
        plt.text(max_point['Size'] +2, max_point['Runtime'] , "fastq.gz", fontsize=14, color="black")
    elif file_type == 'fastq.gz' and tech =='Illumina':
        plt.text(max_point['Size']  +2, max_point['Runtime'] , "fastq.gz", fontsize=14, color="black")
    elif file_type == 'bam' and tech == 'PacBio':
        plt.text(max_point['Size'] +2, max_point['Runtime'] -2, "bam", fontsize=14, color="black")
    elif file_type == 'bam' and tech == 'Illumina':
         plt.text(max_point['Size'] +2, max_point['Runtime'] , "bam", fontsize=14, color="black")
    elif file_type == 'fasta' and tech =="Illumina":
        plt.text(max_point['Size'] +2, max_point['Runtime'] +1, "fasta", fontsize=14, color="black")
    elif file_type == 'fasta' and tech =="PacBio":
        plt.text(max_point['Size'] +2, max_point['Runtime'] -1, "fasta", fontsize=14, color="black")
    elif file_type == 'cram' and tech =="PacBio":
        plt.text(max_point['Size'] +2, max_point['Runtime']  +1.5, "cram", fontsize=14, color="black")
    elif file_type == 'cram' and tech =="Illumina":
        plt.text(max_point['Size'] +2, max_point['Runtime']  , "cram", fontsize=14, color="black")
    elif file_type == 'fastq' and tech == 'Illumina':
        plt.text(max_point['Size'] +2, max_point['Runtime'] +2, "fastq", fontsize=14, color="black")
    elif file_type == 'fastq' and tech == 'PacBio':
        plt.text(max_point['Size'] +2, max_point['Runtime']  , "fastq", fontsize=14, color="black")
    elif file_type == 'fasta.gz' and tech =='Illumina':
        plt.text(max_point['Size'] +2, max_point['Runtime'] , "fasta.gz", fontsize=14, color="black")
    elif file_type == 'fasta.gz' and tech =='PacBio':
        plt.text(max_point['Size'] +2, max_point['Runtime'] , "fasta.gz", fontsize=14, color="black")

# Add labels, title, and legend
plt.xlabel('Dataset length (Gb)', fontsize=26)
plt.ylabel('Runtime (minutes)', fontsize=26)
plt.xlim(-5, 160)
plt.ylim(-5,65)
plt.xticks(fontsize=26)
plt.yticks(fontsize=26)

plt.tight_layout()
plt.grid(False)
plt.savefig("Figure_1f.svg", format="svg", dpi=300)
