import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['svg.fonttype'] = 'none' # preserve text labels as text when saving svg file

# Read and combine data frames with assembly statistics and platform info
assembly_stats = pd.read_csv('assembly_stats_cleaned.tsv', sep = '\t')
assembly_stats = assembly_stats.loc[assembly_stats['Haplotype'] != 'Alternate']
platforms = pd.read_csv('new_sequencing_platform.tsv', sep = '\t')

merged = assembly_stats.merge(platforms)
merged.Platform =  merged.Platform.str.split('_').str[-1]

# Mark assemblies by species class
merged['Class'] = merged['Assembly Name'].str[0]
merged['Class'] = merged['Class'].map({
    'a': 'Amphibian',
    'b': 'Bird',
    'f': 'Fish',
    'm': 'Mammal',
    'r': 'Reptile',
    's': 'Shark'
})

# Get Total Sequence Read Length from SRA metadata
sra_metadata = pd.read_csv('SRA_metadata.tsv', sep = '\t')
sra_metadata = sra_metadata.dropna(axis = 0, subset = 'Type (CLR/HiFi)')
sra_metadata = sra_metadata.loc[((((sra_metadata['Type (CLR/HiFi)'] == 'HiFi') & 
                                           (sra_metadata['LIBRARY_NAME'].str.contains('fastq', na = False))))) | 
                                (sra_metadata['Type (CLR/HiFi)'] == 'CLR')]

total_sequence_length = sra_metadata.groupby(by = ['Sample@acc'], as_index = False)['Statistics@total_bases'].sum()
total_sequence_length.columns = ['SRS', 'Total Sequence Length']

# Calculate coverage from Total Sequence Length and Total Assembly Length
merged = merged.merge(total_sequence_length, on = 'SRS')
merged['Coverage'] = merged['Total Sequence Length']/merged['Total Assembly Length']

# Filter HiFi assemblies
merged_hifi = merged[(merged['Platform'] == 'HiFi')]
merged_hifi['Contig N50 (Mb)'] = merged_hifi['Contig N50']/1000000

# Mark outliers
outliers = ['mThoBot2.hap1', 'mThoBot2.hap2', 'rCycPin1.hap1', 'rCycPin1.hap2', 'bSarPap1.hap1', 'fLatCha1.pri']
merged_hifi['is_outlier'] = merged_hifi['Assembly Name'].isin(outliers)

# Save final dataframe to tsv
merged_hifi.to_csv('Figure_2c_data.tsv', sep = '\t', index = False)

# Set up plotting
plt.style.use('seaborn-v0_8-ticks')
fig, axs = plt.subplots(2, 3, sharex = False, sharey = False, layout = 'constrained')

# Need to generate facets iteratively in order to superimpose regression line on each scatterplot
for ax, species_class in zip(axs.flat, pd.unique(merged_hifi['Class'])):
    ax.set_title(f'{species_class}')
    sns.scatterplot(
        merged_hifi[(merged_hifi['Platform'] == 'HiFi') & (merged_hifi['Class'] == species_class)], hue = 'is_outlier',
            palette = ['#000000', '#C1272D'],
            x = 'Coverage', y = 'Contig N50 (Mb)', legend = False,
            ax = ax
    )
    sns.regplot(
        merged_hifi[(merged_hifi['Platform'] == 'HiFi') & (merged_hifi['Class'] == species_class) & ~merged_hifi['is_outlier']],
            x = 'Coverage', y = 'Contig N50 (Mb)', scatter = False, truncate = False, line_kws = dict(color='#000000'),
            ax = ax
    )
plt.savefig('Figure_2c.svg')