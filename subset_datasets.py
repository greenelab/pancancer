"""
Gregory Way 2016
Modified from https://github.com/cognoma/cancer-data/
PanCancer NF1 Classifier
subset_datasets.py

Usage: Run in command line `python subset_datasets.py`

Output:
Cleaned RNAseq, mutation, and clinical data. Each dataset is subset to the
intersection of samples measured on all three platforms. The script also
outputs a heatmap for proportion of NF1 and RAS mutations across tissues.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure.max_open_warning = 0
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("paper", rc={"font.size": 8, "axes.titlesize": 20,
                             "axes.labelsize": 20, "xtick.labelsize": 12,
                             "ytick.labelsize": 8})

rnaseq_df = pd.read_table('data/raw/pancan_normalized_rnaseq.tsv', index_col=0)
mutation_df = pd.read_table('data/raw/pancan_mutation.maf')
clinical_df = pd.read_table('data/raw/pancan_clinical.tsv', index_col=1,
                            encoding='ISO-8859-1', low_memory=False)

# Assign TCGA_ID from substring barcode
mutation_df = mutation_df.assign(TCGA_ID=(mutation_df.Tumor_Sample_Barcode
                                          .str.slice(start=0, stop=12)))

# Determine which samples to use for RNAseq data
rna_ids = pd.DataFrame(rnaseq_df.columns, columns=['full'])
rna_ids = rna_ids.assign(base=rna_ids.full.str.slice(start=0, stop=12))
rna_ids = rna_ids.assign(tumor=rna_ids.full.str.slice(start=13, stop=15))

# Primary solid tumor and primary blood tumor
rna_ids = rna_ids[rna_ids.tumor.isin(['01', '03'])]

# Some primary samples are measured twice
rna_ids = rna_ids[~rna_ids.base.duplicated()]

# Subset data to common samples
common_samples = list(set(rna_ids.base.unique()) &
                      set(clinical_df.index) &
                      set(mutation_df.TCGA_ID.unique()))
rna_ids = rna_ids[rna_ids.base.isin(common_samples)]

clin_df = clinical_df[clinical_df.index.isin(common_samples)].sort_index()
mut_df = mutation_df[mutation_df.TCGA_ID.isin(common_samples)]
rna_df = rnaseq_df.loc[:, rna_ids.full]
rna_df.columns = rna_ids.base

# Subset gene identifiers
rna_df.index = rna_df.index.map(lambda x: x.split('|')[0])
rna_df = rna_df.drop('?').fillna(0).sort_index(axis=1)

# Total study sample size
study_sample_size = clin_df.acronym.value_counts()

# Filter mutation types
mutations = {
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'In_Frame_Del',
    'In_Frame_Ins',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Nonstop_Mutation',
    'RNA',
    'Splice_Site',
    'Translation_Start_Site',
}

mut_df = (mut_df.query("Variant_Classification in @mutations")
                .groupby(['TCGA_ID', 'Chromosome', 'Hugo_Symbol'])
                .apply(len).reset_index().rename(columns={0: 'mutation'}))

mut_df = (mut_df.pivot_table(index='TCGA_ID', columns='Hugo_Symbol',
                             values='mutation', fill_value=0)
                .astype(bool).astype(int))

# Write subsets to file
mut_df.to_csv('data/mutation_table.tsv', sep='\t', header=True)
rna_df.to_csv('data/rnaseq_data.tsv', sep='\t', header=True)
clin_df.to_csv('data/clinical_data.tsv', sep='\t', header=True)

# Explore mutation counts
mut_df = mut_df.join(clin_df['acronym'], how='inner')
mut_df = mut_df.assign(ALL_RAS=mut_df[['KRAS', 'NRAS', 'HRAS']].max(axis=1))
mutation_heatmap = mut_df[['NF1', 'KRAS', 'NRAS', 'HRAS', 'ALL_RAS',
                           'acronym']].groupby('acronym').sum()

# Show mutation frequency heatmap
heatmap_df = mutation_heatmap.divide(mut_df.acronym.value_counts(sort=False)
                                     .sort_index(), axis=0)
ax = sns.heatmap(heatmap_df, annot=True)
ax.set(xlabel='gene', ylabel='study')
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig('figures/nf1_ras_mutation_summary.png')
plt.close()

mutation_table = mutation_heatmap.join(study_sample_size)
mutation_table.to_csv('NF1_RAS_mutation_summary.tsv', sep='\t')
