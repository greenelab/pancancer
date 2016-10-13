"""
Gregory Way 2016
Modified from https://github.com/cognoma/cancer-data/
PanCancer NF1 Classifier
subset_datasets.py

Usage: Run in command line with required command argument:

        python subset_datasets.py --database {"Xena", "Synapse"}

Output:
Cleaned RNAseq, mutation, and clinical data. Each dataset is subset to the
intersection of samples measured on all three platforms. The script also
outputs a heatmap for proportion of NF1 and RAS mutations across tissues.
"""

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure.max_open_warning = 0
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("paper", rc={"font.size": 8, "axes.titlesize": 20,
                             "axes.labelsize": 20, "xtick.labelsize": 12,
                             "ytick.labelsize": 8})

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--database', help='Which data to subset and clean')
args = parser.parse_args()

data = args.database

if data == 'Synapse':
    rnaseq_df = pd.read_table('data/raw/pancan_normalized_rnaseq.tsv',
                              index_col=0)
    mutation_df = pd.read_table('data/raw/pancan_mutation.maf')
    clinical_df = pd.read_table('data/raw/pancan_clinical.tsv', index_col=1,
                                encoding='ISO-8859-1', low_memory=False)
elif data == 'Xena':
    mutation_df = pd.read_table('data/xena/gbm_download/PANCAN_mutation')

else:
    raise ValueError('Enter either "Synapse" or "Xena" using "--data" flag')

if data == 'Synapse':
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

    # Process synapse mutations
    mut_df = (mut_df.query("Variant_Classification in @mutations")
                    .groupby(['TCGA_ID', 'Chromosome', 'Hugo_Symbol'])
                    .apply(len).reset_index().rename(columns={0: 'mutation'}))

    mut_df = (mut_df.pivot_table(index='TCGA_ID', columns='Hugo_Symbol',
                                 values='mutation', fill_value=0)
                    .astype(bool).astype(int))

    # Write subsets to file
    mut_df.to_csv('data/mutation_table.tsv', sep='\t')
    rna_df.to_csv('data/rnaseq_data.tsv', sep='\t')
    clin_df.to_csv('data/clinical_data.tsv', sep='\t')

elif data == 'Xena':
    # Process Xena mutations
    mutation_df = (mutation_df.query("effect != 'Silent'")
                   .groupby(['#sample', 'chr', 'gene'])
                   .apply(len).reset_index()).rename(columns={0: 'mutation'})

    mut_df = (mutation_df.pivot_table(index='#sample', columns='gene',
              values='mutation', fill_value=0) .astype(bool).astype(int))
    mut_df.to_csv('data/xena/xena_mutation_table.tsv', sep='\t')
