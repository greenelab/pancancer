"""
Gregory Way 2017
PanCancer Classifier
scripts/initialize/process_sample_freeze.py

Takes in sample freeze data that was determined by TCGA PanCancer Atlas
consortium along with raw RNAseq and mutation data. The script will process
the datasets and subset each according to the frozen samples.

Usage: Run once in command line

        python scripts/initialize/process_sample_freeze.py

Output:
RNAseq and mutation data subset by sample freeze
"""

import os
import numpy as np
import pandas as pd

# Input Files
rna_file = os.path.join('data', 'raw', 'pancan_normalized_rnaseq.tsv')
mut_file = os.path.join('data', 'raw', 'mc3.v0.2.8.PUBLIC.maf')
sample_freeze_file = os.path.join('data', 'raw',
                                  'sampleset_freeze_version4_modify.csv')

# Output Files
rna_out_file = os.path.join('data', 'pancan_rnaseq_freeze.tsv')
mut_out_file = os.path.join('data', 'pancan_mutation_freeze.tsv')
freeze_out_file = os.path.join('data', 'sample_freeze.tsv')
burden_out_file = os.path.join('data', 'mutation_burden_freeze.tsv')

# Load Data
rnaseq_df = pd.read_table(rna_file, index_col=0)
mutation_df = pd.read_table(mut_file)
sample_freeze_df = pd.read_csv(sample_freeze_file)

# Process RNAseq file
rnaseq_df.index = rnaseq_df.index.map(lambda x: x.split('|')[0])
rnaseq_df.columns = rnaseq_df.columns.str.slice(start=0, stop=15)
rnaseq_df = rnaseq_df.drop('?').fillna(0).sort_index(axis=1)

# Gene is listed twice in RNAseq data, drop both occurrences
rnaseq_df.drop('SLC35E2', axis=0, inplace=True)
rnaseq_df = rnaseq_df.T

# Determine consistent sample freeze in RNAseq
freeze_barcodes = set(sample_freeze_df.SAMPLE_BARCODE)
freeze_barcodes = freeze_barcodes.intersection(set(rnaseq_df.index))

# Process Mutation File
mutation_df = mutation_df.assign(PATIENT_BARCODE=mutation_df
                                 .Tumor_Sample_Barcode
                                 .str.slice(start=0, stop=12))
mutation_df = mutation_df.assign(SAMPLE_BARCODE=mutation_df
                                 .Tumor_Sample_Barcode
                                 .str.slice(start=0, stop=15))

# Determine consistent sample freeze between RNAseq and mutation
mut_samples = set(mutation_df.SAMPLE_BARCODE.unique())
freeze_barcodes = freeze_barcodes.intersection(mut_samples)
freeze_barcodes = sorted(freeze_barcodes)

# Subset rnaseq data to only barcodes and remove duplicate rows
rnaseq_df = rnaseq_df.loc[freeze_barcodes, :]
rnaseq_df = rnaseq_df[~rnaseq_df.index.duplicated()]
rnaseq_df.to_csv(rna_out_file, sep='\t', index_col=0)

# Filter mutation types and generate binary matrix
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
mut_pivot = (mutation_df.query("Variant_Classification in @mutations")
                        .groupby(['SAMPLE_BARCODE', 'Chromosome',
                                  'Hugo_Symbol'])
                        .apply(len).reset_index()
                        .rename(columns={0: 'mutation'}))

mut_pivot = (mut_pivot.pivot_table(index='SAMPLE_BARCODE',
                                   columns='Hugo_Symbol', values='mutation',
                                   fill_value=0)
                      .astype(bool).astype(int))

# 12 Samples don't have any deleterious mutations
# This command will introduce NAs for these 12 samples, fill them with zeros
mut_pivot = mut_pivot.loc[freeze_barcodes, :]
mut_pivot = mut_pivot.fillna(0)
mut_pivot = mut_pivot.astype(int)
mut_pivot.to_csv(mut_out_file, sep='\t', index_col=0)

# Generate a mutation burden variable (log10 total deleterious mutations)
burden_df = mutation_df[mutation_df['Variant_Classification'].isin(mutations)]
burden_df = burden_df.groupby('SAMPLE_BARCODE').apply(len)
burden_df = np.log10(burden_df)
burden_df = burden_df.loc[freeze_barcodes]
burden_df = burden_df.fillna(0)
burden_df = pd.DataFrame(burden_df, columns=['log10_mut'])
burden_df.to_csv(burden_out_file, sep='\t')

# Write out finalized and subset sample freeze file
sample_freeze_df = sample_freeze_df[sample_freeze_df.SAMPLE_BARCODE
                                                    .isin(freeze_barcodes)]
sample_freeze_df.to_csv(freeze_out_file, sep='\t')
