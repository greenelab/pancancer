"""
Gregory Way 2017
PanCancer Classifier
build_heatmap_data.py

Usage: Run in command line with required command argument:

        python build_heatmap_data.py --genes $GENES

Where GENES is a comma separated string. There are also optional arguments:

    --copy_number       store_true: supplement copy number to define Y
    --remove_hyper      store_true: remove hypermutated samples
    --alt_folder        string: location of where to save output files

Output:
.tsv file of gene events per TCGA disease types saved in the tables/ folder
"""

import os
import warnings
import pandas as pd
import argparse
from tcga_util import integrage_copy_number

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genes',
                    help='Comma separated string of HUGO gene symbols')
parser.add_argument('-u', '--copy_number', action='store_true',
                    help='Supplement Y matrix with copy number events')
parser.add_argument('-v', '--remove_hyper', action='store_true',
                    help='Remove hypermutated samples')
parser.add_argument('-f', '--alt_folder', default='Auto',
                    help='location to save')
args = parser.parse_args()

# Load command arguments
genes = args.genes.split(',')
copy_number = args.copy_number
remove_hyper = args.remove_hyper
alt_folder = args.alt_folder

if alt_folder == 'Auto':
    base_folder = 'tables'
else:
    base_folder = os.path.join(alt_folder, 'tables')

if not os.path.exists(base_folder):
    os.makedirs(base_folder)

gene_output_file = os.path.join(base_folder, '{}_mutations_heatmap_data.csv'
                                .format(args.genes.replace(',', '_')))

# Load Datasets
expr_file = os.path.join('data', 'pancan_rnaseq_freeze.tsv')
mut_file = os.path.join('data', 'pancan_mutation_freeze.tsv')
clin_file = os.path.join('data', 'clinical_data.tsv')
sample_freeze_file = os.path.join('data', 'sampleset_freeze_version3.csv')
copy_loss_file = os.path.join('data', 'copy_number_loss_status.tsv')
copy_gain_file = os.path.join('data', 'copy_number_gain_status.tsv')
cancer_gene_file = os.path.join('data', 'vogelstein_cancergenes.tsv')
mut_burden_file = os.path.join('ddr', 'data', 'mutation-load.txt')

rnaseq_df = pd.read_table(expr_file, index_col=0)
mutation_df = pd.read_table(mut_file, index_col=0)
clinical_df = pd.read_table(clin_file, index_col=0, low_memory=False)
sample_freeze = pd.read_csv(sample_freeze_file)

# Subset Y matrix
common_genes = set(mutation_df.columns).intersection(genes)
common_genes = list(common_genes.intersection(rnaseq_df.columns))

y = mutation_df[common_genes]
missing_genes = set(genes).difference(common_genes)

if len(common_genes) != len(genes):
    warnings.warn('All input genes were not found in data. The missing genes '
                  'are {}'.format(missing_genes), category=Warning)

# Incorporate copy number for gene activation/inactivation
if copy_number:
    copy_loss_df = pd.read_table(copy_loss_file, index_col=0)
    copy_gain_df = pd.read_table(copy_gain_file, index_col=0)

    # Load cancer gene classification table
    cancer_genes = pd.read_table(cancer_gene_file)

    y = integrage_copy_number(y=y, cancer_genes_df=cancer_genes,
                              genes=common_genes, loss_df=copy_loss_df,
                              gain_df=copy_gain_df)

# Remove hypermutated samples from consideration
if remove_hyper:
    mut_burden = pd.read_table(mut_burden_file)
    mut_burden = mut_burden[mut_burden['Silent per Mb'] <
                            2 * mut_burden['Silent per Mb'].std()]
    y = y[y.index.isin(mut_burden.Tumor_Sample_ID)]

# Summarize Y matrix
y_sum = y.reset_index().merge(sample_freeze,
                              how='left').set_index('SAMPLE_BARCODE')

count_df = y_sum.groupby('DISEASE').sum()
prop_df = count_df.divide(y_sum['DISEASE'].value_counts(sort=False)
                                          .sort_index(), axis=0)

out_df = prop_df.merge(count_df, left_index=True, right_index=True)
out_df.to_csv(gene_output_file, sep='\t')
