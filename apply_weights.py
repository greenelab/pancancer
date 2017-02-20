"""
Gregory Way 2017
PanCancer Classifier
apply_weights.py

Usage: Run in command line with required command argument:

        python apply_weights.py --classifier $summary_file

Where summary_file is a string. There are also optional arguments:

    --copy_number       store_true: supplement copy number to define Y

Output:
.tsv file of classifier scores and other covariate info for plotting
"""

import os
import argparse
import warnings
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from tcga_util import integrage_copy_number

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--classifier',
                    help='folder location of classifier file')
parser.add_argument('-u', '--copy_number', action='store_true',
                    help='Supplement Y matrix with copy number events')
args = parser.parse_args()

# Load command arguments
classifier_file = os.path.join(args.classifier, 'classifier_summary.txt')
copy_number = args.copy_number

# Load Constants
rnaseq_file = os.path.join('data', 'pancan_rnaseq_freeze.tsv')
sample_file = os.path.join('data', 'sampleset_freeze_version3.csv')
cancer_gene_file = os.path.join('data', 'vogelstein_cancergenes.tsv')
copy_loss_file = os.path.join('data', 'copy_number_loss_status.tsv')
copy_gain_file = os.path.join('data', 'copy_number_gain_status.tsv')
mutation_burden = os.path.join('ddr', 'data', 'mutation-load.txt')
mut_file = os.path.join('data', 'pancan_mutation_freeze.tsv')

# Generate filenames to save output plots
output_base_file = os.path.dirname(classifier_file)
output_dec_file = os.path.join(output_base_file, 'classifier_decisions.tsv')

# Load Data
rnaseq_df = pd.read_table(rnaseq_file, index_col=0)
mutation_df = pd.read_table(mut_file, index_col=0)
sample_info = pd.read_csv(sample_file, index_col=0)
mut_burden = pd.read_table(mutation_burden)

# Summarize data based on classifier summary file
with open(classifier_file) as class_fh:
    for line in class_fh:
        line = line.strip().split('\t')
        if line[0] == 'Genes:':
            genes = line[1:]
        if line[0] == 'Tissues:':
            diseases = line[1:]
        if line[0] == 'Coefficients:':
            coef_df = pd.read_table(line[1])

# Subset matrix that indicates mutation status (y)
common_genes = set(mutation_df.columns).intersection(genes)
common_genes = list(common_genes.intersection(rnaseq_df.columns))

mut_subset_df = mutation_df[common_genes]
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

    mut_subset_df = integrage_copy_number(y=mut_subset_df,
                                          cancer_genes_df=cancer_genes,
                                          genes=common_genes,
                                          loss_df=copy_loss_df,
                                          gain_df=copy_gain_df)

# Add covariate info to y_matrix
mut_subset_df = mut_subset_df.assign(total_status=mut_subset_df.max(axis=1))
mut_subset_df = mut_subset_df.reset_index().merge(sample_info, how='left')\
                                           .set_index('SAMPLE_BARCODE')
y_burden_matrix = mut_burden.merge(pd.DataFrame(mut_subset_df.total_status),
                                   right_index=True,
                                   left_on='Tumor_Sample_ID')\
                            .set_index('Tumor_Sample_ID')

# Add covariate information and extract Y DataFrame
covar_df = pd.get_dummies(y_burden_matrix['cohort']).astype(int)
covar_df.index = y_burden_matrix.index
covar_df = covar_df.merge(y_burden_matrix, right_index=True, left_index=True)
covar_df.index = y_burden_matrix.index
covar_df = covar_df.drop(['cohort', 'Patient_ID', 'Non-silent per Mb'], axis=1)
y_df = covar_df.total_status

# Subset x matrix to MAD genes, scale expression, and add covariate info
x_df = rnaseq_df.ix[y_df.index, :]
scaled_fit = StandardScaler().fit(x_df)
x_df_update = pd.DataFrame(scaled_fit.transform(x_df),
                           columns=x_df.columns)
x_df_update.index = x_df.index
x_df = x_df_update.merge(covar_df, left_index=True, right_index=True)\
           .drop('total_status', axis=1)

# Reorder x matrix to the same features as the weight coefficients
x_df = x_df[coef_df['feature']]

# Apply a logit transform [y = 1/(1+e^(-wX))] to output probabilities
apply_weights = pd.DataFrame(coef_df['weight'])
apply_weights.index = coef_df.feature
result = apply_weights.T.dot(x_df.T)
result = 1 / (1 + np.exp(-1 * result))

# Investigate decisions
final_pred = y_burden_matrix.join(result.T)
final_pred = final_pred.merge(mut_subset_df.drop('total_status', axis=1),
                              left_index=True, right_index=True)

# Add information (hypermutated samples and if they were used to train model)
final_pred = final_pred.assign(include=0)
final_pred = final_pred.assign(hypermutated=1)
final_pred.loc[(final_pred.cohort.isin(diseases)) &
               (final_pred.hypermutated == 0), 'include'] = 1
final_pred.loc[final_pred['Silent per Mb'] < 2 *
               final_pred['Silent per Mb'].std(), 'hypermutated'] = 0
final_pred.to_csv(output_dec_file, sep='\t')
