"""
Gregory Way 2016
PanCancer NF1 Classifier
process_copynumber.py

Usage: Run in command line to process thresholded copy number data

Output:
Two sample by gene matrices that indicate if the gene has 1) copy number gains
or 2) copy number losses in each sample
"""

import pandas as pd

copy_thresh_df = pd.read_table('data/raw/pancan_GISTIC_threshold.tsv',
                               index_col=0)

copy_thresh_df.drop(['Locus ID', 'Cytoband'], axis=1, inplace=True)
copy_thresh_df.columns = copy_thresh_df.columns.str[0:15]

# Thresholded copy number includes 5 values [-2, -1, 0, 1, 2], which correspond
# to "deep loss", "moderate loss", "no change", "moderate gain", and "deep
# gain", respectively. We are being conservative and consider only "deep gains"
# and "deep losses" to define impactful copy number events.
copy_loss_df = copy_thresh_df.replace(to_replace=[1, 2, -1], value=0)
copy_loss_df.replace(to_replace=-2, value=1, inplace=True)
copy_loss_df.T.to_csv('data/copy_number_loss_status.tsv', sep='\t')

copy_gain_df = copy_thresh_df.replace(to_replace=[-1, -2, 1], value=0)
copy_gain_df.replace(to_replace=2, value=1, inplace=True)
copy_gain_df.T.to_csv('data/copy_number_gain_status.tsv', sep='\t')
