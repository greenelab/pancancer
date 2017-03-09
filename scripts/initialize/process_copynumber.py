"""
Gregory Way 2016
PanCancer Classifier
scripts/process_copynumber.py

Usage: Run by initialize.sh to initialize repository or in command line

Output:
Two sample by gene matrices that indicate if the gene has 1) copy number gains
or 2) copy number losses in each sample
"""

import os
import pandas as pd

copy_input_file = os.path.join('data', 'raw', 'pancan_GISTIC_threshold.tsv')
copy_loss_file = os.path.join('data', 'copy_number_loss_status.tsv')
copy_gain_file = os.path.join('data', 'copy_number_gain_status.tsv')
sample_freeze_file = os.path.join('data', 'sample_freeze.tsv')

# Load data
copy_thresh_df = pd.read_table(copy_input_file, index_col=0)
copy_thresh_df.drop(['Locus ID', 'Cytoband'], axis=1, inplace=True)
copy_thresh_df.columns = copy_thresh_df.columns.str[0:15]
sample_freeze_df = pd.read_table(sample_freeze_file)

# Process and subset data
copy_thresh_df = copy_thresh_df.T
copy_thresh_df = copy_thresh_df.loc[sample_freeze_df['SAMPLE_BARCODE']]
copy_thresh_df = copy_thresh_df.fillna(0)
copy_thresh_df = copy_thresh_df.astype(int)

# Thresholded copy number includes 5 values [-2, -1, 0, 1, 2], which correspond
# to "deep loss", "moderate loss", "no change", "moderate gain", and "deep
# gain", respectively. We are being conservative and consider only "deep gains"
# and "deep losses" to define impactful copy number events.
copy_loss_df = copy_thresh_df.replace(to_replace=[1, 2, -1], value=0)
copy_loss_df.replace(to_replace=-2, value=1, inplace=True)
copy_loss_df.to_csv(copy_loss_file, sep='\t')

copy_gain_df = copy_thresh_df.replace(to_replace=[-1, -2, 1], value=0)
copy_gain_df.replace(to_replace=2, value=1, inplace=True)
copy_gain_df.to_csv(copy_gain_file, sep='\t')
