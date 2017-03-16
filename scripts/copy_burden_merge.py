"""
Gregory Way 2017
PanCancer Classifier
scripts/copy_burden_merge.py

Merge per sample classifier scores with segment based scores

Usage: Run in command line with required command argument:

        python scripts/copy_burden_merge.py --classifier_folder

classifier_folder is a string pointing to the location of the classifier data

Output:
.tsv file of classifier scores merged with segment based copy number scores
"""

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--classifier_folder',
                    help='string of the location of classifier data')
args = parser.parse_args()

# Load command arguments
pred_fild = os.path.join(args.classifier_folder, 'classifier_decisions.tsv')
burden_file = os.path.join('data', 'seg_based_scores.tsv')
out_file = os.path.join(os.path.dirname(pred_fild), 'tables',
                        'copy_burden_predictions.tsv')

# Load and process data
copy_burden_df = pd.read_table(burden_file)
classifier_df = pd.read_table(pred_fild, index_col=0)

combined_df = classifier_df.merge(copy_burden_df, left_index=True,
                                  right_on='Sample')
combined_df.index = combined_df['Sample']
combined_df.to_csv(out_file, sep='\t')
