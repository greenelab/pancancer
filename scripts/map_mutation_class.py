"""
Gregory Way 2017
PanCancer Classifier
scripts/map_mutation_class.py

Merge per sample classifier scores with mutation types present in each sample

Usage: Run in command line with required command argument:

        python map_mutation_class.py --scores $prediction_file --genes $genes

prediction_file is a string and genes is a comma sep string of gene symbols

Output:
.tsv file of classifier scores according to mutation class
"""

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scores',
                    help='string of the location of classifier scores')
parser.add_argument('-g', '--genes',
                    help='string of the genes to extract')
args = parser.parse_args()

# Load command arguments
prediction_file = os.path.join(args.scores, 'classifier_decisions.tsv')
genes = args.genes.split(',')
out_dir = os.path.join(os.path.dirname(prediction_file), 'tables')

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

out_file = os.path.join(out_dir, 'mutation_classification_scores.tsv')
raw_mut_file = os.path.join('data', 'raw', 'mc3.v0.2.8.PUBLIC.maf')

pred_df = pd.read_table(prediction_file, index_col=0)
mut_df = pd.read_table(raw_mut_file)

# Process mutation file
mut_df = mut_df.assign(ID=mut_df.Tumor_Sample_Barcode.str.slice(start=0,
                                                                stop=15))
sub_mut_df = mut_df[mut_df['Hugo_Symbol'].isin(genes)]
sub_mut_df = sub_mut_df[['ID', 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'HGVSc',
                         'HGVSp', 'Variant_Classification', ]]

mapped_mutation_df = pred_df.merge(sub_mut_df, left_index=True,
                                   right_on='ID', how='outer')
mapped_mutation_df.to_csv(out_file, sep='\t')
