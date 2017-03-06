"""
Gregory Way 2017
PanCancer Classifier
scripts/visualize_decisions.py

Usage: Run in command line with required command argument:

        python visualize_decisions.py --scores $prediction_file

Where prediction_file is a string pointing to the location of sample scores

Output:
.tsv file of classifier scores and other covariate info for plotting
"""

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_decision_function(subset_df, filename, title):
    """
    plot and save decision function for classifier

    Arguments:
    subset_df - a pandas dataframe storing prediction and covariate info
    filename - a string holding the filename to save the figure
    title - the title to add to the plot
    """
    sns.set(font_scale=1)
    sns.set_style("white")
    plt.figure(figsize=(3.1, 2.5))
    if subset_df[subset_df.total_status == 1].shape[0] > 0:
        ax = sns.kdeplot(subset_df.ix[subset_df.total_status == 1, :].weight,
                         color='red', label='Deficient', shade=True)
    if subset_df[subset_df.total_status == 0].shape[0] > 0:
        ax = sns.kdeplot(subset_df.ix[subset_df.total_status == 0, :].weight,
                         color='blue', label='Wild-Type', shade=True)
    ax.set(xlabel='Probability', ylabel='Density')
    ax.set_xlim(0, 1.1)
    ax.set_title(title)
    sns.despine()
    plt.axvline(x=0.5, color='k', ls='dashed', linewidth=0.7)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scores',
                    help='folder location of classifier scores')
parser.add_argument('-c', '--custom',
                    help='comma separated list of columns to plot')
args = parser.parse_args()

# Load command arguments
prediction_file = os.path.join(args.scores, 'classifier_decisions.tsv')
custom_columns = args.custom.split(',')

out_dir = os.path.join(os.path.dirname(prediction_file), 'figures')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

pred_df = pd.read_table(prediction_file)
diseases = pred_df['cohort'].unique()

# Visualize decision function
plot_decision_function(subset_df=pred_df[pred_df.include == 1],
                       filename=os.path.join(out_dir, 'total_decisions.pdf'),
                       title='Classifier Decision Function')

# Plot disease type specific decision functions
for disease_type in diseases:
    sub_df = pred_df.ix[pred_df.cohort == disease_type]
    d_file = os.path.join(out_dir, 'decision_plot_{}.pdf'.format(disease_type))
    d_title = 'Classifier Decision Function\n{}'.format(disease_type)
    plot_decision_function(subset_df=sub_df, filename=d_file, title=d_title)

# Visualize decision function for hyper mutated tumors
plot_decision_function(subset_df=pred_df[pred_df.hypermutated == 1],
                       filename=os.path.join(out_dir, 'hyper_mutated.pdf'),
                       title='Hypermutated Tumors')

# Visualize decision function for copy number loss tumors
for col in custom_columns:
    plot_file = os.path.join(out_dir, '{}_decision_plot.pdf'.format(col))
    plot_decision_function(subset_df=pred_df[pred_df[col] == 1],
                           filename=plot_file,
                           title=col)
