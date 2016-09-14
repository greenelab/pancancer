"""
Gregory Way 2016
Modified from https://github.com/cognoma/machine-learning/
PanCancer NF1 Classifier
classifier.py

Usage: Run in command line with required command argument:

        python classifier.py --genes $GENES

Where GENES is a comma separated string. There are also optional arguments:

    --tissues           comma separated string of tissues to build classifier
    --drop              drop the input genes from the X matrix
    --filter_count      int of low count of mutation to include tissue
    --filter_prop       float of low proportion of mutated samples per tissue
    --num_features      int of number of genes to include in classifier
    --alphas            comma separated string of alphas to test in pipeline
    --l1_ratios         comma separated string of l1 parameters to test
    --alt_genes         comma separated string of alternative genes to test

Output:
ROC curves, AUROC across tissues, and classifier coefficients
"""

import warnings

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn import grid_search
from sklearn.linear_model import SGDClassifier
from sklearn.cross_validation import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest
from statsmodels.robust.scale import mad

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genes',
                    help='Comma separated string of HUGO gene symbols')
parser.add_argument('-t', '--tissues', default='Auto',
                    help='Comma separated string of TCGA tissue acronyms. If'
                         'no arguments are passed, filtering will default to '
                         'options given in --filter_count and --filter_prop.')
parser.add_argument('-d', '--drop', action='store_true',
                    help='Decision to drop input genes from X matrix')
parser.add_argument('-c', '--filter_count', default=15,
                    help='Minimum number of mutations in tissue to include')
parser.add_argument('-p', '--filter_prop', default=0.05,
                    help='Minimum proportion of positives to include tissue')
parser.add_argument('-f', '--num_features', default=8000,
                    help='Number of MAD genes to include in classifier')
parser.add_argument('-a', '--alphas', default='0.01,0.1,0.15,0.2,0.5,0.8',
                    help='the alphas for parameter sweep')
parser.add_argument('-l', '--l1_ratios', default='0,0.1,0.15,0.18,0.2,0.3',
                    help='the alphas for parameter sweep')
parser.add_argument('-b', '--alt_genes', default='None',
                    help='alternative genes to test classifier performance')
args = parser.parse_args()

# Load command arguments
genes = args.genes.split(',')
tissues = args.tissues.split(',')
drop = args.drop
filter_count = args.filter_count
filter_prop = args.filter_prop
num_features_kept = args.num_features
alphas = [float(x) for x in args.alphas.split(',')]
l1_ratios = [float(x) for x in args.l1_ratios.split(',')]
alt_genes = args.alt_genes.split(',')
warnings.filterwarnings('ignore',
                        message='Changing the shape of non-C contiguous array')


def fs_mad(x, y):
    """
    Get the median absolute deviation (MAD) for each column of x. The function
    is input as the 'score_func' in sklearn.feature_selection.SelectKBest()

    Arguments:
    x - pandas DataFrame of gene expression values (samples by genes)
    y - pandas Series of mutation outcomes (required by SelectKBest())

    Output:
    Per gene median absolute deviation and p values (according to SelectKBest)
    """
    scores = mad(x)
    return scores, np.array([np.NaN]*len(scores))


def grid_scores_to_df(grid_scores):
    """
    Convert a sklearn.grid_search.GridSearchCV.grid_scores_ attribute to
    a tidy pandas DataFrame where each row is a hyperparam-fold combinatination

    Arguments:
    grid_scores - a list of named tuples with cross valiation results

    Output:
    pandas DataFrame of cross validation scores across each fold
    """
    rows = list()
    for grid_score in grid_scores:
        for fold, score in enumerate(grid_score.cv_validation_scores):
            row = grid_score.parameters.copy()
            row['fold'] = fold
            row['score'] = score
            rows.append(row)
    df = pd.DataFrame(rows)
    return df


def get_threshold_metrics(y_true, y_pred, tissue='all'):
    """
    Retrieve true/false positive rates and auroc for classification predictions

    Arguments:
    y_true - an array of gold standard mutation status
    y_pred - an array of predicted mutation status
    tissue - a string that includes the corresponding TCGA study acronym

    Output:
    A dictionary storing AUROC, a pandas dataframe of ROC data, and tissue
    """
    roc_columns = ['fpr', 'tpr', 'threshold']
    roc_items = zip(roc_columns, roc_curve(y_true, y_pred))
    roc_df = pd.DataFrame.from_items(roc_items)
    auroc = roc_auc_score(y_true, y_pred, average='weighted')
    return {'auroc': auroc, 'roc_df': roc_df, 'tissue': tissue}

# Load Data
rnaseq_df = pd.read_table('data/rnaseq_data.tsv', index_col=0)
mutation_df = pd.read_table('data/mutation_table.tsv', index_col=0)
clinical_df = pd.read_table('data/clinical_data.tsv', index_col=0,
                            low_memory=False)

# Subset data
rnaseq_df.drop('SLC35E2', axis=0, inplace=True)

# Generate file names for output
base_fh = 'tissues_' + args.tissues.replace(',', '_') + '_genes_' + \
          args.genes.replace(',', '_')
if drop:
    base_fh = base_fh + '_DROP_GENE_INPUT_'
    for gene in genes:
        rnaseq_df.drop(gene, axis=0, inplace=True)

mutation_heatmap_fh = 'figures/mutation_heatmap_' + base_fh + '.png'
cv_plot_fh = 'figures/cv_plot_' + base_fh + '.png'
cv_heatmap_fh = 'figures/cv_heatmap_' + base_fh + '.png'
full_roc_fh = 'figures/all_tissue_roc_' + base_fh + '.png'
tissue_roc_fh = 'figures/tissue/classifier_' + base_fh
tissue_summary_fh = 'figures/tissue_summary_' + base_fh + '.png'
classifier_fh = 'classifiers/log_coef_' + base_fh + '.tsv'
alt_gene_tissue_fh = 'figures/alt_gene_' + args.alt_genes.replace(',', '_') + \
                     '_classifier_' + base_fh + '.png'

# Construct data for classifier
y = mutation_df[genes]
y = y.assign(disease=clinical_df['acronym'])

# Filter tissues
mut_count = y.groupby('disease').sum().sum(axis=1).sort_index()
prop = mut_count.divide(y.disease.value_counts(sort=False).sort_index())

if tissues[0] == 'Auto':
    filter_tissue = (mut_count > filter_count) & (prop > filter_prop)
    tissues = filter_tissue.index[filter_tissue].tolist()

y_df = y[y.disease.isin(tissues)].max(axis=1)
x_df = rnaseq_df.ix[:, y_df.index].T
clinical_sub = clinical_df.ix[y_df.index]

# Summary heatmap of mutation counts
mut_df = mutation_df.join(clinical_df['acronym'], how='inner')
mut_df = mut_df.assign(TOTAL=mut_df[genes].max(axis=1))
mutation_heatmap = mut_df[genes + ['acronym']].groupby('acronym').sum()

# Show mutation frequency heatmap
heatmap_df = mutation_heatmap.divide(mut_df.acronym.value_counts(sort=False)
                                     .sort_index(), axis=0)
ax = sns.heatmap(heatmap_df, annot=True)
ax.set(xlabel='gene', ylabel='study')
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(mutation_heatmap_fh)
plt.close()

# Build classifier pipeline
strat = clinical_df.ix[y_df.index]['acronym'].str.cat(y_df.astype(str))
x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.1,
                                                    random_state=0,
                                                    stratify=strat)
feature_select = SelectKBest(fs_mad, k=num_features_kept)


param_fixed = {'loss': 'log', 'penalty': 'elasticnet'}
param_grid = {'alpha': alphas, 'l1_ratio': l1_ratios}

# Include loss='log' in param_grid doesn't work with pipeline somehow
clf = SGDClassifier(random_state=0, class_weight='balanced',
                    loss=param_fixed['loss'], penalty=param_fixed['penalty'])

clf_grid = grid_search.GridSearchCV(estimator=clf, param_grid=param_grid,
                                    n_jobs=-1, scoring='roc_auc')
pipeline = make_pipeline(feature_select, StandardScaler(), clf_grid)

# Fit the model
pipeline.fit(X=x_train, y=y_train)
best_clf = clf_grid.best_estimator_
feature_mask = feature_select.get_support()
cv_score_df = grid_scores_to_df(clf_grid.grid_scores_)

# Visualize cross validation performance
facet_grid = sns.factorplot(x='l1_ratio', y='score', col='alpha',
                            data=cv_score_df, kind='violin', size=4, aspect=1)
facet_grid.set_ylabels('AUROC')
plt.tight_layout()
plt.savefig(cv_plot_fh)
plt.close()

# Cross-validated performance heatmap
cv_score_mat = pd.pivot_table(cv_score_df, values='score', index='l1_ratio',
                              columns='alpha')
ax = sns.heatmap(cv_score_mat, annot=True, fmt='.2%')
ax.set_xlabel('Regularization strength multiplier (alpha)')
ax.set_ylabel('Elastic net mixing parameter (l1_ratio)')
plt.tight_layout()
plt.savefig(cv_heatmap_fh)
plt.close()

# Get predictions
y_pred_train = pipeline.decision_function(x_train)
y_pred_test = pipeline.decision_function(x_test)
metrics_train = get_threshold_metrics(y_train, y_pred_train)
metrics_test = get_threshold_metrics(y_test, y_pred_test)

# Plot ROC
plt.figure()
for label, metrics in ('Training', metrics_train), ('Testing', metrics_test):
    roc_df = metrics['roc_df']
    plt.plot(roc_df.fpr, roc_df.tpr,
             label='{} (AUROC = {:.1%})'.format(label, metrics['auroc']))
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(full_roc_fh)
plt.close()

# Plot tissue specific performance
tissue_metrics = {}
for tissue in clinical_sub.acronym.unique():
    sample_sub = clinical_sub[clinical_sub['acronym'] == tissue].index.values

    y_tissue_train = y_train[y_train.index.isin(sample_sub)]
    y_tissue_pred_train = y_pred_train[y_train.index.isin(sample_sub)]
    y_tissue_test = y_test[y_test.index.isin(sample_sub)]
    y_tissue_pred_test = y_pred_test[y_test.index.isin(sample_sub)]

    metrics_train = get_threshold_metrics(y_tissue_train, y_tissue_pred_train,
                                          tissue=tissue)
    metrics_test = get_threshold_metrics(y_tissue_test, y_tissue_pred_test,
                                         tissue=tissue)
    tissue_metrics[tissue] = [metrics_train, metrics_test]

tissue_auroc = {}
for tissue, metrics_val in tissue_metrics.items():
    met_train, met_test = metrics_val
    tissue_roc_sub_fh = tissue_roc_fh + '_pred_' + tissue + '.png'
    plt.figure()
    auroc = []
    for label, metrics in ('Training', met_train), ('Testing', met_test):
        roc_df = metrics['roc_df']
        auroc.append(metrics['auroc'])
        plt.plot(roc_df.fpr, roc_df.tpr,
                 label='{} (AUROC = {:.1%})'.format(label, metrics['auroc']))
    tissue_auroc[tissue] = auroc
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Predicting Mutation in {}'.format(tissue))
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(tissue_roc_sub_fh)
    plt.close()

tissue_results = pd.DataFrame(tissue_auroc, index=['Train', 'Test']).T
tissue_results = tissue_results.sort_values('Test', ascending=False)
ax = tissue_results.plot(kind='bar', title='Tissue Specific Performance')
ax.set_ylabel('AUROC')
plt.tight_layout()
plt.savefig(tissue_summary_fh)
plt.close()

# Save classifier coefficients
coef_df = pd.DataFrame(best_clf.coef_.transpose(),
                       index=x_train.columns[feature_mask], columns=['weight'])
coef_df['abs'] = coef_df['weight'].abs()
coef_df = coef_df.sort_values('abs', ascending=False)
coef_df.to_csv(classifier_fh, sep='\t')

# Apply classifier to predict alternative genes
if alt_genes[0] is not 'None':
    # Classifying NF1 mutations
    y_alt = pd.DataFrame(mutation_df[alt_genes])

    # Append tissue id
    y_alt = y_alt.assign(disease=clinical_df['acronym'])

    # Filter tissues
    mut_co = y_alt.groupby('disease').sum().sum(axis=1).sort_index()
    prop = mut_co.divide(y_alt.disease.value_counts(sort=False).sort_index())

    filter_tissue = (mut_co > filter_count) & (prop > filter_prop)
    tissues = filter_tissue.index[filter_tissue].tolist()

    # Subset data
    y_alt_df = y_alt[y_alt.disease.isin(tissues)].max(axis=1)
    x_alt_df = rnaseq_df.ix[coef_df.index, y_alt_df.index].T
    clin_sub = clinical_df.ix[y_alt_df.index]

    validation_metrics = {}
    for tissue in tissues:
        sample_sub = clin_sub[clin_sub['acronym'] == tissue].index.values
        y_sub = y_alt_df[sample_sub]
        x_sub = x_alt_df.ix[sample_sub]
        y_pred_nf1 = best_clf.predict(x_sub)
        alt_metrics = get_threshold_metrics(y_sub, y_pred_nf1, tissue=tissue)
        validation_metrics[tissue] = alt_metrics

    tissue_auroc = {}
    for tissue, metrics_val in validation_metrics.items():
        tissue_auroc[tissue] = metrics_val['auroc']

    tissue_results = pd.DataFrame(tissue_auroc, index=['Test']).T
    tissue_results = tissue_results.sort_values('Test', ascending=False)
    ax = tissue_results.plot(kind='bar', title='Alt Gene Performance')
    ax.set_ylim([0, 1])
    ax.set_ylabel('AUROC')
    plt.tight_layout()
    plt.savefig(alt_gene_tissue_fh)
    plt.close()
