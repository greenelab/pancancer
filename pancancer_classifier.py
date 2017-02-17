"""
Gregory Way 2017
Heavily modified from https://github.com/cognoma/machine-learning/
PanCancer Classifier
pancancer_classifier.py

Usage: Run in command line with required command argument:

        python pancancer_classifier.py --genes $GENES

Where GENES is a comma separated string. There are also optional arguments:

    --diseases          comma separated string of disease types for classifier
    --folds             number of cross validation folds (defaults to 5)
    --drop              drop the input genes from the X matrix
    --copy_number       optional flag to supplement copy number to define Y
    --filter_count      int of low count of mutation to include tissue
    --filter_prop       float of low proportion of mutated samples per tissue
    --num_features      int of number of genes to include in classifier
    --alphas            comma separated string of alphas to test in pipeline
    --l1_ratios         comma separated string of l1 parameters to test
    --alt_genes         comma separated string of alternative genes to test
    --alt_tissues       comma separated string of alternative tissues to test
    --alt_filter_count  int of low count of mutations to include alt_tissues
    --alt_filter_prop   float of low proportion of mutated samples alt_tissue
    --alt_folder        string of where to save the classifier figures
    --remove_hyper      store_true: remove hypermutated samples

Output:
ROC curves, AUROC across tissues, and classifier coefficients
"""

import os
import warnings

import numpy as np
import pandas as pd
import csv
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import train_test_split, GridSearchCV, \
                                    cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.robust.scale import mad

from tcga_util import integrage_copy_number

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genes',
                    help='Comma separated string of HUGO gene symbols')
parser.add_argument('-t', '--diseases', default='Auto',
                    help='Comma separated string of TCGA tissue acronyms. If'
                         'no arguments are passed, filtering will default to '
                         'options given in --filter_count and --filter_prop.')
parser.add_argument('-f', '--folds', default='5',
                    help='Number of cross validation folds to perform')
parser.add_argument('-d', '--drop', action='store_true',
                    help='Decision to drop input genes from X matrix')
parser.add_argument('-u', '--copy_number', action='store_true',
                    help='Supplement Y matrix with copy number events')
parser.add_argument('-c', '--filter_count', default=15,
                    help='Minimum number of mutations in tissue to include')
parser.add_argument('-p', '--filter_prop', default=0.05,
                    help='Minimum proportion of positives to include tissue')
parser.add_argument('-n', '--num_features', default=8000,
                    help='Number of MAD genes to include in classifier')
parser.add_argument('-a', '--alphas', default='0.1,0.15,0.2,0.5,0.8,1',
                    help='the alphas for parameter sweep')
parser.add_argument('-l', '--l1_ratios', default='0,0.1,0.15,0.18,0.2,0.3',
                    help='the l1 ratios for parameter sweep')
parser.add_argument('-b', '--alt_genes', default='None',
                    help='alternative genes to test classifier performance')
parser.add_argument('-s', '--alt_tissues', default="Auto",
                    help='The alternative tissues to test classifier'
                         'performance')
parser.add_argument('-i', '--alt_filter_count', default=15,
                    help='Minimum number of mutations in tissue to include')
parser.add_argument('-r', '--alt_filter_prop', default=0.05,
                    help='Minimum proportion of positives to include tissue')
parser.add_argument('-o', '--alt_folder', default='Auto',
                    help='Provide an alternative folder to save results')
parser.add_argument('-v', '--remove_hyper', action='store_true',
                    help='Remove hypermutated samples')
args = parser.parse_args()

# Load command arguments
genes = args.genes.split(',')
tissues = args.diseases.split(',')
folds = int(args.folds)
drop = args.drop
copy_number = args.copy_number
filter_count = int(args.filter_count)
filter_prop = float(args.filter_prop)
num_features_kept = args.num_features
alphas = [float(x) for x in args.alphas.split(',')]
l1_ratios = [float(x) for x in args.l1_ratios.split(',')]
alt_genes = args.alt_genes.split(',')
alt_filter_count = int(args.alt_filter_count)
alt_filter_prop = float(args.alt_filter_prop)
alt_tissues = args.alt_tissues.split(',')
alt_folder = args.alt_folder
remove_hyper = args.remove_hyper

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

# Generate file names for output
base_add = 'synapse'
base_folder = '{}__tissues_{}__genes_{}'.format(base_add,
                                                args.diseases.replace(',', '_'),
                                                args.genes.replace(',', '_'))
if drop:
    base_folder = '{}__drop_gene_input'.format(base_folder)
if copy_number:
    base_folder = '{}__copy_number'.format(base_folder)

base_folder = os.path.join('classifiers', base_folder + '__classifier')

if alt_folder != 'Auto':
    base_folder = alt_folder

if not os.path.exists(base_folder):
    os.makedirs(base_folder)
    os.makedirs(os.path.join(base_folder, 'tissue'))
else:
    warnings.warn('Classifier may have already been built!', category=Warning)

cv_heatmap_file = os.path.join(base_folder, 'cv_heatmap.pdf')
full_roc_file = os.path.join(base_folder, 'all_tissue_roc.pdf')
tissue_roc_file = os.path.join(base_folder, 'tissue', 'classifier_')
tissue_summary_file = os.path.join(base_folder, 'tissue_summary.pdf')
classifier_file = os.path.join(base_folder, 'classifier_coefficients.tsv')
alt_gene_tissue_file = os.path.join(base_folder, 'alt_gene_{}_alt_tissues_'
                                    '{}_classifier.pdf'.format(
                                     args.alt_genes.replace(',', '_'),
                                     args.alt_tissues.replace(',', '_')))

# Load Datasets
expr_file = os.path.join('data', 'pancan_rnaseq_freeze.tsv')
mut_file = os.path.join('data', 'pancan_mutation_freeze.tsv')
sample_freeze_file = os.path.join('data', 'sampleset_freeze_version3.csv')

rnaseq_df = pd.read_table(expr_file, index_col=0)
mutation_df = pd.read_table(mut_file, index_col=0)
sample_freeze = pd.read_csv(sample_freeze_file)

# Construct data for classifier
common_genes = set(mutation_df.columns).intersection(genes)
common_genes = list(common_genes.intersection(rnaseq_df.columns))

y = mutation_df[common_genes]
missing_genes = set(genes).difference(common_genes)

if len(common_genes) != len(genes):
    warnings.warn('All input genes were not found in data. The missing genes '
                  'are {}'.format(missing_genes), category=Warning)

if drop:
    rnaseq_df.drop(common_genes, axis=1, inplace=True)

# Incorporate copy number for gene activation/inactivation
if copy_number:
    # Load copy number matrices
    copy_loss_df = pd.read_table(os.path.join('data',
                                              'copy_number_loss_status.tsv'),
                                 index_col=0)
    copy_gain_df = pd.read_table(os.path.join('data',
                                              'copy_number_gain_status.tsv'),
                                 index_col=0)

    # Load cancer gene classification table
    cancer_genes = pd.read_table(os.path.join('data',
                                              'vogelstein_cancergenes.tsv'))

    y = integrage_copy_number(y=y, cancer_genes_df=cancer_genes,
                              genes=common_genes, loss_df=copy_loss_df,
                              gain_df=copy_gain_df)

# Process y matrix
y = y.assign(total_status=y.max(axis=1))
y = y.reset_index().merge(sample_freeze,
                          how='left').set_index('SAMPLE_BARCODE')
group_count_df = y.groupby('DISEASE').sum()
heatmap_df = group_count_df.divide(y['DISEASE'].value_counts(sort=False)
                                               .sort_index(), axis=0)

# Filter tissues
mut_count = group_count_df['total_status']
prop = heatmap_df['total_status']

if tissues[0] == 'Auto':
    filter_tissue = (mut_count > filter_count) & (prop > filter_prop)
    tissues = filter_tissue.index[filter_tissue].tolist()

# Load mutation burden and process covariates
y = y[y.DISEASE.isin(tissues)].total_status
mut_burden = pd.read_table(os.path.join('ddr', 'data',
                                        'mutation-load.txt'))
if remove_hyper:
    mut_burden = mut_burden[mut_burden['Silent per Mb'] <
                            2 * mut_burden['Silent per Mb'].std()]
y_matrix = mut_burden.merge(pd.DataFrame(y), right_index=True,
                            left_on='Tumor_Sample_ID')\
    .set_index('Tumor_Sample_ID')

# Add covariate information
covar = pd.get_dummies(y_matrix['cohort']).astype(int)
covar.index = y_matrix.index
covar = covar.merge(y_matrix, right_index=True, left_index=True)
covar.index = y_matrix.index
covar = covar.drop(['cohort', 'Patient_ID', 'Non-silent per Mb'], axis=1)
y_df = covar.total_status
strat = y_matrix.cohort.str.cat(y_matrix.total_status.astype(str))[y_df.index]

# Subset x matrix to MAD genes and scale
x_df = rnaseq_df.ix[y_df.index, :]
med_dev = pd.DataFrame(mad(x_df), index=x_df.columns)
mad_genes = med_dev.sort_values(by=0, ascending=False)\
                   .iloc[0:num_features_kept].index.tolist()
x_df = x_df.loc[:, mad_genes]
x_df_update = pd.DataFrame(StandardScaler().fit_transform(x_df),
                           columns=x_df.columns)
x_df_update.index = x_df.index
x_df = x_df_update.merge(covar, left_index=True, right_index=True)\
                  .drop('total_status', axis=1)

# Build classifier pipeline
x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.1,
                                                    random_state=0,
                                                    stratify=strat)

clf_parameters = {'classify__loss': ['log'],
                  'classify__penalty': ['elasticnet'],
                  'classify__alpha': alphas, 'classify__l1_ratio': l1_ratios}

estimator = Pipeline(steps=[('classify', SGDClassifier(random_state=0,
                                                       class_weight='balanced',
                                                       loss='log'))])

cv_pipeline = GridSearchCV(estimator=estimator, param_grid=clf_parameters,
                           n_jobs=-1, cv=folds, scoring='roc_auc')
cv_pipeline.fit(X=x_train, y=y_train)

cv_results = pd.concat([pd.DataFrame(cv_pipeline.cv_results_)
                          .drop('params', axis=1),
                        pd.DataFrame.from_records(cv_pipeline
                                                  .cv_results_['params'])],
                       axis=1)

# Cross-validated performance heatmap
cv_score_mat = pd.pivot_table(cv_results, values='mean_test_score',
                              index='classify__l1_ratio',
                              columns='classify__alpha')
ax = sns.heatmap(cv_score_mat, annot=True, fmt='.1%')
ax.set_xlabel('Regularization strength multiplier (alpha)')
ax.set_ylabel('Elastic net mixing parameter (l1_ratio)')
plt.tight_layout()
plt.savefig(cv_heatmap_file, dpi=600, bbox_inches='tight')
plt.close()

# Get predictions
y_pred_train = cv_pipeline.decision_function(x_train)
y_pred_test = cv_pipeline.decision_function(x_test)
metrics_train = get_threshold_metrics(y_train, y_pred_train)
metrics_test = get_threshold_metrics(y_test, y_pred_test)

# Rerun "cross validation" for the best hyperparameter set to define
# cross-validation disease-specific performance. Each sample prediction is
# based on the fold that the sample was in the testing partition
y_cv = cross_val_predict(cv_pipeline.best_estimator_, X=x_train, y=y_train,
                         cv=folds, method='decision_function')
metrics_cv = get_threshold_metrics(y_train, y_cv)

# Plot ROC
sns.set_style("whitegrid")
plt.figure(figsize=(2.7, 2.4))
total_auroc = {}
colors = ['blue', 'green', 'orange']
idx = 0
for label, metrics in [('Training', metrics_train), ('Testing', metrics_test),
                       ('CV', metrics_cv)]:

    roc_df = metrics['roc_df']
    plt.plot(roc_df.fpr, roc_df.tpr,
             label='{} (AUROC = {:.1%})'.format(label, metrics['auroc']),
             linewidth=2, c=colors[idx])
    total_auroc[label] = metrics['auroc']
    idx += 1
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=8)
plt.ylabel('True Positive Rate', fontsize=8)
plt.title('')
plt.tick_params(labelsize=8)
plt.legend(bbox_to_anchor=(0.2, -0.45, 0.7, .202), loc=0, borderaxespad=0.,
           fontsize=7.5)
plt.tight_layout()
plt.savefig(full_roc_file, dpi=600, bbox_inches='tight')
plt.close()

# Tissue specific performance
tissue_metrics = {}
for tissue in tissues:
    sample_sub = y_matrix[y_matrix.cohort == tissue].index.values

    y_tissue_train = y_train[y_train.index.isin(sample_sub)]
    y_tissue_pred_train = y_pred_train[y_train.index.isin(sample_sub)]
    y_tissue_test = y_test[y_test.index.isin(sample_sub)]
    y_tissue_pred_test = y_pred_test[y_test.index.isin(sample_sub)]
    y_tissue_pred_cv = y_cv[y_train.index.isin(sample_sub)]

    metrics_train_tis = get_threshold_metrics(y_tissue_train,
                                              y_tissue_pred_train,
                                              tissue=tissue)
    metrics_test_tis = get_threshold_metrics(y_tissue_test, y_tissue_pred_test,
                                             tissue=tissue)
    metrics_cv_tis = get_threshold_metrics(y_tissue_train, y_tissue_pred_cv,
                                           tissue=tissue)
    tissue_metrics[tissue] = [metrics_train_tis, metrics_test_tis,
                              metrics_cv_tis]

tissue_auroc = {}
for tissue, metrics_val in tissue_metrics.items():
    met_train, met_test, met_cv = metrics_val
    tissue_roc_sub_file = '{}_pred_{}.pdf'.format(tissue_roc_file, tissue)

    plt.figure(figsize=(2.7, 2.4))
    auroc = []
    idx = 0
    for label, metrics in [('Training', met_train), ('Testing', met_test),
                           ('CV', met_cv)]:
        roc_df = metrics['roc_df']
        plt.plot(roc_df.fpr, roc_df.tpr,
                 label='{} (AUROC = {:.1%})'.format(label, metrics['auroc']),
                 linewidth=2, c=colors[idx])
        auroc.append(metrics['auroc'])
        idx += 1
    tissue_auroc[tissue] = auroc

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=8)
    plt.ylabel('True Positive Rate', fontsize=8)
    plt.title('')
    plt.tick_params(labelsize=8)
    plt.legend(bbox_to_anchor=(0.2, -0.45, 0.7, .202), loc=0, borderaxespad=0.,
               fontsize=7.5)
    plt.tight_layout()
    plt.savefig(tissue_roc_sub_file, dpi=600, bbox_inches='tight')
    plt.close()

tissue_results = pd.DataFrame(tissue_auroc, index=['Train', 'Test',
                                                   'Cross Validation']).T
tissue_results = tissue_results.sort_values('Cross Validation',
                                            ascending=False)
ax = tissue_results.plot(kind='bar', title='Disease Specific Performance')
ax.set_ylabel('AUROC')
plt.tight_layout()
plt.savefig(tissue_summary_file, dpi=600, bbox_inches='tight')
plt.close()

# Save classifier coefficients
final_pipeline = cv_pipeline.best_estimator_
final_classifier = final_pipeline.named_steps['classify']

coef_df = pd.DataFrame.from_items([
    ('feature', x_df.columns),
    ('weight', final_classifier.coef_[0])])

coef_df['abs'] = coef_df['weight'].abs()
coef_df = coef_df.sort_values('abs', ascending=False)
coef_df.to_csv(classifier_file, sep='\t')

# Apply the same classifier previously built to predict alternative genes
if alt_genes[0] is not 'None':
    # Classifying alternative mutations
    y_alt = mutation_df[alt_genes]

    # Add copy number info if applicable
    if copy_number:
        y_alt = integrage_copy_number(y=y_alt, cancer_genes_df=cancer_genes,
                                      genes=alt_genes, loss_df=copy_loss_df,
                                      gain_df=copy_gain_df)
    # Append tissue id
    y_alt = y_alt.assign(total_status=y_alt.max(axis=1))
    y_alt = y_alt.reset_index().merge(sample_freeze,
                                      how='left').set_index('SAMPLE_BARCODE')

    # Filter data
    group_count_df = y_alt.groupby('DISEASE').sum()
    group_prop_df = group_count_df.divide(y_alt['DISEASE']
                                          .value_counts(sort=False)
                                          .sort_index(), axis=0)

    mut_co = group_count_df['total_status']
    prop = group_prop_df['total_status']

    if alt_tissues[0] == 'Auto':
        alt_filter_tis = (mut_co > alt_filter_count) & (prop > alt_filter_prop)
        alt_tissues = alt_filter_tis.index[alt_filter_tis].tolist()

    # Subset data
    y_alt_df = y_alt[y_alt.DISEASE.isin(alt_tissues)].total_status

    # Retrieve cross validation performance for predicting alternative genes
    # The classifier used here is the same classifier built using input genes
    y_alt_cv = cross_val_predict(cv_pipeline.best_estimator_,
                                 X=rnaseq_df.ix[y_alt_df.index],
                                 y=y_alt_df, cv=folds,
                                 method='decision_function')

    alt_metrics_cv = get_threshold_metrics(y_alt_df, y_alt_cv)

    validation_metrics = {}
    validation_x_type = {}
    for tissue in alt_tissues:
        sample_tissues = y_alt[y_alt.DISEASE == tissue].index.values

        # Subset full data if it has not been trained on
        if tissue not in tissues:
            x_sub = rnaseq_df.ix[sample_tissues]
            y_sub = y_alt_df[sample_tissues]
            validation_x_type[tissue] = 'Full'

        # Only subset to the holdout set if data was trained on
        else:
            x_sub = x_test.ix[x_test.index.isin(sample_tissues)]
            y_sub = y_test[y_test.index.isin(sample_tissues)]
            validation_x_type[tissue] = 'Holdout'

        y_pred_alt = cv_pipeline.decision_function(x_sub)
        y_pred_alt_cv = y_alt_cv[y_alt_df.index.isin(y_sub.index)]

        alt_metrics_tis = get_threshold_metrics(y_sub, y_pred_alt,
                                                tissue=tissue)
        alt_metrics_tis_cv = get_threshold_metrics(y_sub, y_pred_alt_cv,
                                                   tissue=tissue)
        validation_metrics[tissue] = [alt_metrics_tis, alt_metrics_tis_cv]

    alt_tissue_auroc = {}
    for tissue, metrics_val in validation_metrics.items():
        met_test, met_cv = metrics_val
        alt_tissue_auroc[tissue] = [met_test['auroc'], met_cv['auroc']]

    alt_tissue_results = pd.DataFrame(alt_tissue_auroc,
                                      index=['Test', 'Cross Validation']).T
    alt_tissue_results = alt_tissue_results.sort_values('Cross Validation',
                                                        ascending=False)
    ax = alt_tissue_results.plot(kind='bar', title='Alt Gene Performance')
    ax.set_ylim([0, 1])
    ax.set_ylabel('AUROC')
    plt.tight_layout()
    plt.savefig(alt_gene_tissue_file)
    plt.close()

# Write a summary for the inputs and outputs of the classifier
with open(os.path.join(base_folder, 'classifier_summary.txt'), 'w') as sum_fh:
    summarywriter = csv.writer(sum_fh, delimiter='\t')

    # Summarize parameters
    summarywriter.writerow(['Parameters:'])
    summarywriter.writerow(['Genes:'] + genes)
    summarywriter.writerow(['Tissues:'] + tissues)
    summarywriter.writerow(['Alternative Genes:'] + alt_genes)
    summarywriter.writerow(['Alternative Tissues:'] + alt_tissues)
    summarywriter.writerow(['Number of Features:', str(x_df.shape[1])])
    summarywriter.writerow(['Drop Gene:', drop])
    summarywriter.writerow(['Copy Number:', copy_number])
    summarywriter.writerow(['Alphas:'] + alphas)
    summarywriter.writerow(['L1_ratios:'] + l1_ratios)
    summarywriter.writerow(['Hypermutated Removed:', str(remove_hyper)])
    summarywriter.writerow([])

    # Summaryize results
    summarywriter.writerow(['Results:'])
    summarywriter.writerow(['Optimal Alpha:',
                            str(cv_pipeline.best_params_['classify__alpha'])])
    summarywriter.writerow(['Optimal L1:', str(cv_pipeline.best_params_
                                               ['classify__l1_ratio'])])
    summarywriter.writerow(['Coefficients:', classifier_file])
    summarywriter.writerow(['Training AUROC:', metrics_train['auroc']])
    summarywriter.writerow(['Testing AUROC:', metrics_test['auroc']])
    summarywriter.writerow(['Cross Validation AUROC', metrics_cv['auroc']])
    summarywriter.writerow(['Disease specific performance:'])
    for tissue, auroc in tissue_auroc.items():
        summarywriter.writerow(['', tissue, 'Training AUROC:', auroc[0],
                                'Testing AUROC:', auroc[1],
                                'Cross Validation AUROC:', auroc[2]])
    if alt_genes[0] is not 'None':
        summarywriter.writerow(['Alternate gene performance:'] + alt_genes)
        summarywriter.writerow(['Alternative gene AUROC:',
                                str(alt_metrics_cv['auroc'])])
        for alt_tissue, alt_auroc in alt_tissue_auroc.items():
            summarywriter.writerow(['', alt_tissue, 'Testing AUROC:',
                                    alt_auroc[0], 'Cross validation AUROC',
                                    alt_auroc[1], 'Data:',
                                    validation_x_type[alt_tissue]])
