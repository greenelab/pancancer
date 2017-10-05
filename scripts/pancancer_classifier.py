"""
Gregory Way 2017
PanCancer Classifier
scripts/pancancer_classifier.py

Usage: Run in command line with required command argument:

        python pancancer_classifier.py --genes $GENES

Where GENES is a comma separated string. There are also optional arguments:

    --diseases          comma separated string of disease types for classifier
                            default: Auto (will pick diseases from filter args)
    --folds             number of cross validation folds
                            default: 5
    --drop              drop the input genes from the X matrix
                            default: False if flag omitted
    --copy_number       optional flag to supplement copy number to define Y
                            default: False if flag omitted
    --filter_count      int of low count of mutation to include disease
                            default: 15
    --filter_prop       float of low proportion of mutated samples per disease
                            default: 0.05
    --num_features      int of number of genes to include in classifier
                            default: 8000
    --alphas            comma separated string of alphas to test in pipeline
                            default: '0.1,0.15,0.2,0.5,0.8,1'
    --l1_ratios         comma separated string of l1 parameters to test
                            default: '0,0.1,0.15,0.18,0.2,0.3'
    --alt_genes         comma separated string of alternative genes to test
                            default: None
    --alt_diseases      comma separated string of alternative diseases to test
                            default: Auto
    --alt_filter_count  int of low count of mutations to include alt_diseases
                            default: 15
    --alt_filter_prop   float of low proportion of mutated samples alt_disease
                            default: 0.05
    --alt_folder        string of where to save the classifier figures
                            default: Auto
    --remove_hyper      store_true: remove hypermutated samples
                            default: False if flag omitted
    --keep_intermediate store_true: keep intermediate roc curve items
                            default: False if flag omitted
    --x_matrix          string of which feature matrix to use
                            default: raw
    --y_matrix          string of which status matrix to use
                            default: default

Output:
ROC curves, AUROC across diseases, and classifier coefficients
"""

import os
import sys
import warnings

import pandas as pd
import csv
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import train_test_split, cross_val_predict
from dask_searchcv import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.robust.scale import mad

sys.path.insert(0, os.path.join('scripts', 'util'))
from tcga_util import get_threshold_metrics, integrate_copy_number

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genes',
                    help='Comma separated string of HUGO gene symbols')
parser.add_argument('-t', '--diseases', default='Auto',
                    help='Comma separated string of TCGA disease acronyms. If'
                         'no arguments are passed, filtering will default to '
                         'options given in --filter_count and --filter_prop.')
parser.add_argument('-f', '--folds', default='5',
                    help='Number of cross validation folds to perform')
parser.add_argument('-d', '--drop', action='store_true',
                    help='Decision to drop input genes from X matrix')
parser.add_argument('-u', '--copy_number', action='store_true',
                    help='Supplement Y matrix with copy number events')
parser.add_argument('-c', '--filter_count', default=15,
                    help='Minimum number of mutations in diseases to include')
parser.add_argument('-p', '--filter_prop', default=0.05,
                    help='Minimum proportion of positives to include disease')
parser.add_argument('-n', '--num_features', default=8000,
                    help='Number of MAD genes to include in classifier')
parser.add_argument('-a', '--alphas', default='0.1,0.15,0.2,0.5,0.8,1',
                    help='the alphas for parameter sweep')
parser.add_argument('-l', '--l1_ratios', default='0,0.1,0.15,0.18,0.2,0.3',
                    help='the l1 ratios for parameter sweep')
parser.add_argument('-b', '--alt_genes', default='None',
                    help='alternative genes to test classifier performance')
parser.add_argument('-s', '--alt_diseases', default="Auto",
                    help='The alternative diseases to test performance')
parser.add_argument('-i', '--alt_filter_count', default=15,
                    help='Minimum number of mutations in disease to include')
parser.add_argument('-r', '--alt_filter_prop', default=0.05,
                    help='Minimum proportion of positives to include disease')
parser.add_argument('-o', '--alt_folder', default='Auto',
                    help='Provide an alternative folder to save results')
parser.add_argument('-v', '--remove_hyper', action='store_true',
                    help='Remove hypermutated samples')
parser.add_argument('-k', '--keep_intermediate', action='store_true',
                    help='Keep intermediate ROC values for plotting')
parser.add_argument('-x', '--x_matrix', default='raw',
                    help='Filename of features to use in model')
parser.add_argument('-y', '--y_matrix', default='default',
                    help='Either `default` or `xena` based on underlying data')
args = parser.parse_args()

# Load command arguments
genes = args.genes.split(',')
diseases = args.diseases.split(',')
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
alt_diseases = args.alt_diseases.split(',')
alt_folder = args.alt_folder
remove_hyper = args.remove_hyper
keep_inter = args.keep_intermediate
x_matrix = args.x_matrix
y_matrix = args.y_matrix

warnings.filterwarnings('ignore',
                        message='Changing the shape of non-C contiguous array')

# Generate file names for output
genes_folder = args.genes.replace(',', '_')
base_folder = os.path.join('classifiers', genes_folder)

if alt_folder != 'Auto':
    base_folder = alt_folder

if not os.path.exists(base_folder):
    os.makedirs(base_folder)
else:
    warnings.warn('Classifier may have already been built! Classifier results'
                  ' will be overwritten!', category=Warning)

disease_folder = os.path.join(base_folder, 'disease')
if not os.path.exists(disease_folder):
    os.makedirs(disease_folder)

count_table_file = os.path.join(base_folder, 'summary_counts.csv')
cv_heatmap_file = os.path.join(base_folder, 'cv_heatmap.svg')
full_roc_file = os.path.join(base_folder, 'all_disease_roc.svg')
full_pr_file = os.path.join(base_folder, 'all_disease_pr.svg')
disease_roc_file = os.path.join(base_folder, 'disease', 'classifier_roc_')
disease_pr_file = os.path.join(base_folder, 'disease', 'classifier_pr_')
dis_summary_auroc_file = os.path.join(base_folder, 'disease_auroc.svg')
dis_summary_aupr_file = os.path.join(base_folder, 'disease_aupr.svg')
classifier_file = os.path.join(base_folder, 'classifier_coefficients.tsv')
roc_results_file = os.path.join(base_folder, 'pancan_roc_results.tsv')

alt_gene_base = 'alt_gene_{}_alt_disease_{}'.format(
                args.alt_genes.replace(',', '_'),
                args.alt_diseases.replace(',', '_'))
alt_count_table_file = os.path.join(base_folder, 'alt_summary_counts.csv')
alt_gene_auroc_file = os.path.join(base_folder,
                                   '{}_auroc_bar.svg'.format(alt_gene_base))
alt_gene_aupr_file = os.path.join(base_folder,
                                  '{}_aupr_bar.svg'.format(alt_gene_base))
alt_gene_summary_file = os.path.join(base_folder,
                                     '{}_summary.tsv'.format(alt_gene_base))

# Load Datasets
if x_matrix == 'raw':
    expr_file = os.path.join('data', 'pancan_rnaseq_freeze.tsv')
else:
    expr_file = x_matrix

if y_matrix == 'default':
    mut_file = os.path.join('data', 'pancan_mutation_freeze.tsv')
    mut_burden_file = os.path.join('data', 'mutation_burden_freeze.tsv')
elif y_matrix == 'xena':
    mut_file = 'https://github.com/greenelab/tybalt/raw/928804ffd3bb3f9d5559796b2221500c303ed92c/data/pancan_mutation.tsv.gz'
    mut_burden_file = 'https://github.com/greenelab/tybalt/raw/87496e23447a06904bf9c07c389584147b87bd65/data/pancan_mutation_burden.tsv'

sample_freeze_file = os.path.join('data', 'sample_freeze.tsv')

rnaseq_full_df = pd.read_table(expr_file, index_col=0)
mutation_df = pd.read_table(mut_file, index_col=0)
sample_freeze = pd.read_table(sample_freeze_file, index_col=0)
mut_burden = pd.read_table(mut_burden_file)

# Construct data for classifier
common_genes = set(mutation_df.columns).intersection(genes)
if x_matrix == 'raw':
    common_genes = list(common_genes.intersection(rnaseq_full_df.columns))
else:
    common_genes = list(common_genes)

y = mutation_df[common_genes]
missing_genes = set(genes).difference(common_genes)

if len(common_genes) != len(genes):
    warnings.warn('All input genes were not found in data. The missing genes '
                  'are {}'.format(missing_genes), category=Warning)

if drop:
    if x_matrix == 'raw':
        rnaseq_full_df.drop(common_genes, axis=1, inplace=True)

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

    y = integrate_copy_number(y=y, cancer_genes_df=cancer_genes,
                              genes=common_genes, loss_df=copy_loss_df,
                              gain_df=copy_gain_df)

# Process y matrix
y = y.assign(total_status=y.max(axis=1))
y = y.reset_index().merge(sample_freeze,
                          how='left').set_index('SAMPLE_BARCODE')
count_df = y.groupby('DISEASE').sum()
prop_df = count_df.divide(y['DISEASE'].value_counts(sort=False).sort_index(),
                          axis=0)

count_table = count_df.merge(prop_df, left_index=True, right_index=True)
count_table.to_csv(count_table_file)

# Filter diseases
mut_count = count_df['total_status']
prop = prop_df['total_status']

if diseases[0] == 'Auto':
    filter_disease = (mut_count > filter_count) & (prop > filter_prop)
    diseases = filter_disease.index[filter_disease].tolist()

# Load mutation burden and process covariates
y_df = y[y.DISEASE.isin(diseases)].total_status
common_samples = list(set(y_df.index) & set(rnaseq_full_df.index))
y_df = y_df.loc[common_samples]
rnaseq_df = rnaseq_full_df.loc[y_df.index, :]
if remove_hyper:
    burden_filter = mut_burden['log10_mut'] < 5 * mut_burden['log10_mut'].std()
    mut_burden = mut_burden[burden_filter]

y_matrix = mut_burden.merge(pd.DataFrame(y_df), right_index=True,
                            left_on='SAMPLE_BARCODE')\
    .set_index('SAMPLE_BARCODE')

# Add covariate information
y_sub = y.loc[y_matrix.index]['DISEASE']
covar_dummy = pd.get_dummies(sample_freeze['DISEASE']).astype(int)
covar_dummy.index = sample_freeze['SAMPLE_BARCODE']
covar = covar_dummy.merge(y_matrix, right_index=True, left_index=True)
covar = covar.drop('total_status', axis=1)

# How cross validation splits will be balanced and stratified
y_df = y_df.loc[y_sub.index]
strat = y_sub.str.cat(y_df.astype(str))
x_df = rnaseq_df.loc[y_df.index, :]

# Subset x matrix to MAD genes and scale
if x_matrix == 'raw':
    med_dev = pd.DataFrame(mad(x_df), index=x_df.columns)
    mad_genes = med_dev.sort_values(by=0, ascending=False)\
                       .iloc[0:num_features_kept].index.tolist()
    x_df = x_df.loc[:, mad_genes]

fitted_scaler = StandardScaler().fit(x_df)
x_df_update = pd.DataFrame(fitted_scaler.transform(x_df),
                           columns=x_df.columns)
x_df_update.index = x_df.index
x_df = x_df_update.merge(covar, left_index=True, right_index=True)

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
plt.savefig(cv_heatmap_file, dpi=600, format='svg', bbox_inches='tight')
plt.close()

# Get predictions
y_predict_train = cv_pipeline.decision_function(x_train)
y_predict_test = cv_pipeline.decision_function(x_test)
metrics_train = get_threshold_metrics(y_train, y_predict_train,
                                      drop_intermediate=keep_inter)
metrics_test = get_threshold_metrics(y_test, y_predict_test,
                                     drop_intermediate=keep_inter)

# Rerun "cross validation" for the best hyperparameter set to define
# cross-validation disease-specific performance. Each sample prediction is
# based on the fold that the sample was in the testing partition
y_cv = cross_val_predict(cv_pipeline.best_estimator_, X=x_train, y=y_train,
                         cv=folds, method='decision_function')
metrics_cv = get_threshold_metrics(y_train, y_cv,
                                   drop_intermediate=keep_inter)

# Decide to save ROC results to file
if keep_inter:
    train_roc = metrics_train['roc_df']
    train_roc = train_roc.assign(train_type='train')
    test_roc = metrics_test['roc_df']
    test_roc = test_roc.assign(train_type='test')
    cv_roc = metrics_cv['roc_df']
    cv_roc = cv_roc.assign(train_type='cv')
    full_roc_df = pd.concat([train_roc, test_roc, cv_roc])
    full_roc_df = full_roc_df.assign(disease='PanCan')

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
plt.savefig(full_roc_file, dpi=600, format='svg', bbox_inches='tight')
plt.close()

# Plot PR
sns.set_style("whitegrid")
plt.figure(figsize=(2.7, 2.4))
total_aupr = {}
colors = ['blue', 'green', 'orange']
idx = 0
for label, metrics in [('Training', metrics_train), ('Testing', metrics_test),
                       ('CV', metrics_cv)]:

    pr_df = metrics['pr_df']
    plt.plot(pr_df.recall, pr_df.precision,
             label='{} (AUPR = {:.1%})'.format(label, metrics['aupr']),
             linewidth=2, c=colors[idx])
    total_aupr[label] = metrics['aupr']
    idx += 1
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall', fontsize=8)
plt.ylabel('Precision', fontsize=8)
plt.title('')
plt.tick_params(labelsize=8)
plt.legend(bbox_to_anchor=(0.2, -0.45, 0.7, .202), loc=0, borderaxespad=0.,
           fontsize=7.5)
plt.tight_layout()
plt.savefig(full_pr_file, dpi=600, format='svg', bbox_inches='tight')
plt.close()

# disease specific performance
disease_metrics = {}
for disease in diseases:
    # Get all samples in current disease
    sample_sub = y_sub[y_sub == disease].index

    # Get true and predicted training labels
    y_disease_train = y_train[y_train.index.isin(sample_sub)]
    y_disease_predict_train = y_predict_train[y_train.index.isin(sample_sub)]

    # Get true and predicted testing labels
    y_disease_test = y_test[y_test.index.isin(sample_sub)]
    y_disease_predict_test = y_predict_test[y_test.index.isin(sample_sub)]

    # Get predicted labels for samples when they were in cross validation set
    # The true labels are y_pred_train
    y_disease_predict_cv = y_cv[y_train.index.isin(sample_sub)]

    # Get classifier performance metrics for three scenarios for each disease
    met_train_dis = get_threshold_metrics(y_disease_train,
                                          y_disease_predict_train,
                                          disease=disease,
                                          drop_intermediate=keep_inter)
    met_test_dis = get_threshold_metrics(y_disease_test,
                                         y_disease_predict_test,
                                         disease=disease,
                                         drop_intermediate=keep_inter)
    met_cv_dis = get_threshold_metrics(y_disease_train,
                                       y_disease_predict_cv,
                                       disease=disease,
                                       drop_intermediate=keep_inter)

    if keep_inter:
        train_roc = met_train_dis['roc_df']
        train_roc = train_roc.assign(train_type='train')
        test_roc = met_test_dis['roc_df']
        test_roc = test_roc.assign(train_type='test')
        cv_roc = met_cv_dis['roc_df']
        cv_roc = cv_roc.assign(train_type='cv')
        full_dis_roc_df = train_roc.append(test_roc).append(cv_roc)
        full_dis_roc_df = full_dis_roc_df.assign(disease=disease)
        full_roc_df = full_roc_df.append(full_dis_roc_df)

    # Store results in disease indexed dictionary
    disease_metrics[disease] = [met_train_dis, met_test_dis, met_cv_dis]

disease_auroc = {}
disease_aupr = {}
for disease, metrics_val in disease_metrics.items():
    met_train, met_test, met_cv = metrics_val
    disease_pr_sub_file = '{}_pred_{}.svg'.format(disease_pr_file, disease)
    disease_roc_sub_file = '{}_pred_{}.svg'.format(disease_roc_file, disease)

    # Plot disease specific PR
    plt.figure(figsize=(2.7, 2.4))
    aupr = []
    idx = 0
    for label, metrics in [('Training', met_train), ('Testing', met_test),
                           ('CV', met_cv)]:
        pr_df = metrics['pr_df']
        plt.plot(pr_df.recall, pr_df.precision,
                 label='{} (AUPR = {:.1%})'.format(label, metrics['aupr']),
                 linewidth=2, c=colors[idx])
        aupr.append(metrics['aupr'])
        idx += 1
    disease_aupr[disease] = aupr

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall', fontsize=8)
    plt.ylabel('Precision', fontsize=8)
    plt.title('')
    plt.tick_params(labelsize=8)
    plt.legend(bbox_to_anchor=(0.2, -0.45, 0.7, .202), loc=0, borderaxespad=0.,
               fontsize=7.5)
    plt.tight_layout()
    plt.savefig(disease_pr_sub_file, dpi=600, format='svg',
                bbox_inches='tight')
    plt.close()

    # Plot disease specific ROC
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
    disease_auroc[disease] = auroc

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=8)
    plt.ylabel('True Positive Rate', fontsize=8)
    plt.title('')
    plt.tick_params(labelsize=8)
    plt.legend(bbox_to_anchor=(0.2, -0.45, 0.7, .202), loc=0, borderaxespad=0.,
               fontsize=7.5)
    plt.tight_layout()
    plt.savefig(disease_roc_sub_file, dpi=600, format='svg',
                bbox_inches='tight')
    plt.close()

disease_auroc_df = pd.DataFrame(disease_auroc, index=['Train', 'Test',
                                                      'Cross Validation']).T
disease_auroc_df = disease_auroc_df.sort_values('Cross Validation',
                                                ascending=False)
ax = disease_auroc_df.plot(kind='bar', title='Disease Specific Performance')
ax.set_ylabel('AUROC')
plt.tight_layout()
plt.savefig(dis_summary_auroc_file, dpi=600, format='svg', bbox_inches='tight')
plt.close()

disease_aupr_df = pd.DataFrame(disease_aupr, index=['Train', 'Test',
                                                    'Cross Validation']).T
disease_aupr_df = disease_aupr_df.sort_values('Cross Validation',
                                              ascending=False)
ax = disease_aupr_df.plot(kind='bar', title='Disease Specific Performance')
ax.set_ylabel('AUPR')
plt.tight_layout()
plt.savefig(dis_summary_aupr_file, dpi=600, format='svg', bbox_inches='tight')
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

if keep_inter:
    full_roc_df.to_csv(roc_results_file, sep='\t')

# Apply the same classifier previously built to predict alternative genes
if alt_genes[0] is not 'None':
    # Classifying alternative mutations
    y_alt = mutation_df[alt_genes]

    # Add copy number info if applicable
    if copy_number:
        y_alt = integrate_copy_number(y=y_alt, cancer_genes_df=cancer_genes,
                                      genes=alt_genes, loss_df=copy_loss_df,
                                      gain_df=copy_gain_df)
    # Append disease id
    y_alt = y_alt.assign(total_status=y_alt.max(axis=1))
    y_alt = y_alt.reset_index().merge(sample_freeze,
                                      how='left').set_index('SAMPLE_BARCODE')

    # Filter data
    alt_count_df = y_alt.groupby('DISEASE').sum()
    alt_prop_df = alt_count_df.divide(y_alt['DISEASE'].value_counts(sort=False)
                                                      .sort_index(), axis=0)

    alt_count_table = alt_count_df.merge(alt_prop_df, left_index=True,
                                         right_index=True)
    alt_count_table.to_csv(alt_count_table_file)

    mut_co = alt_count_df['total_status']
    prop = alt_prop_df['total_status']

    if alt_diseases[0] == 'Auto':
        alt_filter_dis = (mut_co > alt_filter_count) & (prop > alt_filter_prop)
        alt_diseases = alt_filter_dis.index[alt_filter_dis].tolist()

    # Subset data
    y_alt_df = y_alt[y_alt.DISEASE.isin(alt_diseases)].total_status
    common_alt_samples = list(set(y_alt_df.index) & set(rnaseq_full_df.index))

    y_alt_df = y_alt_df.loc[common_alt_samples]
    rnaseq_alt_df = rnaseq_full_df.loc[y_alt_df.index, :]

    y_alt_matrix = mut_burden.merge(pd.DataFrame(y_alt_df), right_index=True,
                                    left_on='SAMPLE_BARCODE')\
                             .set_index('SAMPLE_BARCODE')

    # Add Covariate Info to alternative y matrix
    y_alt_sub = y_alt.loc[y_alt_matrix.index]['DISEASE']
    covar_dummy_alt = pd.get_dummies(sample_freeze['DISEASE']).astype(int)
    covar_dummy_alt.index = sample_freeze['SAMPLE_BARCODE']
    covar_alt = covar_dummy_alt.merge(y_alt_matrix, right_index=True,
                                      left_index=True)
    covar_alt = covar_alt.drop('total_status', axis=1)
    y_alt_df = y_alt_df.loc[y_alt_sub.index]

    # Process alternative x matrix
    x_alt_df = rnaseq_alt_df.loc[y_alt_df.index, :]
    if x_matrix == 'raw':
        x_alt_df = x_alt_df.loc[:, mad_genes]

    x_alt_df_update = pd.DataFrame(fitted_scaler.transform(x_alt_df),
                                   columns=x_alt_df.columns)
    x_alt_df_update.index = x_alt_df.index
    x_alt_df = x_alt_df_update.merge(covar_alt, left_index=True,
                                     right_index=True)

    # Apply the previously fit model to predict the alternate Y matrix
    y_alt_cv = cv_pipeline.decision_function(X=x_alt_df)
    alt_metrics_cv = get_threshold_metrics(y_alt_df, y_alt_cv,
                                           drop_intermediate=keep_inter)

    validation_metrics = {}
    val_x_type = {}
    for disease in alt_diseases:
        sample_dis = y_alt_sub[y_alt_sub == disease].index

        # Subset full data if it has not been trained on
        if disease not in diseases:
            x_sub = x_alt_df.loc[sample_dis]
            y_sub = y_alt_df[sample_dis]
            category = 'Full'

        # Only subset to the holdout set if data was trained on
        else:
            x_sub = x_test.loc[x_test.index.isin(sample_dis)]
            y_sub = y_test[y_test.index.isin(sample_dis)]
            category = 'Holdout'

        neg, pos = y_sub.value_counts()
        val_x_type[disease] = [category, neg, pos]
        y_pred_alt = cv_pipeline.decision_function(x_sub)
        y_pred_alt_cv = y_alt_cv[y_alt_df.index.isin(y_sub.index)]

        alt_metrics_dis = get_threshold_metrics(y_sub, y_pred_alt,
                                                disease=disease,
                                                drop_intermediate=keep_inter)
        alt_metrics_di_cv = get_threshold_metrics(y_sub, y_pred_alt_cv,
                                                  disease=disease,
                                                  drop_intermediate=keep_inter)
        validation_metrics[disease] = [alt_metrics_dis, alt_metrics_di_cv]

    # Compile a summary dataframe
    val_x_type = pd.DataFrame.from_dict(val_x_type)
    val_x_type.index = ['class', 'negatives', 'positives']
    val_x_type.to_csv(alt_gene_summary_file, sep='\t')

    alt_disease_auroc = {}
    alt_disease_aupr = {}
    for disease, metrics_val in validation_metrics.items():
        met_test, met_cv = metrics_val
        alt_disease_auroc[disease] = [met_test['auroc'], met_cv['auroc']]
        alt_disease_aupr[disease] = [met_test['aupr'], met_cv['aupr']]

    # Plot alternative gene cancer-type specific AUROC plots
    alt_disease_auroc_df = pd.DataFrame(alt_disease_auroc,
                                        index=['Hold Out', 'Full Data']).T
    alt_disease_auroc_df = alt_disease_auroc_df.sort_values('Full Data',
                                                            ascending=False)
    ax = alt_disease_auroc_df.plot(kind='bar', title='Alt Gene Performance')
    ax.set_ylim([0, 1])
    ax.set_ylabel('AUROC')
    plt.tight_layout()
    plt.savefig(alt_gene_auroc_file, dpi=600, format='svg',
                bbox_inches='tight')
    plt.close()

    # Plot alternative gene cancer-type specific AUPR plots
    alt_disease_aupr_df = pd.DataFrame(alt_disease_aupr,
                                       index=['Hold Out', 'Full Data']).T
    alt_disease_aupr_df = alt_disease_aupr_df.sort_values('Full Data',
                                                          ascending=False)
    ax = alt_disease_aupr_df.plot(kind='bar', title='Alt Gene Performance')
    ax.set_ylim([0, 1])
    ax.set_ylabel('AUPR')
    plt.tight_layout()
    plt.savefig(alt_gene_aupr_file, dpi=600, format='svg',
                bbox_inches='tight')
    plt.close()

# Write a summary for the inputs and outputs of the classifier
with open(os.path.join(base_folder, 'classifier_summary.txt'), 'w') as sum_fh:
    summarywriter = csv.writer(sum_fh, delimiter='\t')

    # Summarize parameters
    summarywriter.writerow(['Parameters:'])
    summarywriter.writerow(['Genes:'] + genes)
    summarywriter.writerow(['Diseases:'] + diseases)
    summarywriter.writerow(['Alternative Genes:'] + alt_genes)
    summarywriter.writerow(['Alternative Diseases:'] + alt_diseases)
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
    summarywriter.writerow(['Training AUPR:', metrics_train['aupr']])
    summarywriter.writerow(['Testing AUPR:', metrics_test['aupr']])
    summarywriter.writerow(['Cross Validation AUPR:', metrics_cv['aupr']])
    summarywriter.writerow(['Disease specific performance:'])
    for disease, auroc in disease_auroc.items():
        summarywriter.writerow(['', disease, 'Training AUROC:', auroc[0],
                                'Testing AUROC:', auroc[1],
                                'Cross Validation AUROC:', auroc[2]])
    for disease, aupr in disease_aupr.items():
        summarywriter.writerow(['', disease, 'Training AUPR:', aupr[0],
                                'Testing AUPR:', aupr[1],
                                'Cross Validation AUPR:', aupr[2]])
    if alt_genes[0] is not 'None':
        summarywriter.writerow(['Alternate gene performance:'] + alt_genes)
        summarywriter.writerow(['Alternative gene AUROC:',
                                str(alt_metrics_cv['auroc'])])
        summarywriter.writerow(['Alternative gene AUPR:',
                                str(alt_metrics_cv['aupr'])])
        for alt_dis, alt_auroc in alt_disease_auroc.items():
            summarywriter.writerow(['', alt_dis,
                                    'Holdout AUROC:', alt_auroc[0],
                                    'Full Data AUROC:', alt_auroc[1],
                                    'Category:', val_x_type[alt_dis]['class'],
                                    'num_positive:',
                                    str(val_x_type[alt_dis]['positives']),
                                    'num_negatives:',
                                    str(val_x_type[alt_dis]['negatives'])])
        for alt_dis, alt_aupr in alt_disease_aupr.items():
            summarywriter.writerow(['', alt_dis,
                                    'Holdout AUPR:', alt_aupr[0],
                                    'Full Data AUPR:', alt_aupr[1],
                                    'Category:', val_x_type[alt_dis]['class'],
                                    'num_positive:',
                                    str(val_x_type[alt_dis]['positives']),
                                    'num_negatives:',
                                    str(val_x_type[alt_dis]['negatives'])])
