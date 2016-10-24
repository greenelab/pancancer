"""
Gregory Way 2016
Modified from https://github.com/cognoma/machine-learning/
PanCancer NF1 Classifier
pancancer_classifier.py

Usage: Run in command line with required command argument:

        python pancancer_classifier.py --genes $GENES

Where GENES is a comma separated string. There are also optional arguments:

    --tissues           comma separated string of tissues to build classifier
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
    --xena              if added, uses publicly available data instead of
                        controlled access synapse data

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
from sklearn.feature_selection import SelectKBest
from statsmodels.robust.scale import mad

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genes',
                    help='Comma separated string of HUGO gene symbols')
parser.add_argument('-t', '--tissues', default='Auto',
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
parser.add_argument('-a', '--alphas', default='0.01,0.1,0.15,0.2,0.5,0.8',
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
parser.add_argument('-x', '--xena', action='store_true',
                    help='Use Xena data instead of controlled-access data')
args = parser.parse_args()

# Load command arguments
genes = args.genes.split(',')
tissues = args.tissues.split(',')
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
xena = args.xena

warnings.filterwarnings('ignore',
                        message='Changing the shape of non-C contiguous array')


def integrage_copy_number(y, cancer_genes_df, genes, loss_df, gain_df):
    """
    Function to integrate copy number data to define gene activation or gene
    inactivation events. Copy number loss results in gene inactivation events
    and is important for tumor suppressor genes while copy number gain results
    in gene activation events and is important for oncogenes.

    Arguments:
    y - pandas dataframe samples by genes where a 1 indicates event
    cancer_genes_df - a dataframe listing bona fide cancer genes as defined by
                      the 20/20 rule in Vogelstein et al. 2013
    genes - the input list of genes to build the classifier for
    loss_df - a sample by gene dataframe listing copy number loss events
    gain_df - a sample by gene dataframe listing copy number gain events
    """

    # Find if the input genes are in this master list
    cancer_genes_sub = cancer_genes[cancer_genes['Gene Symbol'].isin(genes)]

    # Add status to the Y matrix depending on if the gene is a tumor suppressor
    # or an oncogene. An oncogene can be activated with copy number gains, but
    # a tumor suppressor is inactivated with copy number loss
    tumor_suppressor = cancer_genes_sub[cancer_genes_sub['Classification*'] ==
                                        'TSG']['Gene Symbol']
    oncogene = cancer_genes_sub[cancer_genes_sub['Classification*'] ==
                                'Oncogene']['Gene Symbol']
    copy_loss_sub = loss_df[tumor_suppressor]
    copy_gain_sub = gain_df[oncogene]

    # Append to column names for visualization
    copy_loss_sub.columns = [col + '_loss' for col in copy_loss_sub.columns]
    copy_gain_sub.columns = [col + '_gain' for col in copy_gain_sub.columns]

    # Add columns to y matrix
    y = y.join(copy_loss_sub)
    y = y.join(copy_gain_sub)

    # Fill missing data with zero (measured mutation but not copy number)
    y = y.fillna(0)
    y = y.astype(int)

    return y


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
if xena:
    base_add = 'xena'
else:
    base_add = 'synapse'

base_folder = '{}__tissues_{}__genes_{}'.format(base_add,
                                                args.tissues.replace(',', '_'),
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

mutation_heatmap_file = '{}/mutation_heatmap.png'.format(base_folder)
unfiltered_heatmap_file = '{}/unfiltered_heatmap.png'.format(base_folder)
pos_count_file = '{}/total_mutation_count.tsv'.format(base_folder)
cv_heatmap_file = '{}/cv_heatmap.png'.format(base_folder)
full_roc_file = '{}/all_tissue_roc.png'.format(base_folder)
tissue_roc_file = '{}/tissue/classifier_'.format(base_folder)
tissue_summary_file = '{}/tissue_summary.png'.format(base_folder)
classifier_file = '{}/classifier_coefficients.tsv'.format(base_folder)
alt_gene_tissue_file = '{}/alt_gene_{}_alt_tissues_{}_classifier.png'\
                       .format(base_folder, args.alt_genes.replace(',', '_'),
                               args.alt_tissues.replace(',', '_'))

# Load Datasets
if xena:
    expr_fh = 'data/xena/expression.tsv.bz2'
    mut_fh = 'data/xena/mutation-matrix.bz2'
    clin_fh = 'data/xena/samples.tsv'
    gene_map = pd.read_table('data/xena/HiSeqV2-gene-map.tsv')
    tcga_acronym = pd.read_table('data/tcga_dictionary.tsv')
    ac = dict(zip(tcga_acronym.tissue, tcga_acronym.acronym))
else:
    expr_fh = 'data/rnaseq_data.tsv'
    mut_fh = 'data/mutation_table.tsv'
    clin_fh = 'data/clinical_data.tsv'

rnaseq_df = pd.read_table(expr_fh, index_col=0)
mutation_df = pd.read_table(mut_fh, index_col=0)
clinical_df = pd.read_table(clin_fh, index_col=0, low_memory=False)

# Subset data
if xena:
    gene_map_dict = {str(k): str(v) for v, k in zip(gene_map.symbol,
                                                    gene_map.entrez_gene_id)}
    rnaseq_df = rnaseq_df.rename(index=str, columns=gene_map_dict)
    mutation_df = mutation_df.rename(index=str, columns=gene_map_dict)
    clinical_df = clinical_df.assign(acronym=clinical_df
                                     .disease.replace(to_replace=ac))
else:
    rnaseq_df.drop('SLC35E2', axis=0, inplace=True)
    rnaseq_df = rnaseq_df.T
    base_add = 'synapse'

# Construct data for classifier
y = mutation_df[genes]

if drop:
    for gene in genes:
        rnaseq_df.drop(gene, axis=1, inplace=True)

# Incorporate copy number for gene activation/inactivation
if copy_number:
    # Load copy number matrices
    copy_loss_df = pd.read_table('data/copy_number_loss_status.tsv',
                                 index_col=0)
    copy_gain_df = pd.read_table('data/copy_number_gain_status.tsv',
                                 index_col=0)

    # Load cancer gene classification table
    cancer_genes = pd.read_table('data/vogelstein_cancergenes.tsv')

    y = integrage_copy_number(y=y, cancer_genes_df=cancer_genes, genes=genes,
                              loss_df=copy_loss_df, gain_df=copy_gain_df)

# Process y matrix
y = y.assign(disease=clinical_df['acronym'])
y = y.assign(total_status=y.max(axis=1))
num_positives = y.groupby('disease').apply(lambda x: x['total_status'].sum())
group_count_df = y.groupby('disease').sum()
heatmap_df = group_count_df.divide(y['disease'].value_counts(sort=False)
                                               .sort_index(), axis=0)

# Plot unfiltered heatmap
ax = sns.heatmap(heatmap_df, annot=True)
ax.set(xlabel='gene', ylabel='study')
plt.yticks(rotation=0)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(unfiltered_heatmap_file)
plt.close()

# Filter tissues
mut_count = group_count_df['total_status']
prop = heatmap_df['total_status']

if tissues[0] == 'Auto':
    filter_tissue = (mut_count > filter_count) & (prop > filter_prop)
    tissues = filter_tissue.index[filter_tissue].tolist()

y_df = y[y.disease.isin(tissues)].total_status
x_df = rnaseq_df.ix[y_df.index, :]
clinical_sub = clinical_df.ix[y_df.index]

# Summary heatmap of mutation counts after filtering
filter_heatmap_df = group_count_df[group_count_df.index.isin(tissues)]

# Show mutation frequency heatmap
filter_heatmap_df = filter_heatmap_df.divide(y['disease']
                                             .value_counts(sort=False)
                                             .sort_index(), axis=0).dropna()

ax = sns.heatmap(filter_heatmap_df, annot=True)
ax.set(xlabel='gene', ylabel='study')
plt.yticks(rotation=0)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(mutation_heatmap_file)
plt.close()

# Save total count matrix
pos_counts = group_count_df.join(y['disease'].value_counts())
pos_counts.to_csv(pos_count_file, sep='\t')

# Build classifier pipeline
strat = clinical_df.ix[y_df.index]['acronym'].str.cat(y_df.astype(str))
x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.1,
                                                    random_state=0,
                                                    stratify=strat)
feature_select = SelectKBest(fs_mad, k=num_features_kept)

clf_parameters = {'select__k': [8000], 'classify__loss': ['log'],
                  'classify__penalty': ['elasticnet'],
                  'classify__alpha': alphas, 'classify__l1_ratio': l1_ratios}

estimator = Pipeline(steps=[('select', SelectKBest(fs_mad)),
                            ('standardize', StandardScaler()),
                            ('classify', SGDClassifier(random_state=0,
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
ax = sns.heatmap(cv_score_mat, annot=True, fmt='.2%')
ax.set_xlabel('Regularization strength multiplier (alpha)')
ax.set_ylabel('Elastic net mixing parameter (l1_ratio)')
plt.tight_layout()
plt.savefig(cv_heatmap_file)
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
plt.figure()
total_auroc = {}
for label, metrics in [('Training', metrics_train), ('Testing', metrics_test),
                       ('Cross Validation', metrics_cv)]:
    roc_df = metrics['roc_df']
    plt.plot(roc_df.fpr, roc_df.tpr,
             label='{} (AUROC = {:.1%})'.format(label, metrics['auroc']))
    total_auroc[label] = metrics['auroc']

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(full_roc_file)
plt.close()

# Plot tissue specific performance
tissue_metrics = {}
for tissue in clinical_sub.acronym.unique():
    sample_sub = clinical_sub[clinical_sub['acronym'] == tissue].index.values

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
    tissue_roc_sub_file = tissue_roc_file + '_pred_' + tissue + '.png'

    plt.figure()
    auroc = []
    for label, metrics in [('Training', met_train), ('Testing', met_test),
                           ('Cross Validation', met_cv)]:
        roc_df = metrics['roc_df']
        plt.plot(roc_df.fpr, roc_df.tpr,
                 label='{} (AUROC = {:.1%})'.format(label, metrics['auroc']))
        auroc.append(metrics['auroc'])

    tissue_auroc[tissue] = auroc
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Predicting Event in {}'.format(tissue))
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(tissue_roc_sub_file)
    plt.close()

tissue_results = pd.DataFrame(tissue_auroc, index=['Train', 'Test',
                                                   'Cross Validation']).T
tissue_results = tissue_results.sort_values('Cross Validation',
                                            ascending=False)
ax = tissue_results.plot(kind='bar', title='Disease Specific Performance')
ax.set_ylabel('AUROC')
plt.tight_layout()
plt.savefig(tissue_summary_file)
plt.close()

# Save classifier coefficients
final_pipeline = cv_pipeline.best_estimator_
final_classifier = final_pipeline.named_steps['classify']

select_indices = final_pipeline.named_steps['select'].transform(
    np.arange(len(rnaseq_df.columns)).reshape(1, -1)).tolist()

coef_df = pd.DataFrame.from_items([
    ('feature', rnaseq_df.columns[select_indices]),
    ('weight', final_classifier.coef_[0])])

coef_df['abs'] = coef_df['weight'].abs()
coef_df = coef_df.sort_values('abs', ascending=False)
coef_df.to_csv(classifier_file, sep='\t')

# Apply classifier to predict alternative genes
if alt_genes[0] is not 'None':
    # Classifying alt_gene mutations
    y_alt = mutation_df[alt_genes]

    # Add copy number info if applicable
    if copy_number:
        y_alt = integrage_copy_number(y=y_alt, cancer_genes_df=cancer_genes,
                                      genes=alt_genes, loss_df=copy_loss_df,
                                      gain_df=copy_gain_df)
    # Append tissue id
    y_alt = y_alt.assign(disease=clinical_df['acronym'])
    y_alt = y_alt.assign(total_status=y_alt.max(axis=1))

    # Filter data
    group_count_df = y_alt.groupby('disease').sum()
    group_prop_df = group_count_df.divide(y_alt['disease']
                                          .value_counts(sort=False)
                                          .sort_index(), axis=0)

    mut_co = group_count_df['total_status']
    prop = group_prop_df['total_status']

    if alt_tissues[0] == 'Auto':
        alt_filter_tis = (mut_co > alt_filter_count) & (prop > alt_filter_prop)
        alt_tissues = alt_filter_tis.index[alt_filter_tis].tolist()

    # Subset data
    y_alt_df = y_alt[y_alt.disease.isin(alt_tissues)].total_status
    clin_sub = clinical_df.ix[y_alt_df.index]

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
        sample_tissues = clin_sub[clin_sub['acronym'] == tissue].index.values

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
    summarywriter.writerow(['Number of Features:', str(num_features_kept)])
    summarywriter.writerow(['Drop Gene:', drop])
    summarywriter.writerow(['Copy Number:', copy_number])
    summarywriter.writerow(['Alphas:'] + alphas)
    summarywriter.writerow(['L1_ratios:'] + l1_ratios)
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
        summarywriter.writerow(['Alternate gene performance:'])
        summarywriter.writerow(['Alternative gene AUROC:',
                                str(alt_metrics_cv['auroc'])])
        for alt_tissue, alt_auroc in alt_tissue_auroc.items():
            summarywriter.writerow(['', alt_tissue, 'Testing AUROC:',
                                    alt_auroc[0], 'Cross validation AUROC',
                                    alt_auroc[1], 'Data:',
                                    validation_x_type[alt_tissue]])
