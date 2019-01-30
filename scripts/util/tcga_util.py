"""
Gregory Way 2017
PanCancer Classifier
tcga_util.py

Usage: For import only
"""


def get_args():
    """
    Get arguments for the main pancancer classifier script
    """
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genes',
                        help='Comma separated string of HUGO gene symbols')
    parser.add_argument('-t', '--diseases', default='Auto',
                        help='Comma sep string of TCGA disease acronyms. '
                             'If no arguments are passed, filtering will '
                             'default to options given in --filter_count and '
                             '--filter_prop.')
    parser.add_argument('-f', '--folds', default='5', type=int,
                        help='Number of cross validation folds to perform')
    parser.add_argument('-d', '--drop', action='store_true',
                        help='Decision to drop input genes from X matrix')
    parser.add_argument('-u', '--copy_number', action='store_true',
                        help='Supplement Y matrix with copy number events')
    parser.add_argument('-c', '--filter_count', default=15, type=int,
                        help='Min number of mutations in diseases to include')
    parser.add_argument('-p', '--filter_prop', default=0.05, type=float,
                        help='Min proportion of positives to include disease')
    parser.add_argument('-n', '--num_features', default=8000, type=int,
                        help='Number of MAD genes to include in classifier')
    parser.add_argument('-a', '--alphas', default='0.1,0.15,0.2,0.5,0.8,1',
                        help='the alphas for parameter sweep')
    parser.add_argument('-l', '--l1_ratios', default='0,0.1,0.15,0.18,0.2,0.3',
                        help='the l1 ratios for parameter sweep')
    parser.add_argument('-b', '--alt_genes', default='None',
                        help='alternative genes to test performance')
    parser.add_argument('-s', '--alt_diseases', default="Auto",
                        help='The alternative diseases to test performance')
    parser.add_argument('-i', '--alt_filter_count', default=15, type=int,
                        help='Min number of mutations in disease to include')
    parser.add_argument('-r', '--alt_filter_prop', default=0.05, type=float,
                        help='Min proportion of positives to include disease')
    parser.add_argument('-o', '--alt_folder', default='Auto',
                        help='Provide an alternative folder to save results')
    parser.add_argument('-v', '--remove_hyper', action='store_true',
                        help='Remove hypermutated samples')
    parser.add_argument('-k', '--keep_intermediate', action='store_true',
                        help='Keep intermediate ROC values for plotting')
    parser.add_argument('-x', '--x_matrix', default='raw',
                        help='Filename of features to use in model')
    parser.add_argument('-e', '--shuffled', action='store_true',
                        help='Shuffle the input gene exprs matrix alongside')
    parser.add_argument('--shuffled_before_training', action='store_true',
                        help='Shuffle the gene exprs matrix before training')
    parser.add_argument('-m', '--no_mutation', action='store_false',
                        help='Remove mutation data from y matrix')
    parser.add_argument('-z', '--drop_rasopathy', action='store_true',
                        help='Decision to drop rasopathy genes from X matrix')
    parser.add_argument('-q', '--drop_expression', action='store_true',
                        help='Decision to drop gene expression values from X')
    parser.add_argument('-j', '--drop_covariates', action='store_true',
                        help='Decision to drop covariate information from X')

    args = parser.parse_args()
    return args


def get_threshold_metrics(y_true, y_pred, drop_intermediate=False,
                          disease='all'):
    """
    Retrieve true/false positive rates and auroc/aupr for class predictions

    Arguments:
    y_true - an array of gold standard mutation status
    y_pred - an array of predicted mutation status
    disease - a string that includes the corresponding TCGA study acronym

    Output:
    dict of AUROC, AUPR, pandas dataframes of ROC and PR data, and cancer-type
    """
    import pandas as pd
    from sklearn.metrics import roc_auc_score, roc_curve
    from sklearn.metrics import precision_recall_curve, average_precision_score

    roc_columns = ['fpr', 'tpr', 'threshold']
    pr_columns = ['precision', 'recall', 'threshold']

    if drop_intermediate:
        roc_items = zip(roc_columns,
                        roc_curve(y_true, y_pred, drop_intermediate=False))
    else:
        roc_items = zip(roc_columns, roc_curve(y_true, y_pred))

    roc_df = pd.DataFrame.from_dict(dict(roc_items))

    prec, rec, thresh = precision_recall_curve(y_true, y_pred)
    pr_df = pd.DataFrame.from_records([prec, rec]).T
    pr_df = pd.concat([pr_df, pd.Series(thresh)], ignore_index=True, axis=1)
    pr_df.columns = pr_columns

    auroc = roc_auc_score(y_true, y_pred, average='weighted')
    aupr = average_precision_score(y_true, y_pred, average='weighted')

    return {'auroc': auroc, 'aupr': aupr, 'roc_df': roc_df,
            'pr_df': pr_df, 'disease': disease}


def integrate_copy_number(y, cancer_genes_df, genes, loss_df, gain_df,
                          include_mutation=True):
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
    include_mutation - boolean to decide to include mutation status
    """

    # Find if the input genes are in this master list
    genes_sub = cancer_genes_df[cancer_genes_df['Gene Symbol'].isin(genes)]

    # Add status to the Y matrix depending on if the gene is a tumor suppressor
    # or an oncogene. An oncogene can be activated with copy number gains, but
    # a tumor suppressor is inactivated with copy number loss
    tumor_suppressor = genes_sub[genes_sub['Classification*'] == 'TSG']
    oncogene = genes_sub[genes_sub['Classification*'] == 'Oncogene']

    copy_loss_sub = loss_df[tumor_suppressor['Gene Symbol']]
    copy_gain_sub = gain_df[oncogene['Gene Symbol']]

    # Append to column names for visualization
    copy_loss_sub.columns = [col + '_loss' for col in copy_loss_sub.columns]
    copy_gain_sub.columns = [col + '_gain' for col in copy_gain_sub.columns]

    # Add columns to y matrix
    y = y.join(copy_loss_sub)
    y = y.join(copy_gain_sub)

    # Fill missing data with zero (measured mutation but not copy number)
    y = y.fillna(0)
    y = y.astype(int)

    if not include_mutation:
        y = y.drop(genes, axis=1)
    return y


def shuffle_columns(gene):
    """
    To be used in an `apply` pandas func to shuffle columns around a datafame
    Import only
    """
    import numpy as np
    return np.random.permutation(gene.tolist())
