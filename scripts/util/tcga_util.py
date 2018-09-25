"""
Gregory Way 2017
PanCancer Classifier
tcga_util.py

Usage: For import only
"""


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
