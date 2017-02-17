"""
Gregory Way 2017
PanCancer Classifier
tcga_util.py

Usage: For import only
"""


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

    return y
