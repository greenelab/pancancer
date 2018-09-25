
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import roc_auc_score, average_precision_score


# In[2]:


# Get the current working directory
cwd = os.getcwd()

# Ensure that the path is starting in the scripts directory
if not cwd.split('/')[-1] == 'scripts':
    os.chdir(os.path.join(cwd, 'scripts'))


# In[3]:


def get_gene_auroc(x, w):
    score = roc_auc_score(x, w, average='weighted')
    return(score)

def get_gene_auprc(x, w):
    score = average_precision_score(x, w, average='weighted')
    return(score)


# In[4]:


ras_folder = os.path.join('..', 'classifiers', 'RAS')


# In[5]:


# Load Datasets
mut_file = os.path.join('..', 'data', 'pancan_mutation_freeze.tsv.gz')
sample_freeze_file = os.path.join('..', 'data', 'sampleset_freeze.csv')
copy_loss_file = os.path.join('..', 'data', 'copy_number_loss_status.tsv.gz')
copy_gain_file = os.path.join('..', 'data', 'copy_number_gain_status.tsv.gz')

mutation_df = pd.read_table(mut_file, index_col=0)
sample_freeze = pd.read_table(sample_freeze_file, index_col=0)
copy_loss_df = pd.read_table(copy_loss_file, index_col=0)
copy_gain_df = pd.read_table(copy_gain_file, index_col=0)


# In[6]:


# Load Ras Pathway Genes
ras_genes_file = os.path.join('..', 'data', 'ras_genes.csv')
ras_genes_df = pd.read_table(ras_genes_file)


# In[7]:


# Load classifier weights
ras_decision_file = os.path.join(ras_folder, 'classifier_decisions.tsv')
ras_decisions_df = pd.read_table(ras_decision_file)
ras_decisions_df.head()


# In[8]:


ras_mutations_df = mutation_df[ras_genes_df['genes']]

# Add status to the Y matrix depending on if the gene is a tumor suppressor
# or an oncogene. An oncogene can be activated with copy number gains, but
# a tumor suppressor is inactivated with copy number loss
oncogene = ras_genes_df[ras_genes_df['og_tsg'] == 'OG']
tumor_suppressor = ras_genes_df[ras_genes_df['og_tsg'] == 'TSG']

# Subset copy number information
ras_copy_gain_sub_df = copy_gain_df[oncogene['genes']]
ras_copy_loss_sub_df = copy_loss_df[tumor_suppressor['genes']]

# Combine Copy Number data
ras_copy_df = pd.concat([ras_copy_gain_sub_df, ras_copy_loss_sub_df], axis=1)


# In[9]:


ras_status_df = ras_mutations_df + ras_copy_df
ras_status_df[ras_status_df == 2] = 1


# In[10]:


subset_columns = ['SAMPLE_BARCODE', 'DISEASE', 'weight', 'total_status', 'log10_mut',
                  'hypermutated', 'include']
ras_decisions_subset_df = ras_decisions_df[subset_columns]
ras_full_status_df = ras_status_df.merge(ras_decisions_subset_df, left_index=True,
                                         right_on='SAMPLE_BARCODE')
ras_full_status_df.index = ras_full_status_df['SAMPLE_BARCODE']


# In[11]:


# Remove hyper mutated samples
burden_filter = ras_full_status_df['hypermutated'] == 0
burden_filter = burden_filter & ras_full_status_df['log10_mut'] < 5 * ras_full_status_df['log10_mut'].std()
ras_full_status_df = ras_full_status_df[burden_filter]
ras_full_status_df.head(3)


# In[12]:


full_auroc = (
    ras_full_status_df[ras_genes_df['genes']]
    .apply(lambda x: get_gene_auroc(x, ras_full_status_df['weight']))
    )

full_auprc = (
    ras_full_status_df[ras_genes_df['genes']]
    .apply(lambda x: get_gene_auprc(x, ras_full_status_df['weight']))
    )


# In[13]:


# Remove Ras positive samples, and recalculate metrics
remove_ras_status = ras_full_status_df[ras_full_status_df['total_status'] == 0]
remove_ras_status_df = remove_ras_status[ras_genes_df['genes']]
remove_ras_status_df = remove_ras_status_df.drop(['KRAS', 'HRAS', 'NRAS'], axis=1)
full_auroc_remove = remove_ras_status_df.apply(lambda x: get_gene_auroc(x, w=remove_ras_status['weight']))
full_auprc_remove = remove_ras_status_df.apply(lambda x: get_gene_auprc(x, w=remove_ras_status['weight']))


# In[14]:


# Get output metrics for Ras classification
output_ras_metrics = pd.concat([full_auroc, full_auroc_remove], axis=1)
output_ras_metrics = output_ras_metrics * 100  # To get percent
output_ras_metrics = output_ras_metrics - 50  # Subtract 50 from AUROC only

# Combine with AUPRC
output_ras_metrics = pd.concat([output_ras_metrics, full_auprc * 100,
                                full_auprc_remove * 100], axis=1)
output_ras_metrics.columns = ['ras_auroc', 'no_ras_auroc', 'ras_auprc', 'no_ras_auprc']

# Fill removed Ras metrics with included metrics
output_ras_metrics['no_ras_auroc'] = (
    output_ras_metrics['no_ras_auroc'].fillna(output_ras_metrics['ras_auroc'])
    )
output_ras_metrics['no_ras_auprc'] = (
    output_ras_metrics['no_ras_auprc'].fillna(output_ras_metrics['ras_auprc'])
    )

# Write results to file
tables_folder = os.path.join(ras_folder, 'tables')
if not os.path.exists(tables_folder):
    os.makedirs(tables_folder)

ras_metric_file = os.path.join(ras_folder, 'tables', 'ras_metrics_pathwaymapper.txt')
output_ras_metrics.to_csv(ras_metric_file, sep='\t')

output_ras_metrics.head()


# In[15]:


# Display Ras pathway metrics
all_samples_ras_pathway_status = ras_full_status_df[ras_genes_df['genes']].max(axis=1)
print('Ras Pathway Performance Summary: All Ras Genes')
print('AUROC:')
print(roc_auc_score(all_samples_ras_pathway_status,
                    ras_full_status_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(all_samples_ras_pathway_status,
                              ras_full_status_df['weight'], average='weighted'))


# In[16]:


print('Ras Pathway Performance Summary: KRAS, NRAS, HRAS')
print('AUROC:')
print(roc_auc_score(ras_full_status_df['total_status'],
                    ras_full_status_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(ras_full_status_df['total_status'],
                              ras_full_status_df['weight'], average='weighted'))


# In[17]:


print('Ras Pathway Performance Summary: Held Out Samples')
held_out_ras_df = ras_full_status_df[ras_full_status_df['include'] == 0]
print('AUROC:')
print(roc_auc_score(held_out_ras_df['total_status'],
                    held_out_ras_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(held_out_ras_df['total_status'],
                              held_out_ras_df['weight'], average='weighted'))


# # Visualize Distribution of AUROC and AUPRC for all genes

# In[18]:


# Subset mutation file by samples
sub_full_mutation_df = mutation_df[burden_filter]
low_mutation_count_filter = (
    sub_full_mutation_df.sum()
    [sub_full_mutation_df.sum() >= 10].sort_values(ascending=False).index
    )
sub_full_mutation_df = sub_full_mutation_df[low_mutation_count_filter]
sub_full_mutation_df.head()


# In[19]:


# Get Metrics for All Genes
all_auprc = sub_full_mutation_df.apply(lambda x: get_gene_auprc(x, w = ras_full_status_df['weight']))
all_auroc = sub_full_mutation_df.apply(lambda x: get_gene_auroc(x, w = ras_full_status_df['weight']))


# In[20]:


# Process file and save results
all_gene_metrics_file = os.path.join(ras_folder, 'tables', 'all_gene_metric_ranks.tsv')

all_genes_auprc_df = pd.DataFrame(all_auprc.sort_values(ascending=False), columns=['auprc'])
all_genes_auroc_df = pd.DataFrame(all_auroc.sort_values(ascending=False), columns=['auroc'])

all_genes_auprc_df = all_genes_auprc_df.assign(auprc_rank = list(range(0, all_genes_auprc_df.shape[0])))
all_genes_auroc_df = all_genes_auroc_df.assign(auroc_rank = list(range(0, all_genes_auprc_df.shape[0])))

all_genes_auprc_df = all_genes_auprc_df.assign(ras = 0)
all_genes_auprc_df.ix[all_genes_auprc_df.index.isin(ras_genes_df['genes']), 'ras'] = 1

all_genes_metrics_df = all_genes_auprc_df.reset_index().merge(all_genes_auroc_df,
                                                              left_on='index', right_index=True)

all_genes_metrics_df.columns = ['Gene', 'AUPRC', 'AUPRC Rank', 'ras', 'AUROC', 'AUROC Rank']
all_genes_metrics_df.to_csv(all_gene_metrics_file, sep='\t', index=False)
all_genes_metrics_df.head(10)
