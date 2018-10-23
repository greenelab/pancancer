"""
Gregory Way 2018
PanCancer Classifier
scripts/compression.py

Perform a series of compression algorithms on the PanCanAtlas RNAseq data. The
compressed features (z) are input as features for gene alteration predictions.

Usage:
Run in command line

        python scripts/compression.py

Output:
Compressed gene expression features in the `feature_comparison/data` folder
"""

import os
import numpy as np
import pandas as pd
from statsmodels.robust.scale import mad

from tybalt.data_models import DataModel

np.random.seed(123)

# Load constants
num_genes_kept = 8000
num_components = 100

vae_epochs = 100
vae_batch_size = 150
vae_lr = 0.001

dae_epochs = 100
dae_batch_size = 50
dae_lr = 0.0005
dae_noise = 0
dae_sparsity = 0

algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']

# Load and Process Data
expr_file = os.path.join('..', 'data', 'pancan_rnaseq_freeze.tsv.gz')
rnaseq_df = pd.read_table(expr_file, index_col=0)

# Subset x matrix to MAD genes
med_dev = pd.DataFrame(mad(rnaseq_df), index=rnaseq_df.columns)
mad_genes = (
    med_dev.sort_values(by=0, ascending=False)
    .iloc[0:num_genes_kept]
    .index
    .tolist()
)

rnaseq_df = rnaseq_df.loc[:, mad_genes]

# Initialize DataModel class with PanCanAtlas RNAseq
dm = DataModel(df=rnaseq_df)

# Transform the input matrix into a range between zero and one
dm.transform(how='zeroone')

# Fit models
dm.pca(n_components=num_components)
dm.ica(n_components=num_components)
dm.nmf(n_components=num_components)

dm.nn(n_components=num_components,
      model='tybalt',
      loss='binary_crossentropy',
      epochs=int(vae_epochs),
      batch_size=int(vae_batch_size),
      learning_rate=float(vae_lr),
      separate_loss=False,
      verbose=False)

dm.nn(n_components=num_components,
      model='adage',
      loss='binary_crossentropy',
      epochs=int(dae_epochs),
      batch_size=int(dae_batch_size),
      learning_rate=float(dae_lr),
      noise=float(dae_noise),
      sparsity=float(dae_sparsity),
      verbose=False)

# Output compressed features to files
pca_file = os.path.join('data', 'pca_pancanatlas_z100.tsv.gz')
dm.pca_df.to_csv(pca_file, sep='\t', compression='gzip')

ica_file = os.path.join('data', 'ica_pancanatlas_z100.tsv.gz')
dm.ica_df.to_csv(ica_file, sep='\t', compression='gzip')

nmf_file = os.path.join('data', 'nmf_pancanatlas_z100.tsv.gz')
dm.nmf_df.to_csv(nmf_file, sep='\t', compression='gzip')

dae_file = os.path.join('data', 'dae_pancanatlas_z100.tsv.gz')
dm.adage_df.to_csv(dae_file, sep='\t', compression='gzip')

vae_file = os.path.join('data', 'vae_pancanatlas_z100.tsv.gz')
dm.tybalt_df.to_csv(vae_file, sep='\t', compression='gzip')
