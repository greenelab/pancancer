"""
Gregory Way 2017=8
PanCancer Classifier
feature_comparison/analysis_pipeline.py

Master script to submit all pancancer classifiers trained on compressed gene
expression features.

Usage:

    python feature_comparison/analysis_pipeline.py

Output:

    Classifier results in the `feature_comparison/classifiers` folder
"""

import os
import subprocess

# Location of the data
raw = ('raw', 'data/pancan_rnaseq_freeze.tsv.gz')
pca = ('pca', 'feature_comparison/data/pca_pancanatlas_z100.tsv.gz')
ica = ('ica', 'feature_comparison/data/ica_pancanatlas_z100.tsv.gz')
nmf = ('nmf', 'feature_comparison/data/nmf_pancanatlas_z100.tsv.gz')
dae = ('dae', 'feature_comparison/data/dae_pancanatlas_z100.tsv.gz')
vae = ('vae', 'feature_comparison/data/vae_pancanatlas_z100.tsv.gz')

# Save Constants
genes = ['TP53', 'KRAS,HRAS,NRAS', 'NF1']
l1_ratios = '0.15,0.155,0.16,0.2,0.25,0.3,0.4'
alphas = '0.1,0.13,0.15,0.18,0.2,0.25,0.3'
alt_base = 'feature_comparison/classifiers'

for gene in genes:
    for alg, x_matrix in [raw, pca, ica, nmf, dae, vae]:
        alt_dir = os.path.join(alt_base, gene, alg)

        if gene == 'KRAS,HRAS,NRAS':
            alt_dir = os.path.join(alt_base, 'RAS', alg)

        command = ['python', 'scripts/pancancer_classifier.py',
                   '--genes', gene,
                   '--x_matrix', x_matrix,
                   '--alt_folder', alt_dir,
                   '--alphas', alphas,
                   '--l1_ratios', l1_ratios,
                   '--drop',
                   '--copy_number',
                   '--remove_hyper',
                   '--keep_inter']

        print(command)

        subprocess.call(command)
