#!/bin/bash

# Comparing classifier performance in predicting NF1 loss of function with features transformed by various algorithms.
# NOTE: This analysis was performed using publicly available data

# Specify location of feature matrices by versioned GitHub URLs
raw='https://github.com/greenelab/tybalt/raw/928804ffd3bb3f9d5559796b2221500c303ed92c/data/pancan_scaled_rnaseq.tsv.gz'
shuffled='https://github.com/gwaygenomics/pancan_viz/raw/7725578eaefe3eb3f6caf2e03927349405780ce5/data/pancan_scaled_rnaseq_shuffled_ids.tsv.gz'
pca='https://github.com/gwaygenomics/pancan_viz/raw/7725578eaefe3eb3f6caf2e03927349405780ce5/data/pca_rnaseq.tsv.gz'
ica='https://github.com/gwaygenomics/pancan_viz/raw/7725578eaefe3eb3f6caf2e03927349405780ce5/data/ica_rnaseq.tsv.gz'
nmf='https://github.com/gwaygenomics/pancan_viz/raw/7725578eaefe3eb3f6caf2e03927349405780ce5/data/nmf_rnaseq.tsv.gz'
adage='https://github.com/greenelab/tybalt/raw/87496e23447a06904bf9c07c389584147b87bd65/data/encoded_adage_features.tsv'
vae='https://github.com/greenelab/tybalt/raw/87496e23447a06904bf9c07c389584147b87bd65/data/encoded_rnaseq_onehidden_warmup_batchnorm.tsv'
vae_twolayer='../tybalt/data/encoded_rnaseq_twohidden_100model.tsv.gz'
vae_twolayer300='../tybalt/data/encoded_rnaseq_twohidden_300model.tsv.gz'

# Set consistent arguments for both TP53/NF1
nf1_args='--genes NF1 --alphas 0.1,0.13,0.15,0.2,0.3 --l1_ratios 0.15,0.155,0.16,0.2,0.25 --drop --copy_number --remove_hyper --keep_inter --y_matrix xena'

# Perform the NF1 Classification Analysis
python scripts/pancancer_classifier.py $nf1_args --x_matrix $raw --alt_folder 'feature_comparison/NF1/raw'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $shuffled --alt_folder 'feature_comparison/NF1/raw_shuffled'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $pca --alt_folder 'feature_comparison/NF1/pca'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $ica --alt_folder 'feature_comparison/NF1/ica'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $nmf --alt_folder 'feature_comparison/NF1/nmf'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $adage --alt_folder 'feature_comparison/NF1/adage'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $vae --alt_folder 'feature_comparison/NF1/tybalt'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $vae_twolayer --alt_folder 'feature_comparison/NF1/tybalt_twohidden'
python scripts/pancancer_classifier.py $nf1_args --x_matrix $vae_twolayer300 --alt_folder 'feature_comparison/NF1/tybalt_twohidden300'

# Visualize the results
Rscript feature_comparison/feature_comparison.R

