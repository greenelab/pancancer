#!/bin/bash

###################
# NF1 Classifier - PanCancer Project
# Gregory Way 2016
###################

# Step 0: Download Data
# NOTE: this requires controlled access to synapse data
./download_data.sh

# NOTE: UCSC Xena data is publicly available
pyton download_xena_data.py

# Step 1: Describe the data
python subset_datasets.py --database 'Synapse'
python subset_datasets.py --database 'Xena'

# Step 2: Example Pancancer Classifiers
mkdir -p 'figures/tissue' 'classifiers'

# Analysis 1 - PanCancer RAS classifier to predict NF1 in GBM
python pancancer_classifier.py \
        --genes 'KRAS,NRAS,HRAS' \
        --tissues 'BLCA,COAD,HNSC,LUAD,PCPG,SKCM,CESC,LAML,PAAD,READ,TGCT,THCA' \
        --alt_genes 'NF1'

# Analysis 2 - PanCancer RAS classifier to predict NF1 in GBM - Auto tissue
python pancancer_classifier.py \
        --genes 'KRAS,NRAS,HRAS' \
        --alt_genes 'NF1'

# Analysis 3 - PanCancer NF1 classifier
python pancancer_classifier.py \
        --genes 'NF1'

# Analysis 4 - PanCancer NF1 classifier - Remove NF1 from expression matrix
python pancancer_classifier.py \
        --genes 'NF1' \
        --drop

# Analysis 5 - PanCancer NF1/RAS combined classifier
python pancancer_classifier.py \
        --genes 'NF1,KRAS,NRAS,HRAS' \
        --alt_genes 'NF1'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyses using UCSC Xena Data
# The flag "--xena" will use publicly available data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Analysis 6 - PanCancer NF1 classifier using UCSC Xena Data
python pancancer_classifier.py \
        --genes 'NF1' \
        --tissues 'BRCA,GBM,LGG,LUSC,BLCA,COAD,HNSC,LUAD,PCPG,SKCM' \
        --xena

# Analysis 7 - PanCancer NF1/RAS combined classifier using UCSC Xena Data
python pancancer_classifier.py \
        --genes 'NF1,KRAS,NRAS,HRAS' \
        --alt_genes 'NF1' \
        --alt_filter_count '10' \
        --alt_filter_prop '0.05' \
        --xena

# Analysis 8 - PanCancer NF1/RAS classifier predicts GBM NF1 with preselected tissues
python pancancer_classifier.py \
        --genes 'KRAS,HRAS,NRAS,NF1' \
        --tissues 'CESC,LAML,PAAD,READ,TGCT,THCA,BRCA,LGG,LUSC,BLCA,COAD,HNSC,LUAD,PCPG' \
        --alt_filter_count '10' \
        --alt_filter_prop '0.05' \
        --alt_tissues 'GBM' \
        --alt_genes 'NF1' \
        --xena
