#!/bin/bash

###################
# NF1 Classifier - PanCancer Project
# Gregory Way 2016
###################

# Step 0: Download Data
# NOTE: this requires controlled access to synapse data
./download_data.sh

# Step 1: Describe the data
python subset_datasets.py

# Step 2: Pancancer classifier performance
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
