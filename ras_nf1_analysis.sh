#!/bin/bash
#
# Pipeline to reproduce RAS/NF1 classifier
#
# Usage: bash ras_nf1_analysis.sh
#
# Output: will run all specified classifiers which will output performance plots
#         and summarize how a machine learning classifier can detect aberrantly
#         regulated genes using RNAseq gene expression.

alphas='0.1,0.13,0.15,0.18,0.2,0.3,0.4,0.6,0.7'
l1_mixing='0.1,0.125,0.15,0.2,0.25,0.3,0.35'
nf1_diseases='BLCA,COAD,GBM,LGG,LUAD,LUSC,OV,PCPG,SARC,SKCM,STAD,UCEC'
ras_diseases='BLCA,CESC,COAD,ESCA,HNSC,LUAD,LUSC,OV,PAAD,PCPG,READ,SKCM,STAD,TGCT,THCA,UCEC'

# 1. PanCancer NF1 Classification
python scripts/pancancer_classifier.py --genes 'NF1' --drop --copy_number \
        --diseases $nf1_diseases --alphas $alphas --l1_ratios $l1_mixing \
        --remove_hyper

# 2. PanCancer RAS Classification and predict NF1 using RAS classifier
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --diseases $ras_diseases --drop \
        --remove_hyper --copy_number --alphas $alphas --l1_ratios $l1_mixing \
        --alt_genes 'NF1' --alt_diseases $nf1_diseases

# 3. Within Tissue NF1 Classification
python scripts/within_tissue_analysis.py --genes 'NF1' \
        --diseases $nf1_diseases --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing

# 4. Within Tissue RAS Classification
python scripts/within_tissue_analysis.py --genes 'KRAS,HRAS,NRAS' \
        --diseases $ras_diseases --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing
