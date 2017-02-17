#!/bin/bash
#
# Pipeline to reproduce pancancer pathways machine learning classifier
#
# Usage: bash pancancer_pathways_analysis.sh
#
# Output: will run all specified classifiers which will output performance plots
#         and summarize how a machine learning classifier can detect aberrantly
#         regulated genes using RNAseq gene expression.

mixing_params='0,0.1,0.2,0.3,0.35,0.4,0.5'
nf1_diseases='BLCA,COAD,GBM,LGG,LUAD,LUSC,OV,PCPG,SARC,SKCM,STAD,UCEC'
ras_diseases='BLCA,CESC,COAD,ESCA,HNSC,LUAD,LUSC,OV,PAAD,PCPG,READ,SKCM,STAD,TGCT,THCA,UCEC'

# 1. PanCancer NF1 Classification
python pancancer_classifier.py --genes 'NF1' --drop --copy_number --tissues $nf1_diseases --l1_ratios $mixing_params

# 2. PanCancer RAS Classification and predict NF1 using RAS classifier
python pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop --alt_genes 'NF1' --copy_number --tissues $ras_diseases --l1_ratios $mixing_params

# 3. Within Tissue NF1 Classification
python within_tissue_analysis.py --genes 'NF1' --diseases $nf1_diseases --l1_ratios $mixing_params

# 4. Within Tissue RAS Classification
python within_tissue_analysis.py --genes 'KRAS,HRAS,NRAS' --diseases $ras_diseases --l1_ratios $mixing_params
