#!/bin/bash

# Pipeline to reproduce pancancer pathways TP53 machine learning classifier
#
# Usage: bash tp53_analysis.sh
#
# Output: will run all specified classifiers which will output performance plots
#         and summarize how a machine learning classifier can detect aberrant
#         TP53 activity RNAseq, copy number, and gene expression.

# Set Constants
tp53_diseases='BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KICH,LGG,LIHC,LUAD,LUSC,OV,'\
'PAAD,PRAD,READ,SARC,SKCM,STAD,UCEC'
alphas='0.1,0.13,0.15,0.18,0.2,0.3,0.4,0.6,0.7'
l1_mixing='0.1,0.125,0.15,0.2,0.25,0.3,0.35'
tp53_dir='classifiers/TP53'

# 1. PanCancer TP53 classification
python scripts/pancancer_classifier.py --genes 'TP53' --diseases $tp53_diseases \
        --drop --copy_number --remove_hyper --alt_folder $tp53_dir \
        --alphas $alphas --l1_ratios $l1_mixing

# 2. Within disease type TP53 classification
python scripts/within_tissue_analysis.py --genes 'TP53' \
        --diseases $tp53_diseases --remove_hyper \
        --alt_folder $tp53_dir'/within_disease' \
        --alphas $alphas --l1_ratios $l1_mixing

# 3. Prepare heatmap figure describing TP53 mutations across the set
python scripts/build_heatmap_data.py --genes 'TP53' --copy_number \
        --remove_hyper --alt_folder $tp53_dir

# 4. Apply PanCan classifier to all samples and output scores for each sample
python scripts/apply_weights.py --classifier $tp53_dir --copy_number

# 5. Summarize and visualize performance of classifiers
python scripts/visualize_decisions.py --scores $tp53_dir --custom 'TP53_loss'
python scripts/map_mutation_class.py --scores $tp53_dir --genes 'TP53'
Rscript --vanilla scripts/ddr_summary_figures.R
Rscript --vanilla scripts/compare_within_models.R \
        --within_dir $tp53_dir'/within_disease' --pancan_summary $tp53_dir 

