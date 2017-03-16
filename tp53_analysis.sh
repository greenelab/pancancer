#!/bin/bash

# Pipeline to reproduce pancancer pathways TP53 machine learning classifier
#
# Usage: bash tp53_analysis.sh
#
# Output: Will train a pan cancer model to detect TP53 aberration. Will also
#         train a unique classifier within each specific cancer type

# Set Constants
tp53_diseases='BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KICH,LGG,LIHC,LUAD,LUSC,'\
'PAAD,PRAD,READ,SARC,SKCM,STAD,UCEC'
alphas='0.1,0.13,0.15,0.18,0.2,0.3,0.4,0.6,0.7'
l1_mixing='0.1,0.125,0.15,0.2,0.25,0.3,0.35'
tp53_dir='classifiers/TP53'

# Pan Cancer TP53 classification
python scripts/pancancer_classifier.py \
        --genes 'TP53' \
        --diseases $tp53_diseases \
        --drop \
        --copy_number \
        --remove_hyper \
        --alt_folder $tp53_dir \
        --alphas $alphas \
        --l1_ratios $l1_mixing

# Within Disease type TP53 classification
python scripts/within_tissue_analysis.py \
        --genes 'TP53' \
        --diseases $tp53_diseases \
        --remove_hyper \
        --alt_folder $tp53_dir'/within_disease' \
        --alphas $alphas \
        --l1_ratios $l1_mixing

