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

###############
# Step 1. Pan Cancer TP53 classification
###############
python scripts/pancancer_classifier.py \
        --genes 'TP53' \
        --diseases $tp53_diseases \
        --drop \
        --copy_number \
        --remove_hyper \
        --alt_folder $tp53_dir \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --keep_intermediate \
        --shuffled

###############
# Step 2. Within Disease type TP53 classification
###############
python scripts/within_tissue_analysis.py \
        --genes 'TP53' \
        --diseases $tp53_diseases \
        --remove_hyper \
        --alt_folder $tp53_dir'/within_disease' \
        --alphas $alphas \
        --l1_ratios $l1_mixing

###############
# Step 3. Get scores for all samples and visualize distribution of scores
###############
python scripts/apply_weights.py \
        --classifier $tp53_dir \
        --copy_number

python scripts/visualize_decisions.py \
        --scores $tp53_dir \
        --custom 'TP53_loss'

python scripts/map_mutation_class.py \
        --scores $tp53_dir \
        --genes 'TP53'

python scripts/copy_burden_merge.py \
        --classifier_folder $tp53_dir

###############
# Step 4. Plot additional TP53 results
###############
# Summary Figures
Rscript --vanilla scripts/viz/ddr_summary_figures.R
Rscript --vanilla scripts/compare_within_models.R \
        --within_dir $tp53_dir'/within_disease' \
        --pancan_summary $tp53_dir

# Mutation classification stratified by cancer-Type
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute scripts/tp53_phenocopy.ipynb

# Mutation classification stratified by phenocopying variant
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute scripts/tp53_ddr_cancertype_subtypes.ipynb

###############
# Step 5. SNAPTRON exon-exon junction analysis
###############
cd scripts/snaptron
bash dna_damage_repair_tp53exon.sh
cd ../..

# Copy burden analysis requires snaptron results
Rscript --vanilla scripts/copy_burden_figures.R
