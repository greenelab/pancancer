#!/bin/bash

# Pipeline to reproduce figures for the DNA Damage Repair manuscript 
#
# Usage: bash scripts/tp53_ddr_figures.sh
#
# Output: summarizes the results of the TP53 classifier and outputs
#         several figures and tables

tp53_dir='classifiers/TP53'

# 1. Apply PanCan classifier to all samples and output scores for each sample
python scripts/apply_weights.py --classifier $tp53_dir --copy_number

# 2. Summarize and visualize performance of classifiers
python scripts/visualize_decisions.py --scores $tp53_dir --custom 'TP53_loss'
python scripts/map_mutation_class.py --scores $tp53_dir --genes 'TP53'
Rscript --vanilla scripts/viz/ddr_summary_figures.R
Rscript --vanilla scripts/compare_within_models.R \
        --within_dir $tp53_dir'/within_disease' --pancan_summary $tp53_dir 

# 3. Perform Snaptron analysis
# NOTE: Snaptron setup must be performed first. See `pancancer/scripts/snaptron/`
bash dna_damage_repair_tp53exon.sh

# 4. Perform copy number burden analysis
python scripts/copy_burden_merge.py --classifier_folder $tp53_dir
Rscript --vanilla scripts/copy_burden_figures.R

