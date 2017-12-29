#!/bin/bash
#
# Pipeline to reproduce RAS/NF1 classifier
#
# Usage: bash ras_analysis.sh
#
# Output: Results from all classifiers - performance tracking plots and various results

alphas='0.1,0.13,0.15,0.18,0.2,0.25,0.3'
l1_mixing='0.15,0.155,0.16,0.2,0.25,0.3,0.4'
nf1_diseases='BLCA,COAD,GBM,LGG,LUAD,LUSC,OV,PCPG,SARC,SKCM,STAD,UCEC'
ras_diseases='BLCA,CESC,COAD,ESCA,HNSC,LUAD,LUSC,OV,PAAD,PCPG,READ,SKCM,STAD,TGCT,THCA,UCEC'

# 1. PanCancer NF1 Classification
python scripts/pancancer_classifier.py --genes 'NF1' --drop --copy_number \
        --diseases $nf1_diseases --alphas $alphas --l1_ratios $l1_mixing \
        --remove_hyper --shuffled --alt_folder 'classifiers/NF1_shuffled' \
        --keep_intermediate

# 2. PanCancer RAS Classification and predict NF1 using RAS classifier
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
        --diseases $ras_diseases --copy_number --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing \
        --shuffled --alt_genes 'NF1' --alt_diseases $nf1_diseases \
        --alt_folder 'classifiers/RAS_shuffled' \
        --keep_intermediate

# 3. Within cancer-type NF1 Classification
python scripts/within_tissue_analysis.py --genes 'NF1' \
        --diseases $nf1_diseases --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing \
        --alt_folder 'classifiers/NF1/within_disease'

# 4. Within cancer-type RAS Classification
python scripts/within_tissue_analysis.py --genes 'KRAS,HRAS,NRAS' \
        --diseases $ras_diseases --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing \
        --alt_folder 'classifiers/RAS/within_disease'

# 5. Compare within disease type classification
Rscript scripts/compare_within_models.R --pancan_summary 'classifiers/NF1/' \
        --within_dir 'classifiers/NF1/within_disease/'
Rscript scripts/compare_within_models.R --pancan_summary 'classifiers/RAS/' \
        --within_dir 'classifiers/RAS/within_disease/' --alt_gene 'classifiers/NF1'

# 6. Visualize scores
python scripts/apply_weights.py --classifier 'classifiers/NF1' --copy_number
python scripts/visualize_decisions.py --scores 'classifiers/NF1'

python scripts/apply_weights.py --classifier 'classifiers/RAS' --copy_number
python scripts/visualize_decisions.py --scores 'classifiers/RAS'

# 7. Map Mutations to RAS pathway scores
python scripts/map_mutation_class.py --scores 'classifiers/RAS' \
        --genes 'classifiers/RAS/ras_genes.csv'
python scripts/alternative_genes_pathwaymapper.py

# 8. Rerun Ras classifier without THCA and SKCM (BRAFV600E in THCA was not predicted)
ras_no_thca_skcm=${ras_diseases/SKCM,}
ras_no_thca_skcm=${ras_no_thca_skcm/THCA,}

python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
        --remove_hyper --copy_number --alphas $alphas --l1_ratio $l1_mixing \
        --diseases $ras_no_thca_skcm --shuffled --alt_folder 'classifiers/RAS_noTHCASKCM_shuffled' \
        --keep_intermediate

python scripts/apply_weights.py --classifier 'classifiers/RAS_noTHCASKCM' --copy_number
python scripts/map_mutation_class.py --scores 'classifiers/RAS_noTHCASKCM' \
        --genes 'classifiers/RAS/ras_genes.csv'

# 9. Plot Ras, NF1, and BRAF results
python scripts/ras_count_heatmaps.py
Rscript --vanilla scripts/viz/ras_summary_figures.R
Rscript --vanilla scripts/viz/nf1_summary_figures.R
Rscript --vanilla scripts/viz/braf_summary_figures.R

# 10. Additional Benchmarking Analysis
# Randomly shuffle input RNAseq features and build a classifier
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
        --diseases $ras_diseases --copy_number --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing \
        --shuffled_before_training \
        --keep_intermediate --alt_folder 'classifiers/RAS_shuffled_before_training'

# Do not include copy number in the classifier construction
# Note that the shuffled flag here makes classifier predictions on shuffled RNAseq data
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
        --diseases $ras_diseases \
        --alphas $alphas --l1_ratios $l1_mixing --remove_hyper \
        --shuffled --alt_folder 'classifiers/RAS_nocopy' --keep_intermediate

# Do not include mutation in the classifier construction
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
        --diseases $ras_diseases --copy_number --no_mutation \
        --alphas $alphas --l1_ratios $l1_mixing --remove_hyper \
        --shuffled --alt_folder 'classifiers/RAS_nomutation' --keep_intermediate

# Drop all Rasopathy genes
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop --drop_rasopathy \
         --diseases $ras_diseases --copy_number --remove_hyper \
         --alphas $alphas --l1_ratios $l1_mixing --shuffled \
         --alt_folder 'classifiers/RAS_droprasopathy'

# Use only covariate information
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
         --diseases $ras_diseases --copy_number --remove_hyper \
         --alphas $alphas --l1_ratios $l1_mixing \
         --shuffled --keep_intermediate --drop_expression \
         --alt_folder 'classifiers/RAS_onlycovariate'

# Use only gene expression information
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
         --diseases $ras_diseases --copy_number --remove_hyper \
         --alphas $alphas --l1_ratios $l1_mixing \
         --shuffled --keep_intermediate --drop_covariate \
         --alt_folder 'classifiers/RAS_onlyexpression'

# Drop no genes
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
         --diseases $ras_diseases --copy_number --remove_hyper \
         --alphas $alphas --l1_ratios $l1_mixing \
         --shuffled --keep_intermediate \
         --alt_folder 'classifiers/RAS_nodrop' \

