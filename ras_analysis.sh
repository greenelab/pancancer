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

################
# Step 1. PanCancer NF1 Classification
################
# Train using RNAseq matrix
python scripts/pancancer_classifier.py --genes 'NF1' --drop --copy_number \
        --diseases $nf1_diseases --alphas $alphas --l1_ratios $l1_mixing \
        --remove_hyper --alt_folder 'classifiers/NF1' --keep_intermediate \
        --shuffled

################
# Step 2. PanCancer RAS Classification
################
# Train using RNAseq matrix and predict NF1 using Ras classifier
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
        --copy_number --diseases $ras_diseases --alphas $alphas \
        --l1_ratios $l1_mixing --remove_hyper --alt_folder 'classifiers/RAS' \
        --keep_intermediate --shuffled \
        --alt_genes 'NF1' --alt_diseases $nf1_diseases

# Train using shuffled RNAseq matrix
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
        --copy_number --diseases $ras_diseases --alphas $alphas \
        --l1_ratios $l1_mixing --remove_hyper --shuffled_before_training \
        --keep_intermediate --alt_folder 'classifiers/RAS_shuffled'

################
# Step 3. Within Cancer-Type Classification Comparison
################
# NF1
python scripts/within_tissue_analysis.py --genes 'NF1' \
        --diseases $nf1_diseases --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing \
        --alt_folder 'classifiers/NF1/within_disease'

# Summarize NF1 within cancer-type performance
Rscript scripts/compare_within_models.R --pancan_summary 'classifiers/NF1' \
        --within_dir 'classifiers/NF1/within_disease'

# Ras
python scripts/within_tissue_analysis.py --genes 'KRAS,HRAS,NRAS' \
        --diseases $ras_diseases --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing \
        --alt_folder 'classifiers/RAS/within_disease'

# Summarize Ras within cancer-type performance
Rscript scripts/compare_within_models.R --pancan_summary 'classifiers/RAS' \
        --within_dir 'classifiers/RAS/within_disease/' \
        --alt_gene 'classifiers/NF1'

###############
# Step 4. Get scores for all samples and visualize distribution of scores
###############
# NF1
python scripts/apply_weights.py --classifier 'classifiers/NF1' --copy_number
python scripts/visualize_decisions.py --scores 'classifiers/NF1'

# Ras
python scripts/apply_weights.py --classifier 'classifiers/RAS' --copy_number
python scripts/visualize_decisions.py --scores 'classifiers/RAS'

###############
# Step 5. Map mutations to Ras pathway scores
###############
python scripts/map_mutation_class.py --scores 'classifiers/RAS' \
        --genes 'data/ras_genes.csv'

jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute scripts/alternative_genes_pathwaymapper.ipynb

###############
# Step 6. Rerun Ras classifier without THCA and SKCM and perform analysis
#         (BRAFV600E in THCA was not predicted)
###############
ras_no_thca_skcm=${ras_diseases/SKCM,}
ras_no_thca_skcm=${ras_no_thca_skcm/THCA,}

python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' --drop \
        --remove_hyper --copy_number --alphas $alphas --l1_ratio $l1_mixing \
        --diseases $ras_no_thca_skcm --shuffled \
        --alt_folder 'classifiers/RAS_noTHCASKCM' \
        --keep_intermediate

python scripts/apply_weights.py --classifier 'classifiers/RAS_noTHCASKCM' \
       --copy_number
python scripts/map_mutation_class.py --scores 'classifiers/RAS_noTHCASKCM' \
        --genes 'data/ras_genes.csv'

###############
# Step 7. Perform Some Additional Benchmarking Analysis
###############
# Randomly shuffle input RNAseq features and build a classifier
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
        --diseases $ras_diseases --copy_number --remove_hyper \
        --alphas $alphas --l1_ratios $l1_mixing \
        --shuffled_before_training --keep_intermediate \
        --alt_folder 'classifiers/RAS_shuffled_before_training'

# Do not include copy number in the classifier construction
# The shuffled flag here makes classifier predictions on shuffled RNAseq data
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS'\
        --drop \
        --diseases $ras_diseases \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --remove_hyper \
        --shuffled \
        --alt_folder 'classifiers/RAS_nocopy' \
        --keep_intermediate

# Do not include mutation in the classifier construction
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
        --drop \
        --diseases $ras_diseases \
        --copy_number \
        --no_mutation \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --remove_hyper \
        --shuffled \
        --alt_folder 'classifiers/RAS_nomutation' \
        --keep_intermediate

# Drop all Rasopathy genes
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
         --drop \
         --drop_rasopathy \
         --diseases $ras_diseases \
         --copy_number \
         --remove_hyper \
         --alphas $alphas \
         --l1_ratios $l1_mixing \
         --shuffled \
         --alt_folder 'classifiers/RAS_droprasopathy' \
         --keep_intermediate

# Use only covariate information
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
         --drop \
         --diseases $ras_diseases \
         --copy_number \
         --remove_hyper \
         --alphas $alphas \
         --l1_ratios $l1_mixing \
         --drop_expression \
         --alt_folder 'classifiers/RAS_onlycovariate' \
         --keep_intermediate

# Use only gene expression information
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
         --drop \
         --diseases $ras_diseases \
         --copy_number --remove_hyper \
         --alphas $alphas \
         --l1_ratios $l1_mixing \
         --shuffled \
         --drop_covariate \
         --alt_folder 'classifiers/RAS_onlyexpression' \
         --keep_intermediate

# Drop no genes
python scripts/pancancer_classifier.py --genes 'KRAS,HRAS,NRAS' \
         --diseases $ras_diseases \
         --copy_number \
         --remove_hyper \
         --alphas $alphas \
         --l1_ratios $l1_mixing \
         --shuffled \
         --alt_folder 'classifiers/RAS_nodrop' \
         --keep_intermediate

###############
# Step 8. Plot additional Ras, NF1, and BRAF results
###############
# Plot Ras pathway heatmaps
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute scripts/ras_count_heatmaps.ipynb

# Visualize CCLE predictions
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute scripts/ras_cell_line_predictions.ipynb

# Plot summary figures
Rscript --vanilla scripts/viz/ras_summary_figures.R
Rscript --vanilla scripts/viz/ras_ccle_pharmacology.R
Rscript --vanilla scripts/viz/ras_benchmarking_figures.R

Rscript --vanilla scripts/viz/nf1_summary_figures.R
Rscript --vanilla scripts/viz/braf_summary_figures.R
