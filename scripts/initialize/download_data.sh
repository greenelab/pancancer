#!/bin/bash

mkdir -p 'data/raw/'

# Download synapse data from Data Freeze Version 1.3
# https://www.synapse.org/#!Synapse:syn4557014
# Must be an approved Synapse user

# Normalized RNAseq data (version 3 - modified on 2016-03-24)
synapse get syn4976369
mv EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv \
data/raw/pancan_normalized_rnaseq.tsv

# Clinical sample freeze info (downloaded 7 March 2017)
# Google doc located here
# https://docs.google.com/spreadsheets/d/1Z1H3mXdO_sk9nc0v8df7VNS_XzXiu6vKGJRbK1qYoh4/edit?usp=gmail
# Move to `data/raw/sampleset_freeze_version4_modify.csv`
# Modify sample freeze internally:
# Replace '01' to '02' in the tumor code following four samples:
# 'TCGA-06-0152-02', 'TCGA-06-0171-02', 'TCGA-06-0221-02', 'TCGA-14-0736-02'

# MC3 Mutation data (Must be TCGA Jamboree User - NIHEXT)
# Accessible throgh sftp with USER@tcgaftps.nci.nih.gov
# /tcgajamboree/mc3/mc3.v0.2.8.PUBLIC.maf.gz

# Copy Number data - Thresholded gain/loss calls (Version 1)
synapse get syn5049520
mv all_thresholded.by_genes_whitelisted.tsv data/raw/pancan_GISTIC_threshold.tsv

# Mutation Burden data - number of mutations per Mb for each sample
synapse get syn7994727
mv mutation-load.txt data/mutation-load.txt

# Check md5sums of downloaded files
md5sum -c scripts/initialize/md5sums.txt
