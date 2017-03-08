#!/bin/bash

mkdir -p 'data/raw/'

# Download synapse data from Data Freeze Version 1.3
# https://www.synapse.org/#!Synapse:syn4557014
# Must be an approved Synapse user

# Normalized RNAseq data (version 3 - modified on 2016-03-24)
synapse get syn4976369
mv EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv \
data/raw/pancan_normalized_rnaseq.tsv

# Clinical sample freeze info (downloaded 11 Jan 2017)
# Google doc located here
# https://docs.google.com/spreadsheets/d/1Z1H3mXdO_sk9nc0v8df7VNS_XzXiu6vKGJRbK1qYoh4/edit?usp=gmail

# MC3 Mutation data (Must be TCGA Jamboree User - NIHEXT)
# Accessible throgh sftp with USER@tcgaftps.nci.nih.gov
# /tcgajamboree/mc3/mc3.v0.2.8.PUBLIC.maf.gz

# Copy Number data - Thresholded gain/loss calls (Version 1)
synapse get syn5049520
mv all_thresholded.by_genes_whitelisted.tsv data/raw/pancan_GISTIC_threshold.tsv

# Check md5sums of downloaded files
md5sum -c md5sums.txt
