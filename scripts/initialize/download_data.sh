#!/bin/bash

mkdir -p 'data/raw/'

####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download synapse data from Data Freeze Version 1.3
# https://www.synapse.org/#!Synapse:syn4557014
# Must be an approved Synapse user
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################

####################################
# RNAseq data
####################################
# Normalized RNAseq data (version 3 - modified on 2016-03-24)
synapse get syn4976369
mv EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv \
data/raw/pancan_normalized_rnaseq.tsv

####################################
# Sample freeze info
####################################
# Clinical sample freeze info (downloaded 7 March 2017)
# Google doc located here
# https://docs.google.com/spreadsheets/d/1Z1H3mXdO_sk9nc0v8df7VNS_XzXiu6vKGJRbK1qYoh4/edit?usp=gmail
# Move to `data/raw/sampleset_freeze_version4_modify.csv`

# Modify sample freeze internally:
# Replace '01' to '02' in the tumor code following four samples:
# 'TCGA-06-0152-02', 'TCGA-06-0171-02', 'TCGA-06-0221-02', 'TCGA-14-0736-02'
# These samples are listed as '02', indicating that they are "recurrent" in the
# RNAseq file. In order to match them to their appropriate mutation calls, they
# must be changed

####################################
# Mutation data
####################################
# MC3 Mutation data (Must be TCGA Jamboree User - NIHEXT)
# Accessible through sftp
# NOTE: user must login using: USER@tcgaftps.nci.nih.gov

# The data is located here, with checksums in scripts/initialize/md5sums.txt
# /tcgajamboree/mc3/mc3.v0.2.8.PUBLIC.maf.gz

####################################
# Copy number data
####################################
# Copy Number data - Thresholded gain/loss calls (Version 1)
synapse get syn5049520
mv all_thresholded.by_genes_whitelisted.tsv data/raw/pancan_GISTIC_threshold.tsv

# Segment based scores - measurement of total copy number burden
synapse get syn7415024
mv seg_based_scores.tsv data/seg_based_scores.tsv

####################################
# Checksums
####################################
# Confirm the integrity of downloaded data sets
md5sum -c scripts/initialize/md5sums.txt

