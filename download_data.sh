#!/bin/bash

mkdir -p 'data/raw/'

# Download synapse data - must be approved user:
# https://www.synapse.org/#!Synapse:syn4557014 

# Normalized RNAseq data
synapse get syn4976369
mv EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv \
data/raw/pancan_normalized_rnaseq.tsv

# Clinical data
synapse get syn4983466
mv clinical_PANCAN_patient_with_followup.tsv data/raw/pancan_clinical.tsv

# Mutation data
synapse get syn6140557
gunzip pancan.merged.v0.2.4.filtered.maf.gz 
mv pancan.merged.v0.2.4.filtered.maf data/raw/pancan_mutation.maf
