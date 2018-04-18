#!/bin/bash

mkdir -p 'data/raw/'

####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download TCGA PanCanAtlas Data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################

####################################
# RNAseq data
####################################
# Normalized RNAseq data (version 3 - modified on 2016-03-24)
# Synapse syn4976369.3 moved to figshare after embargo lifted
wget --output-document=data/raw/pancan_normalized_rnaseq.tsv \
    https://ndownloader.figshare.com/files/11099351

####################################
# MC3 Mutation data
#
# See Ellrott et al. 2018 for more details
# https://doi.org/10.1016/j.cels.2018.03.002
####################################
# MC3 Data Page: https://gdc.cancer.gov/about-data/publications/pancanatlas
# md5sum of the uncompressed file: 91f16c866b069a0b607638fdb0f111de
wget --output-document=data/raw/mc3.v0.2.8.PUBLIC.maf.gz \
    http://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc

####################################
# Copy number data
####################################
# Copy Number data - Thresholded gain/loss calls (Version 1)
wget --output-document=data/raw/pancan_GISTIC_threshold.tsv \
    https://ndownloader.figshare.com/files/11095412

# Segment based scores - measurement of total copy number burden
wget --output-document=data/seg_based_scores.tsv \
    https://ndownloader.figshare.com/files/11095427

####################################
# Checksums
####################################
# Confirm the integrity of downloaded data sets
md5sum -c scripts/initialize/md5sums.txt
