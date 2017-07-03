#!/bin/bash

# Run once to intialize repository and prepare for usage

# Download controlled access data
# Note - currently requires synapse login and dbGAP controlled access
bash scripts/initialize/download_data.sh

# Install correct R package versions
Rscript --vanilla scripts/install.R

# Process RNAseq and Mutation matrix
python scripts/initialize/process_sample_freeze.py

# Process Copy Number matrix
python scripts/initialize/process_copynumber.py

