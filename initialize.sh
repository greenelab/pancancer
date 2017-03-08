#!/bin/bash

# Run once to intialize repository and prepare for usage

# Download controlled access data
# Note - currently requires synapse login and dbGAP controlled access
bash scripts/initialize/download_data.sh

# Process RNAseq and Mutation matrix
python scripts/initialze/process_sample_freeze.py

# Process Copy Number matrix
python scripts/initialize/process_copy_number.py

