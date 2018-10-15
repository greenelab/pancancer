#!/bin/bash

# Process the data for correct plotting format
python process_tp53_junctions.py

# Extract data for TP53 exon enrichment
jupyter nbconvert --to=html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute process-tp53-exons.ipynb 

# Plot summary figures and significance test
Rscript investigate_silent_junctions.R
