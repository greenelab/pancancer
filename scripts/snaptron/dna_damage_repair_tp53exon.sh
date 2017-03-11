#!/bin/bash

# Process the data for correct plotting format
python process_tp53_junctions.py

# Plot summary figures and significance test
Rscript investigate_silent_junctions.R

