# Exon-Exon Junction Analyses using Snaptron Data

**Snaptron experiments by [Wilks _et al.](https://doi.org/10.1101/097881)**

*located at http://snaptron.cs.jhu.edu/*

This pipeline is specific to the DNA Damage Response analysis, but it can
be modified for other use cases

## Setup

In order to setup the downstream analyses, we first must obtain the data.

### Step 1

Clone and navigate to the github repository

```bash
git clone https://github.com/ChristopherWilks/snaptron-experiments.git
cd snaptron-experiments
```

More detailed instructions and specific examples can be found
[here](http://snaptron.cs.jhu.edu/)

### Step 2

Query all TP53 junctions in the TCGA dataset and store in a file

```bash
# From the top directory of "snaptron-experiments"
./qs --region "TP53" --datasrc tcga > tp53_junctions.txt
```

*NOTE:* Move `tp53_junctions.txt` into this folder, 
`pancancer/scripts/snaptron/tp53_junctions.txt`

### Step 3

Download `rail_id` and sample dictionary

```bash
wget http://snaptron.cs.jhu.edu/data/tcga/samples.tsv
```

*NOTE:* Move `samples.tsv` into this folder,
`pancancer/scripts/snaptron/samples.tsv`

## Analysis

We first need to process the data, and then generate the summary plot. Once
all the data is downloaded, to reproduce the analysis, run:

```bash
bash dna_damage_repair_tp53exon.sh
```

