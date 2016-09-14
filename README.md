# Building gene expression classifiers using TCGA Pan-Cancer Atlas data

**Gregory Way and Casey Greene**

## Detecting system-wide changes in whole transcriptomes

We have previously described the ability of a machine learning classifier to detect an
NF1 inactivation signature using Glioblastoma data. We applied an ensemble of logistic
regression classifiers to the problem, but the solutions were unstable and overfit. To
address these issues, we posited that we could leverage data from diverse tissue-types
to build a pancancer NF1 classifier. We also predicted that an NF1/RAS classifier would
better key in on NF1 inactivation since NF1 directly inhibits RAS activity and there
are many more examples of samples with RAS mutations.

## Controlled Access Data

All data used in this analysis are under controlled access by the The National Institutes
of Health (NIH) and The Cancer Genome Atlas (TCGA). All data are downloaded from
[synapse](http://synapse.org), which requires login and access credentials. To request
access contact _SynapseInfo@sagebase.org_ for specific details and instructions.

## Usage

To run the entire pipeline enter the following in the command line:

```sh
# Login to synapse to download controlled-access data
synapse login

# Create and activate conda environment
conda env create --quiet --force --file environment.yml
source activate pancancer-classifier

# Reproduce pipeline
./run_pipeline.sh
```

This will output summary statistics and reproduce the analysis for a pancancer RAS
classifier performance in identifying NF1 inactivated Glioblastoma
