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
# Note, publicly available Xena data is also available for download
synapse login

# Create and activate conda environment
conda env create --quiet --force --file environment.yml
source activate pancancer-classifier

# Reproduce pipeline
./run_pipeline.sh
```

For custom analyses, use the `pancancer_classifier.py` script with command line
arguments.

```
python pancancer_classifier.py ...
```

| Flag | Abbreviation | Required | Description |
| :--: | :----------: | :------: | ----------- |
| `genes` | `-g` | yes |  Build a classifier for the input gene symbols |
| `tissues` | `-t` |  | The tissues to use in building the classifier |
| `drop` | `-d` | | Decision to drop input genes from expression matrix |
| `filter_count` | `-c` |  | Default options to filter tissues if none are specified |
| `filter_prop` | `-p` |  | Default options to filter tissues if none are specified |
| `num_features` | `-f` |  | Number of MAD genes used to build classifier |
| `alphas` | `-a` |  | The alphas to search over in parameter sweep |
| `l1_ratios` | `-l` |  | The l1 ratios to search over in parameter sweep |
| `alt_genes` | `-b` |  | Alternative genes to test classifier performance |
| `alt_tissues` | `-s` |  | Alternative tissues to test classifier performance |
| `alt_tissue_count` | `-i` |  | Filtering used for alternative tissue classification |
| `alt_filter_prop` | `-r` |  | Filtering used for alternative tissue classification |
| `xena` | `-x` |  | If present, use publicly available data for building classifier |

