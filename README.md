# Gene expression machine learning classifiers from TCGA PanCancerAtlas

**Gregory Way and Casey Greene**

## Detecting system-wide changes in whole transcriptomes

A transcriptome can describe the total state of a tumor at a snapshot
in time. In this repository, we use cancer transcriptomes from The Cancer
Genome Atlas PanCancerAtlas project to interrogate gene expression states
induced by deleterious mutations and copy number alterations.

The code in this repository is flexible and can build a Pan-Cancer classifier
for any combination of genes and cancer-types using gene expression, mutation,
and copy number data. In this repository, we provide examples for building
classifiers to detect aberration in _TP53_ and Ras signalling.

### TP53

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.400250.svg)](https://doi.org/10.5281/zenodo.400250)

We are interested in building a classifier to detect _TP53_ inactivation.
_TP53_ is the most highly mutated gene in cancer and regulates several important
tumorigenic processes such as apoptosis and DNA damage response (DDR). We also
include a pipeline to build and evaluate a machine learning _TP53_ classifier.

### Ras Signaling

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1000876.svg)](https://doi.org/10.5281/zenodo.1000876)

The Ras signaling pathway is a major player in cancer development and
treatment resistance. We observed that nearly 60% of all tumors in TCGA
have mutations or copy number alterations in at least one of 38 core pathway
genes. We applied our approach to detect Ras pathway activation using
_KRAS_, _HRAS_, and _NRAS_ gain of function mutations and copy number gains
to define our gold standard Ras hyperactivation events. We train a supervised
classifier to detect when a tumor has activated Ras.

#### Ras signaling classifier identifies phenocopying NF1 loss of function events

We have previously described the ability of a machine learning classifier to
detect an _NF1_ inactivation signature using Glioblastoma data
([Way _et al._ 2016](http://doi.org/10.1186/s12864-017-3519-7)). There, we
applied an ensemble of logistic regression classifiers to the problem, but the
solutions were unstable and overfit. To address these issues, we posited that
we could leverage data from diverse cancer types to build a pancancer _NF1_
classifier. We also hypothesized that a Ras classifier would be able to detect
tumors with _NF1_ inactivation since _NF1_ directly inhibits RAS activity.

## Controlled Access Data

All data used in this analysis are under controlled access by the The National
Institutes of Health (NIH) and The Cancer Genome Atlas (TCGA). All data are
downloaded from [synapse](http://synapse.org) or
[dbGaP](https://www.ncbi.nlm.nih.gov/gap), which require login and access
credentials. To request access contact _SynapseInfo@sagebase.org_ for specific
details and instructions. Additionally, the mutation data requires a TCGA
Jamboree and an eRA commons account.

Eventually, all of the controlled access data used in this pipeline will be
made public. **We will update this database when the data is officially
released.**

## Usage

### Initialization

The pipeline must be initialized before use. Initialization will download and
process data and setup computational environment for R and python.

To initialize, enter the following in the command line:

```sh
# Login to synapse to download controlled-access data
synapse login

# Create and activate conda environment
conda env create --quiet --force --file environment.yml
source activate pancancer-classifier

# Initialize script
./initialize.sh
```

### Example Scripts

We provide two distinct example pipelines for predicting _TP53_ and _NF1_/RAS
loss of function.

1. _TP53_ loss of function (see [tp53_analysis.sh](tp53_analysis.sh))
2. Ras signaling hyperactivation (see [ras_analysis.sh](ras_analysis.sh))

### Customization

For custom analyses, use the
[scripts/pancancer_classifier.py](scripts/pancancer_classifier.py) script with
command line arguments.

```
python scripts/pancancer_classifier.py ...
```

| Flag | Required/Default | Description |
| ---- | :--------------: | ----------- |
| `--genes` | Required |  Build a classifier for the input gene symbols |
| `--diseases` | `Auto` | The disease types to use in building the classifier |
| `--folds` | `5` | Number of cross validation folds |
| `--drop` |  `False` | Decision to drop input genes from expression matrix |
| `--copy_number` |  `False` | Integrate copy number data to gene event |
| `--filter_count` |  `15` | Default options to filter diseases if none are specified |
| `--filter_prop` |  `0.05` | Default options to filter diseases if none are specified |
| `--num_features` |  `8000` | Number of MAD genes used to build classifier |
| `--alphas` | `0.1,0.15,0.2,0.5,0.8,1` | The alpha grid to search over in parameter sweep |
| `--l1_ratios` | `0,0.1,0.15,0.18,0.2,0.3` | The l1 ratio grid to search over in parameter sweep |
| `--alt_genes` | `None` | Alternative genes to test classifier performance |
| `--alt_diseases` |  `Auto` | Alternative diseases to test classifier performance |
| `--alt_filter_count` | `15` | Filtering used for alternative disease classification |
| `--alt_filter_prop` |  `0.05` | Filtering used for alternative disease classification |
| `--alt_folder` | `Auto` | Location to save all classifier figures |
| `--remove_hyper` | `False` | Decision to remove hyper mutated tumors |
| `--keep_intermediate` | `False` | Decision to keep intermediate ROC curve metrics |


