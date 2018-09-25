# Gene expression machine learning classifiers from TCGA PanCancerAtlas

**Gregory Way and Casey Greene**

## Detecting system-wide changes in whole transcriptomes

A transcriptome can describe the total state of a tumor at a snapshot in time.
In this repository, we use cancer transcriptomes from The Cancer Genome Atlas PanCancerAtlas project to interrogate gene expression states induced by deleterious mutations and copy number alterations.

The code in this repository is flexible and can build a Pan-Cancer classifier for any combination of genes and cancer-types using gene expression, mutation,
and copy number data.
In this repository, we provide examples for building classifiers to detect aberration in _TP53_ and Ras signalling.

### TP53

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.400250.svg)](https://doi.org/10.5281/zenodo.400250)

We are interested in building a classifier to detect _TP53_ inactivation.
_TP53_ is the most highly mutated gene in cancer and regulates several important oncogenic processes such as apoptosis and DNA damage response (DDR).
We include a pipeline to build and evaluate a machine learning _TP53_ classifier.
See [`tp53_analysis.sh`](tp53_analysis.sh) for more details.

### Ras Signalling

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1000876.svg)](https://doi.org/10.5281/zenodo.1000876)

The Ras signalling pathway is a major player in cancer development and treatment resistance.
We observed that nearly 60% of all tumors in TCGA have mutations or copy number alterations in at least one of 38 core pathway genes ([Sanchez-Vega et al. 2018](https://doi.org/10.1016/j.cell.2018.03.035)).

We applied our approach to detect Ras pathway activation using _KRAS_, _HRAS_, and _NRAS_ gain of function mutations and copy number gains to define our gold standard Ras hyperactivation events.
We train a supervised classifier to detect when a tumor has activated Ras.

#### Ras signalling classifier identifies phenocopying NF1 loss of function events

We have previously described the ability of a machine learning classifier to detect an _NF1_ inactivation signature using Glioblastoma data ([Way _et al._ 2016](http://doi.org/10.1186/s12864-017-3519-7)).
There, we applied an ensemble of logistic regression classifiers to the problem, but the solutions were unstable and overfit.
To address these issues, we posited that we could leverage data from diverse cancer types to build a pancancer _NF1_ classifier.
We also hypothesized that a Ras classifier would be able to detect tumors with _NF1_ inactivation since _NF1_ directly inhibits RAS activity.

## Open Access Data

All data was released by the TCGA PanCancerAtlas project.
The compendium of papers is described [here](https://www.cell.com/pb-assets/consortium/pancanceratlas/pancani3/index.html).
Supplementary data from these papers can be downloaded from the [NCI](https://gdc.cancer.gov/about-data/publications/pancanatlas).
The specific data used in the analyses presented here are archived on Zenodo
[Gene expression](https://figshare.com/articles/TCGA_PanCanAtlas_Gene_Expression_Data/6146519)) and [copy number](https://figshare.com/articles/TCGA_PanCanAtlas_Copy_Number_Data/6144122) data can be accessed here.

See [`scripts/initialize/download_data.sh`](scripts/initialize/download_data.sh) for more details.

## Usage

### Initialization

The repository must be cloned onto local machine before analyses can proceed.

```sh
# Make sure git-lfs (https://git-lfs.github.com/) is installed before cloning
# If not, run `git lfs install`
git clone git@github.com:greenelab/pancancer.git

cd pancancer
```

### Example Scripts

We provide two distinct example pipelines for predicting _TP53_ and _NF1_/RAS loss of function.

1. _TP53_ loss of function (see [`tp53_analysis.sh`](tp53_analysis.sh))
2. Ras signaling hyperactivation (see [`ras_analysis.sh`](ras_analysis.sh))

### Customization

For custom analyses, use the [`scripts/pancancer_classifier.py`](scripts/pancancer_classifier.py) script with command line arguments.

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
| `--x_matrix` | `raw` | if not "raw", then the filename storing the features |
| `--shuffled` | `False` | Shuffle the X matrix for better training |
| `--shuffled_before_training` | `False` | Remove correlational structure in the data |
| `--no_mutation` | `True` | Decision to remove mutation data from the input matrix |
| `--drop_rasopathy` | `False` | Decision to drop all rasopathy genes from the X matrix |
| `--drop_covariates` | `False` | Decision to drop all covariate information from the X matrix|
