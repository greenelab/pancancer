# Building gene expression classifiers using TCGA Pan-Cancer Atlas data

**Gregory Way and Casey Greene**

## Detecting system-wide changes in whole transcriptomes

We have previously described the ability of a machine learning classifier to
detect an NF1 inactivation signature using Glioblastoma data
([Way _et al._ 2016](http://doi.org/10.1186/s12864-017-3519-7)). We applied an
ensemble of logistic regression classifiers to the problem, but the solutions were
unstable and overfit. To address these issues, we posited that we could leverage
data from diverse tissue-types to build a pancancer NF1 classifier. We also
predicted that an RAS classifier would be able to detect tumors NF1 inactivation
since NF1 directly inhibits RAS activity and there are many more examples of
samples with RAS mutations.

The code in this repository is flexible and can build a Pan-Cancer classifier
for any combination of genes and cancer types using gene expression, mutation,
and copy number data. Currently, we build classifiers to detect NF1/RAS
aberration and TP53 inactivation.

## Controlled Access Data

All data used in this analysis are under controlled access by the The National
Institutes of Health (NIH) and The Cancer Genome Atlas (TCGA). All data are
downloaded from [synapse](http://synapse.org), which requires login and access
credentials. To request access contact _SynapseInfo@sagebase.org_ for specific
details and instructions. Additionally, the mutation data requires a TCGA
Jamboree and an eRA commons account.

Eventually, all of the controlled access data used in this pipeline will be
made public. **We will update this database when the data is officially
released.**

We also provide the ability to build classifiers using publicly available
data retrieved from [UCSC Xena](xena.ucsc.ecu).

## Cancer Genes

Note that in order to use the copy number integration feature, an additional
file must be downloaded. The file is `Supplementary Table S2` of
[Vogelstein _et al._ 2013]("http://doi.org/10.1126/science.1235122"). It can
be downloaded [here]("science.sciencemag.org/highwire/filestream/594203/field_highwire_adjunct_files/1/1235122TablesS1-4.xlsx").

Once the file is downloaded, move it to `data/vogelstein_cancergenes.tsv`

To download file:

```bash
cancer_genes_url="http://science.sciencemag.org/highwire/filestream/594203/field_highwire_adjunct_files/1/1235122TablesS1-4.xlsx"
wget $cancer_genes_url --output-document="data/vogelstein_cancergenes.xlsx"
```

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

| Flag | Abbreviation | Required/Default | Description |
| ---- | :----------: | :------: | ----------- |
| `genes` | `-g` | REQUIRED |  Build a classifier for the input gene symbols |
| `tissues` | `-t` | `Auto` | The tissues to use in building the classifier |
| `folds` | `-f` | `5` | Number of cross validation folds |
| `drop` | `-d` | `False` | Decision to drop input genes from expression matrix |
| `copy_number` | `-u` | `False` | Integrate copy number data to gene event |
| `filter_count` | `-c` | `15` | Default options to filter tissues if none are specified |
| `filter_prop` | `-p` | `0.05` | Default options to filter tissues if none are specified |
| `num_features` | `-n` | `8000` | Number of MAD genes used to build classifier |
| `alphas` | `-a` | `0.01,0.1,0.15,0.2,0.5,0.8` | The alpha grid to search over in parameter sweep |
| `l1_ratios` | `-l` | `0,0.1,0.15,0.18,0.2,0.3` | The l1 ratio grid to search over in parameter sweep |
| `alt_genes` | `-b` | `None` | Alternative genes to test classifier performance |
| `alt_tissues` | `-s` | `Auto` | Alternative tissues to test classifier performance |
| `alt_tissue_count` | `-i` | `15` | Filtering used for alternative tissue classification |
| `alt_filter_prop` | `-r` | `0.05` | Filtering used for alternative tissue classification |
| `alt_folder` | `-o` | `Auto` | Location to save all classifier figures |
| `xena` | `-x` | `False` | If present, use publicly available data for building classifier |

