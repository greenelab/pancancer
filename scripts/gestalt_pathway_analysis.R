# Gregory Way 2017
# PanCancer Classifier
# scripts/gestalt_pathway_analysis.R
#
# Performs Overrepresentation Pathway Analyses (ORA)
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/gestalt_pathway_analysis.R
#
# Output:
# ORA results for Ras Classifier Genes (all genes and positive/negative genes)

library(WebGestaltR)
library(dplyr)

run_webgestalt <- function(genes, output_name, background_file) {
  # Function to run a WebGestalt Pathway Analysis
  #
  # Arguments:
  #    genes: A character of genes to perform ORA test on
  #    output_name: The name of the sub folder to save results in
  #    background_file: A filepath pointing to a .txt file with a single column
  #                     indicating the background to perform the ORA against
  #
  # Output:
  #    Return overrepresented pathways as a dataframe and write full results
  #    to the specified folder (`gestalt/output_name`)

  webgestalt_output <- WebGestaltR(enrichMethod = "ORA",
                                   organism = "hsapiens",
                                   interestGene = genes,
                                   interestGeneType = "genesymbol",
                                   minNum = 4,
                                   fdrMethod = "BH",
                                   is.output = TRUE,
                                   outputDirectory = "gestalt",
                                   referenceGeneFile = background_file,
                                   referenceGeneType = "genesymbol",
                                   projectName = output_name)
  return(webgestalt_output)
}

# Specify filenames and load data
base_path <- file.path("classifiers", "RAS")
bg_file <- file.path(base_path, "functional_analysis_background_genes.txt")
ras_file <- file.path(base_path, "classifier_coefficients.tsv")
ras_gene_df <- readr::read_tsv(ras_file)

# Separate genes
all_genes <- ras_gene_df %>% dplyr::filter(abs > 0) %>% select(feature)
pos_genes <- ras_gene_df %>% dplyr::filter(weight > 0) %>% select(feature)
neg_genes <- ras_gene_df %>% dplyr::filter(weight < 0) %>% select(feature)

# Perform the analysis
all_output <- run_webgestalt(all_genes$feature, 'all_genes_ras', bg_file)
pos_output <- run_webgestalt(pos_genes$feature, 'pos_genes_ras', bg_file)
neg_output <- run_webgestalt(neg_genes$feature, 'neg_genes_ras', bg_file)
