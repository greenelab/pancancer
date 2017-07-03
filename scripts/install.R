# Gregory Way 2017
# PanCancer Classifier
# install.R
#
# Initializes all the R package versioning for use in the PanCancer pipeline
#
# Usage: Run once at the onset of development or reproduction
#
#     Rscript --vanilla scripts/install.R

library("methods")
library("checkpoint")

dir.create(".checkpoint")
checkpoint("2017-06-01", checkpointLocation = ".")

library("pacman")

cran_packages <- c(
  "cowplot",
  "dplyr",
  "ggplot2",
  "gridExtra",
  "magrittr",
  "optparse",
  "pacman",
  "pheatmap",
  "readr",
  "reshape2",
  "tibble"
)

pacman::p_load(cran_packages, install = FALSE, character.only = TRUE)

sink("sessionInfo.txt")
sessionInfo()
sink()
