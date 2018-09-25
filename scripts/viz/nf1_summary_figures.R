# Gregory Way 2017
# PanCancer Classifier
# scripts/viz/nf1_summary_figures.R
#
# Visualize summary for NF1 classification scores
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/viz/nf1_summary_figures.R
#
# Output:
# Heatmap for included cancer-types in NF1 classifier

library(dplyr)
library(pheatmap)
library(ggplot2)
library(readr)
library(cowplot)
library(gridExtra)
source(file.path("scripts", "util", "pancancer_util.R"))

results_folder <- file.path("classifiers", "NF1")
results <- parse_summary(file.path(results_folder, "classifier_summary.txt"))

# 1) Heatmap of the distribution of aberrant events across tumors
heatmap_plot_file <- file.path(results_folder, "figures", "nf1_heatmap.pdf")
heat_file <- file.path(results_folder, "summary_counts.csv")
heat_df <- readr::read_csv(heat_file)

prop_matrix <- as.matrix(heat_df[, c('NF1_loss_y', 'NF1_y')])
rownames(prop_matrix) <- heat_df$DISEASE
colnames(prop_matrix) <- c("Loss", "Mutation")

# All diseases that are used in building the classifier
nf1_dis <- results[["Diseases"]]

# Build a vector for heatmap labels
classifier <- c()
for (disease in rownames(prop_matrix)) {
  if (disease %in% nf1_dis) {
    classifier <- c(classifier, "Training")
  } else {
    classifier <- c(classifier, "Dropped")
  }
}

classifier <- data.frame(classifier)
rownames(classifier) <- rownames(prop_matrix)
classifier$classifier <- factor(classifier$classifier,
                                levels = c("Training", "Dropped"))

prop_matrix <- prop_matrix[names(sort(prop_matrix[,2], decreasing = TRUE)), ]

# Plot and save heatmap
pheatmap(t(prop_matrix * 100), scale = "none", cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.0f", fontsize_number = 8,
         number_color = "black", annotation_col = classifier,
         annotation_names_col = FALSE, legend = FALSE,
         filename = heatmap_plot_file,
         width = 8, height = 2)
