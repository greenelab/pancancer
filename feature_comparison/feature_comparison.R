# Comparing NF1 classification performance using various feature encodings
# determined by different dimensionality reduction algorithms
# Gregory Way 2017
#
# Output:
# ROC curves for predicting NF1 inactivation status (mutation or deep copy
# number loss) using different features

library(dplyr)
library(ggplot2)

source(file.path("scripts", "util", "pancancer_util.R"))

getFeatureResults <- function(results_directory, features, gene) {
  # Process results - will obtain a dataframe for outputing an ROC curve,
  # the cross validation AUROC, and will list the number of nonzero feature
  #
  # Arguments:
  # results_directory - the location of the base directory where results are
  # features - a string indicating the the dimensionality reduction algorithm
  # gene - a string indicating gene name
  #
  # Output:
  # A list with ROC results and a dataframe with nonzero count and algorithm,
  # and cross validation AUROC

  roc_file <- file.path(results_directory, "pancan_roc_results.tsv")
  feature_file <- file.path(results_directory, "classifier_coefficients.tsv")
  summary_file <- file.path(results_directory, "classifier_summary.txt")
  
  # Obtain the ROC for the given file and process dataframe
  roc_ <- readr::read_tsv(roc_file)
  roc_ <- roc_ %>% dplyr::mutate(Features = features)
  roc_ <- roc_ %>% dplyr::mutate(gene = gene)
  
  # Obtain the feature coefficients file
  coef_ <- readr::read_tsv(feature_file)
  num_non_zero <- coef_ %>% dplyr::filter(abs > 0) %>% dplyr::tally()
  colnames(num_non_zero) <- "nonzero"
  num_non_zero$Features <- features

  # Get the classifier summary
  classifier_summary <- parse_summary(summary_file)
  auroc <- classifier_summary$`Cross Validation AUROC`
  auroc <- as.data.frame(round(as.numeric(auroc) * 100, 1))
  colnames(auroc) <- "AUROC (%)"
  auroc$Features <- features
  
  return(list(roc_, num_non_zero, auroc))
}

# Set file path objects
base_file <- file.path("feature_comparison", "NF1")

raw_dir <- file.path(base_file, "raw")
shu_dir <- file.path(base_file, "raw_shuffled")
pca_dir <- file.path(base_file, "pca")
ica_dir <- file.path(base_file, "ica")
nmf_dir <- file.path(base_file, "nmf")
adage_dir <- file.path(base_file, "adage")
vae_dir <- file.path(base_file, "tybalt")
vae_2h_dir <- file.path(base_file, "tybalt_twohidden")
vae_2h300_dir <- file.path(base_file, "tybalt_twohidden300")

# Process algorithm results
raw_results <- getFeatureResults(raw_dir, "Raw", "NF1")
shu_results <- getFeatureResults(shu_dir, "Shuffled", "NF1")
pca_results <- getFeatureResults(pca_dir, "PCA", "NF1")
ica_results <- getFeatureResults(ica_dir, "ICA", "NF1")
nmf_results <- getFeatureResults(nmf_dir, "NMF", "NF1")
adage_results <- getFeatureResults(adage_dir, "ADAGE", "NF1")
vae_results <- getFeatureResults(vae_dir, "Tybalt", "NF1")
vae_2h_results <- getFeatureResults(vae_2h_dir, "VAE (100)", "NF1")
vae_2h300_results <- getFeatureResults(vae_2h300_dir, "VAE (300)", "NF1")

# Combine DataFrames
roc_df <- dplyr::bind_rows(raw_results[[1]],
                           shu_results[[1]],
                           pca_results[[1]],
                           ica_results[[1]],
                           nmf_results[[1]],
                           adage_results[[1]],
                           vae_results[[1]],
                           vae_2h_results[[1]],
                           vae_2h300_results[[1]])

count_df <- dplyr::bind_rows(raw_results[[2]],
                             shu_results[[2]],
                             pca_results[[2]],
                             ica_results[[2]],
                             nmf_results[[2]],
                             adage_results[[2]],
                             vae_results[[2]],
                             vae_2h_results[[2]],
                             vae_2h300_results[[2]])

auroc_df <- dplyr::bind_rows(raw_results[[3]],
                             shu_results[[3]],
                             pca_results[[3]],
                             ica_results[[3]],
                             nmf_results[[3]],
                             adage_results[[3]],
                             vae_results[[3]],
                             vae_2h_results[[3]],
                             vae_2h300_results[[3]])

# Reorder auroc_df
auroc_df <- auroc_df %>% dplyr::select(Features, dplyr::everything())

# Order Factors
factor_levels <-  c("Raw", "Shuffled", "PCA", "ICA", "NMF", "ADAGE",
                    "Tybalt", "VAE (100)", "VAE (300)")

count_df$Features <- factor(count_df$Features, levels = factor_levels)
roc_df$Features <- factor(roc_df$Features, levels = factor_levels)
auroc_df$Features <- factor(auroc_df$Features, levels = factor_levels)

# Subset only to pancancer and also only cross validation performance
roc_df <- roc_df %>%
  dplyr::filter(disease == "PanCan") %>%
  dplyr::filter(train_type == "cv")

# Set plotting colors
manual_colors <- c("#ff7f00", "#cab2d6", "#a6cee3", "#1f78b4", "#b2df8a",
                   "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f")

# Plot ROC curves and feature bar plot
table_theme <- gridExtra::ttheme_default(base_size = 6,
                                         padding = unit(c(2.2, 2.2), "mm"))
table_gg <- gridExtra::tableGrob(auroc_df, rows = NULL, theme = table_theme)
feat_gg <- ggplot(roc_df, aes(x = fpr, y = tpr, color = Features)) +
  geom_step() +
  geom_segment(aes(x = 0 , y = 0, xend = 1, yend = 1),
               linetype = "dashed", color = "black") +
  annotation_custom(table_gg,
                    xmin = 0.75, xmax = 0.85, ymin = 0.25, ymax = 0.35) +
  scale_color_manual(values = manual_colors) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  coord_fixed() + 
  theme_bw() +
  theme(axis.text = element_text(size = rel(0.8)),
        axis.title = element_text(size = rel(1.1)),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0.2, 0.2), "cm"))

coef_gg <- ggplot(count_df, aes(y = nonzero, x = Features,
                                          fill = Features)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = manual_colors) +
  xlab("") +
  ylab("Nonzero Feature Count") + 
  geom_text(aes(label = nonzero, y = nonzero), stat = "identity",
            vjust = -0.2, size = rel(1.9)) +
  theme_bw() +
  ylim(c(0,  max(count_df$nonzero) + 25)) +
  theme(axis.text = element_text(size = rel(0.8)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(0.7)),
        legend.title = element_text(size = rel(0.9)),
        legend.position = "none",
        strip.text = element_text(size = rel(0.9)),
        plot.margin = unit(c(0.25, 0, 0, 0), "cm"))

feature_plot <- cowplot::plot_grid(feat_gg, coef_gg,
                                   labels = c("A", "B"),
                                   scale = c(1, 0.85))

figure_file <- file.path("feature_comparison", "nf1_prediction_comparison.png")
cowplot::save_plot(figure_file, feature_plot, base_aspect_ratio = 1.1,
                   base_height = 3.5, base_width = 7)
