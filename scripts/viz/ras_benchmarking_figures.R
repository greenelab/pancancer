# Gregory Way 2017
# PanCancer Classifier
# scripts/viz/ras_benchmarking_figures.R
#
# Compares performance of Ras classifier in 3 different scenarios
#
# 1) Mutation vs. Copy Number vs. Both
# 2) Dropping Ras vs. Dropping Rasopathy vs. no Drop
# 3) Expression vs. Covariates vs. Both
#
# Also compares classifier coefficients across for comparisons 1 and 2
# (none of the coefficients are the same for 3)
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/viz/ras_benchmarking_figures.R
#
# Output:
# 5 Plots given above

library(dplyr)
library(ggplot2)
library(ggrepel)

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
  coef_ <- coef_ %>% dplyr::mutate(Features = features)
  coef_ <- coef_ %>% dplyr::mutate(gene = gene)

  # Get the classifier summary
  classifier_summary <- parse_summary(summary_file)
  train_auroc <- classifier_summary$`Training AUROC`
  test_auroc <- classifier_summary$`Testing AUROC`
  train_auroc <- as.data.frame(round(as.numeric(train_auroc) * 100, 1))
  test_auroc <- as.data.frame(round(as.numeric(test_auroc) * 100, 1))
  auroc <- cbind(train_auroc, test_auroc)
  colnames(auroc) <- c("Train", "Test")
  auroc$Features <- features
  
  return(list(roc_, coef_, auroc))
}

# Set file path objects
base_folder <- 'classifiers'
base_ras <- file.path(base_folder, 'RAS_shuffled')

ras_mutation_only <- file.path(base_folder, 'RAS_nocopy')
ras_copy_only <- file.path(base_folder, 'RAS_nomutation')

ras_expression_only <- file.path(base_folder, 'RAS_onlyexpression')
ras_covariate_only <- file.path(base_folder, 'RAS_onlycovariate')

ras_nodrop <- file.path(base_folder, "RAS_nodrop")
ras_droprasopathy <- file.path(base_folder, "RAS_droprasopathy")

set.seed(123)

# Process algorithm results
raw_results <- getFeatureResults(base_ras, "Combined", "Ras")
mutation_only_results <- getFeatureResults(ras_mutation_only, 'Mutation Only',
                                           'Ras')
copy_only_results <- getFeatureResults(ras_copy_only, "Copy Only", "Ras")

ras_expression_only_results <- getFeatureResults(ras_expression_only,
                                                 "Exprs Only", "Ras")
ras_covariate_only_results <- getFeatureResults(ras_covariate_only,
                                                "Cov Only", "Ras")

ras_nodrop_results <- getFeatureResults(ras_nodrop, "No Gene Drop", "Ras")
ras_droprasopathy_results <- getFeatureResults(ras_droprasopathy,
                                               "Drop Rasopathy", "Ras")

# Combine DataFrames
mut_v_copy_roc_df <- dplyr::bind_rows(raw_results[[1]],
                                      mutation_only_results[[1]],
                                      copy_only_results[[1]])
mut_v_copy_coef_df <- raw_results[[2]] %>%
  dplyr::full_join(mutation_only_results[[2]], by = 'feature',
                   suffix = c("_combined", "_mutation")) %>%
  dplyr::full_join(copy_only_results[[2]], by = 'feature',
                   suffix = c("", "_copy"))
mut_v_copy_auroc_df <- dplyr::bind_rows(raw_results[[3]],
                                        mutation_only_results[[3]],
                                        copy_only_results[[3]])

exp_v_cov_roc_df <- dplyr::bind_rows(raw_results[[1]],
                                     ras_expression_only_results[[1]],
                                     ras_covariate_only_results[[1]])
exp_v_cov_coef_df <- raw_results[[2]] %>%
  dplyr::full_join(ras_expression_only_results[[2]], by = 'feature',
                   suffix = c("_combined", "_exprs")) %>%
  dplyr::full_join(ras_covariate_only_results[[2]], by = 'feature',
                   suffix = c("", "_cov"))
exp_v_cov_auroc_df <- dplyr::bind_rows(raw_results[[3]],
                                       ras_expression_only_results[[3]],
                                       ras_covariate_only_results[[3]])

drop_genes_roc_df <- dplyr::bind_rows(raw_results[[1]],
                                      ras_nodrop_results[[1]],
                                      ras_droprasopathy_results[[1]])
drop_genes_coef_df <- raw_results[[2]] %>%
  dplyr::full_join(ras_nodrop_results[[2]], by = 'feature',
                   suffix = c("_combined", "_nodrop")) %>%
  dplyr::full_join(ras_droprasopathy_results[[2]], by = 'feature',
                   suffix = c("", "_rasopathy"))
drop_genes_auroc_df <- dplyr::bind_rows(raw_results[[3]],
                                        ras_nodrop_results[[3]],
                                        ras_droprasopathy_results[[3]])

# Reorder auroc_df
mut_v_copy_auroc_df <- mut_v_copy_auroc_df %>%
  dplyr::select(Features, dplyr::everything())
exp_v_cov_auroc_df <- exp_v_cov_auroc_df %>%
  dplyr::select(Features, dplyr::everything())
drop_genes_auroc_df <- drop_genes_auroc_df %>%
  dplyr::select(Features, dplyr::everything())

# Order Factors
mut_copy_levels <- c('Combined', "Mutation Only", "Copy Only")
exp_cov_levels <- c('Combined', "Exprs Only", "Cov Only")
drop_gene_levels <- c('Drop Ras', 'No Gene Drop', 'Drop Rasopathy')
train_type_levels <- c('train', 'test', 'shuffled')

mut_v_copy_roc_df$Features <- factor(mut_v_copy_roc_df$Features,
                                     levels = mut_copy_levels)
mut_v_copy_roc_df$train_type <- factor(mut_v_copy_roc_df$train_type,
                                     levels = train_type_levels)
mut_v_copy_auroc_df$Features <- factor(mut_v_copy_auroc_df$Features,
                                       levels = mut_copy_levels)

exp_v_cov_roc_df$Features <- factor(exp_v_cov_roc_df$Features,
                                     levels = exp_cov_levels)
exp_v_cov_roc_df$train_type <- factor(exp_v_cov_roc_df$train_type,
                                       levels = train_type_levels)
exp_v_cov_auroc_df$Features <- factor(exp_v_cov_auroc_df$Features,
                                       levels = exp_cov_levels)

drop_genes_roc_df$Features <- drop_genes_roc_df$Features %>%
  dplyr::recode('Combined' = 'Drop Ras')
drop_genes_auroc_df$Features <- drop_genes_auroc_df$Features %>%
  dplyr::recode('Combined' = 'Drop Ras')
drop_genes_roc_df$Features <- factor(drop_genes_roc_df$Features,
                                     levels = drop_gene_levels)
drop_genes_roc_df$train_type <- factor(drop_genes_roc_df$train_type,
                                       levels = train_type_levels)
drop_genes_auroc_df$Features <- factor(drop_genes_auroc_df$Features,
                                       levels = drop_gene_levels)

# Subset only to pancancer and also only cross validation performance
mut_v_copy_roc_df <- mut_v_copy_roc_df %>%
  dplyr::filter(disease == "PanCan") %>%
  dplyr::filter(train_type %in% train_type_levels)

exp_v_cov_roc_df <- exp_v_cov_roc_df %>%
  dplyr::filter(disease == "PanCan") %>%
  dplyr::filter(train_type %in% train_type_levels)

drop_genes_roc_df <- drop_genes_roc_df %>%
  dplyr::filter(disease == "PanCan") %>%
  dplyr::filter(train_type %in% train_type_levels)

# Plot ROC curves and feature bar plot
table_theme <- gridExtra::ttheme_default(base_size = 6,
                                         padding = unit(c(0.65, 0.65), "mm"))

plot_roc_info <- function(roc_df, auroc_df, custom_colors) {
  table_gg <- gridExtra::tableGrob(auroc_df, rows = NULL, theme = table_theme)
  
  feat_gg <- ggplot(roc_df, aes(x = fpr, y = tpr, color = Features)) +
    geom_step(aes(linetype = train_type), size = 0.3) +
    geom_segment(aes(x = 0 , y = 0, xend = 1, yend = 1),
                 linetype = "dashed", color = "black", size = 0.2) +
    annotation_custom(table_gg,
                      xmin = 0.6, xmax = 0.7, ymin = 0.08, ymax = 0.18) +
    scale_color_manual(values = custom_colors) +
    scale_linetype_manual('Data Type', labels = c('Train', 'Test', 'Shuffled'),
                          values = c('solid', 'dashed', 'dotted')) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    coord_fixed() + 
    theme_bw() +
    theme(axis.text = element_text(size = rel(0.5)),
          axis.title = element_text(size = rel(0.6)),
          axis.title.y = element_text(margin =
                                        margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(margin =
                                        margin(t = 3, r = 0, b = 0, l = 0)),
          legend.text = element_text(size = rel(0.6)),
          legend.title = element_text(size = rel(0.7)),
          legend.key = element_rect(size = 0.5),
          legend.position = 'right',
          legend.key.size = unit(0.7, 'lines'),
          legend.margin = margin(l = -0.3, unit = 'cm'))
  
  return(feat_gg)
}

# Set plotting colors
mut_v_copy_colors <- c("Combined" = "#984ea3", "Copy Only" = "#a65628",
                       "Mutation Only" = "#326F32")
mut_v_copy_fig <- plot_roc_info(mut_v_copy_roc_df, mut_v_copy_auroc_df,
                                mut_v_copy_colors)

base_fig_folder <- file.path(base_folder, 'RAS', 'figures')
mut_v_copy_fig_file <- file.path(base_fig_folder, 'mut_v_copy_fig.pdf')
ggplot2::ggsave(mut_v_copy_fig_file, plot = mut_v_copy_fig, dpi = 600,
                width = 3, height = 3)

exp_v_cov_colors <- c("Combined" = "#984ea3", "Exprs Only" = "#4daf4a",
                      "Cov Only" = "#ff7f00")
exp_v_cov_fig <- plot_roc_info(exp_v_cov_roc_df,
                               exp_v_cov_auroc_df,
                               exp_v_cov_colors)
exp_v_cov_fig_file <- file.path(base_fig_folder, 'exp_v_cov_fig.pdf')
ggplot2::ggsave(exp_v_cov_fig_file, plot = exp_v_cov_fig, dpi = 600,
                width = 3, height = 3)

drop_genes_colors <- c("Drop Ras" = "#984ea3", "No Gene Drop" = "#999999",
                       "Drop Rasopathy" = "#e41a1c")

drop_fig <- plot_roc_info(drop_genes_roc_df,
                          drop_genes_auroc_df,
                          drop_genes_colors)
drop_gene_fig_file <- file.path(base_fig_folder, 'drop_gene_fig.pdf')
ggplot2::ggsave(drop_gene_fig_file, plot = drop_fig, dpi = 600, width = 3,
                height = 3)

# Plot gene coefficients scatter
ggplot(mut_v_copy_coef_df, aes(x = weight_mutation, y = weight,
                               color = weight_combined)) +
  geom_point(alpha = 0.8, size = 0.1) +
  scale_color_gradient2('Combined', low = "blue", mid = "grey", high = "red") +
  xlab("Mutation Only - Gene Weight") +
  ylab("Copy Number Only - Gene Weight") +
  geom_text_repel(data = subset(mut_v_copy_coef_df,
                                (weight > 0.075 | weight < -0.015) |
                                  (weight_mutation > 0.1 |
                                     weight_mutation < -0.0625)),
                  arrow = arrow(length = unit(0.02, 'npc')),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  box.padding = 0.17,
                  point.padding = 0.1,
                  size = 1.8,
                  fontface = 'italic',
                  aes(x = weight_mutation, y = weight, label = feature)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.title = element_text(size = rel(0.6)),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 3, r = 0, b = 0, l = 0)),
        legend.text = element_text(size = rel(0.4)),
        legend.title = element_text(size = rel(0.5)),
        legend.key = element_rect(size = 0.2),
        legend.position = 'right',
        legend.key.size = unit(0.4, 'lines'),
        legend.margin = margin(l = -0.3, unit = 'cm'))

mut_v_copy_weights_file <- file.path(base_fig_folder, 'mut_v_copy_weights.pdf')
ggplot2::ggsave(mut_v_copy_weights_file, dpi = 600, width = 3, height = 2)

# Set dropped genes that were removed to zero weight
drop_genes_coef_df$weight[is.na(drop_genes_coef_df$weight)] <- 0
drop_genes_coef_df$weight_combined[is.na(drop_genes_coef_df$weight_combined)] <- 0
drop_genes_coef_df$weight_nodrop[is.na(drop_genes_coef_df$weight_nodrop)] <- 0

ras_genes <- c('KRAS', 'HRAS', 'NRAS')

ggplot(drop_genes_coef_df, aes(x = weight_nodrop, y = weight,
                               color = weight_combined)) +
  geom_point(alpha = 0.8, size = 0.1) +
  scale_color_gradient2('Drop Ras Genes', low = "blue", mid = "lightgrey",
                        high = "red") +
  xlab("No Gene Drop - Gene Weight") +
  ylab("Drop Rasopathy - Gene Weight") +
  coord_fixed() + 
  geom_text_repel(data = subset(drop_genes_coef_df,
                                (weight_combined > 0.09 |
                                   weight_combined < -0.06) |
                                  feature %in% ras_genes),
                  arrow = arrow(length = unit(0.02, 'npc')),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  size = 1.8,
                  fontface = 'italic',
                  aes(x = weight_nodrop, y = weight, label = feature)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.title = element_text(size = rel(0.6)),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 3, r = 0, b = 0, l = 0)),
        legend.text = element_text(size = rel(0.4)),
        legend.title = element_text(size = rel(0.5)),
        legend.key = element_rect(size = 0.2),
        legend.position = 'right',
        legend.key.size = unit(0.4, 'lines'),
        legend.margin = margin(l = -0.3, unit = 'cm'))

drop_rasopathy_weights_file <- file.path(base_fig_folder,
                                         'drop_rasopathy_weights.pdf')
ggplot2::ggsave(drop_rasopathy_weights_file, dpi = 600, width = 3, height = 2)
