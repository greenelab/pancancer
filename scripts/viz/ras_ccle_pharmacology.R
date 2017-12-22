# Gregory Way 2017
# PanCancer Classifier
# scripts/viz/ras_ccle_pharmacology.R
#
# Visualize correlations between Ras classifier score and two MEK inhibitor
# activity (AZD6244 (Selumetinib) and PD-0325901)
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/viz/ras_ccle_pharmacology.R
#
# Output:
# Overall and tissue-specific correlation plots between drug activity and
# Ras classifier scores

library(dplyr)
library(ggplot2)
library(ggpmisc)

# Input Files
pharm_file <- file.path("data", "pharmacology_predictions_ccle.tsv")
pharm_full_df <- readr::read_tsv(pharm_file)

# Output Files
base_path <- file.path("figures", "cell_line")
selum_file <- file.path(base_path, "selumetinib.pdf")
selum_braf_file <- file.path(base_path, "selumetinib_braf.pdf")
selum_tissue_file <- file.path(base_path, "selumetinib_tissues.pdf")
selum_tissue_braf_file <- file.path(base_path, "selumetinib_tissues_braf.pdf")

pd_file <- file.path(base_path, "pd_0325901.pdf")
pd_braf_file <- file.path(base_path, "pd_0325901_braf.pdf")
pd_tissue_file <- file.path(base_path, "pd_0325901_tissues.pdf")
pd_tissue_braf_file <- file.path(base_path, "pd_0325901_tissues_braf.pdf")

# Define Plotting Function
plot_drug <- function(pharm_df, compound, tissues = NULL, include_braf = FALSE, 
                      facet_tissue = TRUE, se = FALSE) {
  # Output scatter plots with correlations, visualizing drug activity
  # compared to Ras classifier Score
  #
  # Arguments:
  # pharm_df - dataframe of compound activity by cell line with Ras/BRAF status
  # compound - a specific compound to visualize
  # tissues - a list of tissues to consider plotting specifically in facets
  # include_braf - boolean to include BRAF in considering mutation status
  # facet_tissue - boolean of tissues to determine to plot in facet_wrap
  # se - boolean to plot standard error intervals in geom_smooth
  #
  # Output:
  # Scatter plot with correlation information

  pharm_subset_df <- pharm_df[pharm_df$Compound == compound, ]
  if (!is.null(tissues)) {
    pharm_subset_df <- pharm_subset_df %>%
      dplyr::filter(tissue %in% focus_tissues)
  }
  if (include_braf) {
    pharm_subset_df$ras_status[pharm_subset_df$BRAF_MUT == 1] <- 1
    legend_label <- "Ras/BRAF Status"
  } else {
    legend_label <- "Ras Status"
  }
  
  if (compound == "AZD6244") {
    compound <- "Selumetinib"
  }

  formula <- y ~ x

  p <- ggplot(pharm_subset_df, aes(x = weight, y = ActArea,
                                   color = as.factor(ras_status),
                                   fill = as.factor(ras_status))) +
    geom_point(alpha = 0.5, size = 2) +
    scale_x_continuous(breaks = c(0, 0.5, 1),
                       limits = c(-0.1, 1.1)) +
    geom_smooth(method = "lm", se = se) +
    geom_segment(aes(x = 0.5, y = -0.1, xend = 0.5, yend = 6),
                 linetype = "dashed", color = "grey") +
    scale_fill_manual(values = c("#377eb8", "#ff7f00"),
                      name = legend_label,
                      breaks = c(0, 1),
                      labels = c("Wild-Type", "Mutant")) +
    scale_color_manual(values = c("#377eb8", "#ff7f00"),
                       name = legend_label,
                       breaks = c(0, 1),
                       labels = c("Wild-Type", "Mutant")) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 label.x.npc = 0.17, label.y.npc = 0.92,
                 formula = formula,
                 parse = TRUE, size = 4, na.rm = TRUE,
                 rr.digits = 1) +
    stat_fit_glance(method = "lm", geom = "text",
                    label.x.npc = 0.8, label.y.npc = 0.97,
                    method.args = list(formula = formula), size = 4,
                    aes(label = paste("P = ",
                                      signif(..p.value.., digits = 1),
                                      sep = ""))) +
    xlab("Ras Classifier Score") +
    ylab("Activity Area") +
    ggtitle(compound, subtitle = "CCLE Response") + 
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  if (facet_tissue) {
    p <- p + facet_wrap("tissue")
  }
  
  return(p)
}

# Set constants - tissues are set based on sample size and class size balance
focus_tissues <- c("LUNG", "SKIN", "BREAST", "OVARY", "LARGE INTESTINE",
                   "STOMACH")

# Plot and output figures
# Selumetinib all tissues - Ras only
sel_fig <- plot_drug(pharm_full_df, "AZD6244", facet_tissue = FALSE, se = TRUE)
ggplot2::ggsave(selum_file, plot = sel_fig, dpi = 600, height = 4, width = 5.25)

# Selumetinib all tissues - Ras and BRAF
sel_braf_fig <- plot_drug(pharm_full_df, "AZD6244", facet_tissue = FALSE,
                          include_braf = TRUE, se = TRUE)
ggplot2::ggsave(selum_braf_file, plot = sel_braf_fig, dpi = 600, height = 4,
                width = 5.25)

# PD-0325901 all tissues - Ras only
p_fig <- plot_drug(pharm_full_df, "PD-0325901", facet_tissue = FALSE, se = TRUE)
ggplot2::ggsave(pd_file, plot = p_fig, dpi = 600, height = 4,
                width = 5.25)

# PD-0325901 all tissues - Ras and BRAF
p_braf_fig <- plot_drug(pharm_full_df, "PD-0325901", facet_tissue = FALSE,
                        include_braf = TRUE, se = TRUE)
ggplot2::ggsave(pd_braf_file, plot = p_braf_fig, dpi = 600,
                height = 4.5, width = 8)

# Selumetinib stratified by tissue - Ras only
sel_tis_fig <- plot_drug(pharm_full_df, "AZD6244", tissues = focus_tissues)
ggplot2::ggsave(selum_tissue_file, plot = sel_tis_fig, dpi = 600,
                height = 4.5, width = 8)

# Selumetinib stratified by tissue - Ras and BRAF
sel_tis_braf_fig <- plot_drug(pharm_full_df, "AZD6244",
                              tissues = focus_tissues, include_braf = TRUE)
ggplot2::ggsave(selum_tissue_braf_file, plot = sel_tis_braf_fig, dpi = 600,
                height = 4.5, width = 8)

# PD-0325901 stratified by tissue - Ras only
p_tis_fig <- plot_drug(pharm_full_df, "PD-0325901", tissues = focus_tissues)
ggplot2::ggsave(pd_tissue_file, plot = p_tis_fig, dpi = 600,
                height = 4.5, width = 8)

# PD-0325901 stratified by tissue - Ras and BRAF
p_tis_braf_fig <- plot_drug(pharm_full_df, "PD-0325901",
                              tissues = focus_tissues, include_braf = TRUE)
ggplot2::ggsave(pd_tissue_braf_file, plot = p_tis_braf_fig,
                dpi = 600, height = 4.5, width = 8)
