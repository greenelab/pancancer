# Gregory Way 2017
# PanCancer Classifier
# scripts/viz/ras_differential_expression_figure.R
#
# Visualize differences between differential expression analysis of Ras
# wild-type to Ras mutant tumors and the machine learning derived gene
# importance scores
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/viz/ras_differential_expression_figure.R
#
# Output:
# 1) Differentially expressed genes vs. classifier coefficients (Figure S3)
# 2) Volcano plot of differentially expressed genes (Figure S1)
# 3) All differentially expressed genes resource (Data S1)

library(dplyr)
library(limma)
library(ggplot2)
library(ggpmisc)
library(ggrepel)

set.seed(123)

# Load Data
rnaseq_file <- file.path("data", "RNAseq_scaled_all_genes.tsv")
ras_status_file <- file.path("data", "Ras_sample_status.tsv")
coef_file <- file.path("classifiers", "RAS", "classifier_coefficients.tsv")

rnaseq_df <- readr::read_tsv(rnaseq_file)
ras_status_df <- readr::read_tsv(ras_status_file)
coef_df <- readr::read_tsv(coef_file)

mad_genes <- intersect(coef_df$feature, colnames(rnaseq_df))
gene_coef_df <- coef_df %>% dplyr::filter(feature %in% mad_genes)
nonzero_coef_df <- coef_df %>% dplyr::filter(abs > 0)

# Process design matrix to adjust for cancer-type
ras_design <- model.matrix(~ 0 + ras_status_df$total_status +
                             ras_status_df$DISEASE)
colnames(ras_design) <- c("ras_status", sort(unique(ras_status_df$DISEASE)))

# Determine differentially expressed genes
fit <- lmFit(t(rnaseq_df[, 2:ncol(rnaseq_df)]), ras_design)
ras_contrast <- makeContrasts(ras_status, levels = ras_design)
fit_contrasts <- contrasts.fit(fit, ras_contrast)
fit2 <- eBayes(fit_contrasts)

# Begin processing DEGs
deg_df <- data.frame(cbind(fit2$coefficients, fit2$p.value,
                           -log10(fit2$p.value)))
colnames(deg_df) <- c("stat", "p", "logp")
deg_df$gene <- rownames(deg_df)

data_s1 <- dplyr::full_join(deg_df, gene_coef_df, by = c("gene" = "feature"))
data_s1 <- data_s1[order(data_s1$logp, decreasing = TRUE), ]
data_s1_file <- file.path("classifiers" ,"RAS", "tables",
                          "differentially_expressed_genes.tsv")

# Output DEG analysis results
write.table(data_s1, data_s1_file, row.names = FALSE, sep = "\t")

# Output Plots
# 1) Differentially expressed genes vs. classifier coefficients
color_logic <- (data_s1$weight > 0.06 | data_s1$weight < -0.05)

ggplot(data_s1, aes(x = weight, y = stat)) +
  geom_point(alpha = 0.5, color = ifelse(color_logic, "red", "grey50")) +
  xlab("Ras Classifier Weight") +
  geom_text_repel(data = subset(data_s1, color_logic),
                  arrow = arrow(length = unit(0.01, "npc")),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  size = 1.4,
                  fontface = "italic",
                  point.padding = 0.1,
            aes(x = weight, y = stat, label = gene)) +
  ylab("Diff Exprs Fold Change") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(0.65)),
        axis.title = element_text(size = rel(0.8)),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 3, r = 0, b = 0, l = 0)),
        plot.margin = margin(r = 0.3, l = 0.1, unit = "cm"))

diff_exp_fig <- file.path("classifiers", "RAS", "figures", "diff_exprs.pdf")
ggplot2::ggsave(diff_exp_fig, dpi = 600, width = 2.3, height = 2)

# 2) Differentially expressed genes volcano plot
volcano_color_logic <- (data_s1$logp > 20) | 
  (data_s1$logp > 13 & data_s1$stat < 0) |
  (data_s1$gene %in% c("KRAS", "HRAS", "NRAS"))

ggplot(data_s1, aes(x = stat, y = logp)) +
  geom_point(alpha = 0.5, color = ifelse(volcano_color_logic, "red", "grey50")) +
  xlab("Fold Change") +
  geom_text_repel(data = subset(data_s1, volcano_color_logic),
                  arrow = arrow(length = unit(0.01, "npc")),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  size = 1.75,
                  fontface = "italic",
                  aes(x = stat, y = logp, label = gene)) +
  ylab("-log10 p value") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(0.8)),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 3, r = 0, b = 0, l = 0)),
        plot.margin = margin(r = 0.3, l = 0.1, unit = "cm"))

volcano_fig <- file.path("classifiers", "RAS", "figures",
                         "diff_exprs_volcano.pdf")
ggplot2::ggsave(volcano_fig, dpi = 600, width = 5, height = 4.5)
