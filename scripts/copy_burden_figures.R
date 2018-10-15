# Gregory Way 2017
# PanCancer Classifier
# scripts/copy_burden_figures.R
#
# Generate figures for visualizing copy burden across different samples
# stratified by TP53 mutation status
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/copy_burden_figures.R
#
# Output:
# Two figures summarizing copy burden across TCGA Pan Can samples

library(ggplot2)

# Set File Names
base_file <- file.path("classifiers", "TP53")
burden_file <- file.path(base_file, "tables", "copy_burden_predictions.tsv")
snaptron_file <- file.path("scripts", "snaptron",
                           "junctions_with_mutations.csv.gz")
frac_alt_plot <- file.path(base_file, "figures", "fraction_altered_plot.pdf")
violin_plot <- file.path(base_file, "figures", "seg_altered_violin_plot.pdf")

# Load Files
copy_burden <- readr::read_tsv(burden_file)
junc_df <- readr::read_csv(snaptron_file)
junc_df <- junc_df[,-1]
junc_df <- junc_df[!duplicated(junc_df), ]

# Location of the silent mutation and truncation
junc_exon_df = junc_df[junc_df$start == "7675237", ]
silent_junc <- junc_exon_df[junc_exon_df$Variant_Classification == "Silent", ]
silent_junc <- silent_junc[silent_junc$snaptron_id == "13945701", ]
silent_junc <- silent_junc[silent_junc$TP53 %in% 0, ]
silent_junc <- silent_junc[silent_junc$include %in% 1, ]

ggplot(copy_burden, aes(weight, frac_altered, color = factor(TP53))) +
  geom_point(alpha = 0.6, size = 0.3) + theme_bw() +
  xlab("TP53 Inactivation Probability") +
  ylab("CNV Burden (Fraction Altered)") +
  labs(color = "TP53 Status")
ggsave(frac_alt_plot, width = 5, height = 4)

# Build and Process Copy Burden DataFrame
copy_burden$silent <- 0
copy_burden[copy_burden$Sample %in% silent_junc$tcga_id, "silent"] <- 1
silent_and_junc <-  copy_burden[copy_burden$silent == 1, ]
silent_and_junc$TP53 <- "c.375G>T Mutation"

copy_burden[copy_burden$total_status == 0, "TP53"] = "Wild-Type"
copy_burden[copy_burden$total_status == 1, "TP53"] = "TP53 Loss of Function"

plot_ready <- copy_burden[, c("frac_altered", "TP53")]
plot_ready <- rbind(plot_ready, silent_and_junc[, c("frac_altered", "TP53")])

false_negatives <- copy_burden[(copy_burden$total_status == 1) &
                               (copy_burden$weight < 0.5), ]
false_negatives$TP53 <- "False Negative"
plot_ready <- rbind(plot_ready, false_negatives[, c("frac_altered", "TP53")])

false_positives <- copy_burden[(copy_burden$total_status == 0) &
                               (copy_burden$weight >= 0.5), ]
false_positives$TP53 <- "False Positive"
plot_ready <- rbind(plot_ready, false_positives[, c("frac_altered", "TP53")])

predicted_neg <- copy_burden[copy_burden$weight < 0.5, ]
predicted_neg$TP53 <- "Predicted Wild-Type"
plot_ready <- rbind(plot_ready, predicted_neg[, c("frac_altered", "TP53")])

predicted_pos <- copy_burden[copy_burden$weight >= 0.5, ]
predicted_pos$TP53 <- "Predicted Loss"
plot_ready <- rbind(plot_ready, predicted_pos[, c("frac_altered", "TP53")])

plot_levels <- c("c.375G>T Mutation", "False Positive",
                 "Predicted Loss", "TP53 Loss of Function",
                 "False Negative", "Predicted Wild-Type", "Wild-Type")

plot_ready$TP53 <- factor(plot_ready$TP53, levels = plot_levels)

# Build violin plots for copy number alterations comparison
ggplot(plot_ready, aes(x = TP53, y = frac_altered)) +
  ylab("CNV Burden (Fraction Altered)") + xlab("TP53 Status") +
  labs(fill = "") + geom_violin(aes(fill = TP53), size = 0.3, alpha = 0.3,
                                adjust = 0.7, trim = TRUE) +
  geom_boxplot(aes(fill = TP53), size = 0.3, width = 0.1, outlier.size = 0.3) +
  coord_flip() + geom_hline(yintercept = 0.5, linetype = "dashed",
                            color = "red") +
  theme(legend.position = c(1.4, 0.7), axis.text.y = element_blank(),
        axis.text.x = element_text(size = rel(0.7)),
        axis.title = element_text(size = rel(0.7)),
        legend.text = element_text(size = rel(0.45)),
        legend.key = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.background = element_rect(fill = alpha("white", 0)),
        panel.grid.major = element_line(color = "white", size = 0.3),
        panel.grid.minor = element_line(color = "white", size = 0.3),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1, 2.6, 0.2, 0.2),"cm"),
        panel.border = element_rect(fill = NA, size = 0.4)) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1), color = FALSE)

ggsave(violin_plot, height = 2.25, width = 2.5)
