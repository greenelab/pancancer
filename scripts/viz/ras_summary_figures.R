# Gregory Way 2017
# PanCancer Classifier
# scripts/viz/ras_summary_figures.R
#
# Visualize summary for Ras Classifier Scores
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/viz/ras_summary_figures.R
#
# Output:
# Several figures to summarize Ras findings

checkpoint::checkpoint("2017-06-01", checkpointLocation = ".")

library(dplyr)
library(pheatmap)
library(ggplot2)
library(readr)
library(cowplot)
library(gridExtra)
library(Hmisc)
source(file.path("scripts", "util", "pancancer_util.R"))

results_folder <- file.path("classifiers", "RAS")
results <- parse_summary(file.path(results_folder, "classifier_summary.txt"))

# 1) Heatmap of the distribution of aberrant events across tumors
heatmap_plot_file <- file.path(results_folder, "figures", "ras_heatmap.pdf")
ras_heatmap_file <- file.path(results_folder, "figures", "all_ras_heatmap.pdf")
heat_file <- file.path(results_folder, "summary_counts.csv")
heat_df <- readr::read_csv(heat_file)

heat_comb_df <- heat_df %>% dplyr::mutate(RAS = HRAS_y + NRAS_y + KRAS_y) %>%
  dplyr::mutate(RAS_gain = HRAS_gain_y + NRAS_gain_y + KRAS_gain_y)

prop_matrix <- as.matrix(heat_comb_df[, c("RAS_gain", "RAS")])
rownames(prop_matrix) <- heat_comb_df$DISEASE
colnames(prop_matrix) <- c("Gain", "Mutation")

# All diseases that are used in building the classifier
ras_dis <- results[["Diseases"]]

# Build a vector for heatmap labels
classifier <- c()
for (disease in rownames(prop_matrix)) {
  if (disease %in% ras_dis) {
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

# Plot heatmap without collapsing Ras genes
heat_ras_df <- heat_df %>% dplyr::select(c("NRAS_gain_y", "HRAS_gain_y",
                                           "KRAS_gain_y", "NRAS_y", "HRAS_y",
                                           "KRAS_y"))
colnames(heat_ras_df) <- c("NRAS Gain", "HRAS Gain", "KRAS Gain",
                           "NRAS", "HRAS", "KRAS")
heat_ras_df <- as.data.frame(heat_ras_df)
rownames(heat_ras_df) <- heat_comb_df$DISEASE
heat_ras_df <- heat_ras_df[order(heat_ras_df$KRAS, decreasing = TRUE), ]

# Plot and save heatmap
pheatmap(t(heat_ras_df * 100), scale = "none", cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.0f", fontsize_number = 8,
         number_color = "black", annotation_col = classifier,
         annotation_names_col = FALSE, legend = FALSE,
         filename = ras_heatmap_file,
         width = 8, height = 2)

# 2) Coefficients contributing to the model
coef_plot_file <- file.path(results_folder, "figures", "ras_coef_plot.svg")
coef_df <- results[["Coefficients"]]
coef_df <- coef_df[, -1]
coef_df <- coef_df[order(coef_df$weight, decreasing = FALSE), ]
coef_df$rank <- 1:nrow(coef_df)

p <- ggplot(coef_df, aes(x = 1:nrow(coef_df), y = weight)) +
  geom_point(fill = "black", size = 0.01) +
  base_theme + theme(axis.line.x = element_line(),
                     axis.line.y = element_line(),
                     axis.ticks = element_line(),
                     axis.title = element_text(size = rel(1.5)),
                     plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "cm")) +
  labs(list(x = "Rank", y = "Weight")) +
  scale_y_continuous(breaks = seq(-0.25, 0.25, 0.05)) +
  scale_x_continuous(breaks = seq(0, 8000, 2000)) +
  geom_segment(aes(x = 0, y = 0, yend = 0, xend = nrow(coef_df)),
               colour = "red", linetype = "dashed", size = 0.2)

p <- add_arrow_label(p = p, x = 1050, y = -0.205, label = "CDK13",
                     offset = c(0, -0.003, -685, -.0002))
p <- add_arrow_label(p = p, x = 1150, y = -0.190, label = "PDLIM4",
                     offset = c(40, -0.003, -450, 0.006))
p <- add_arrow_label(p = p, x = 1480, y = -0.172, label = "PDE5A",
                     offset = c(80, -0.001, -690, .0003))
p <- add_arrow_label(p = p, x = 1650, y = -0.155, label = "PURB",
                     offset = c(80, -0.001, -680, -.00013))
p <- add_arrow_label(p = p, x = 1700, y = -0.14, label = "FADS3",
                     offset = c(80, -0.001, -710, .0002))
p <- add_arrow_label(p = p, x = 2000, y = -0.12, label = "SEPP1",
                offset = c(90, -0.0001, -700, -.00015))
p <- add_arrow_label(p = p, x = 1800, y = -0.1, label = "ALDOC",
                     offset = c(80, 0, -750, -.0002))
p <- add_arrow_label(p = p, x = 2100, y = -0.085, label = "PAPLN",
                     offset = c(80, 0, -800, -.0002))
p <- add_arrow_label(p = p, x = 1800, y = -0.07, label = "CLU",
                     offset = c(80, 0, -550, -.0006))
p <- add_arrow_label(p = p, x = 1900, y = -0.055, label = "CUL1",
                     offset = c(80, 0, -690, -.0002))

p <- add_arrow_label(p = p, x = 6800, y = 0.15, label = "PBX3",
                     offset = c(-80, .0004, 600, -.004))
p <- add_arrow_label(p = p, x = 6500, y = 0.12, label = "SPRY2",
                     offset = c(-80, 0, 760, -.0006))
p <- add_arrow_label(p = p, x = 6200, y = 0.1, label = "PPP1R3B",
                     offset = c(-80, .0004, 1000, -.001))
p <- add_arrow_label(p = p, x = 6050, y = 0.084, label = "C15orf52",
                     offset = c(-80, .0004, 950, -.0015))
p <- add_arrow_label(p = p, x = 5850, y = 0.062, label = "MLPH",
                     offset = c(-80, .0004, 700, -.0015))
p <- add_arrow_label(p = p, x = 5350, y = 0.041, label = "ERRFI1",
                     offset = c(-80, .0004, 790, -.0015))
p <- add_arrow_label(p = p, x = 6750, y = 0.032, label = "CMAS",
                     offset = c(-80, -.001, 270, .01))

p <- add_arrow_label(p = p, x = 6500, y = 0.016, label = "log10_mut",
                     offset = c(-80, -0.0142, 900, 0.005))

svg(coef_plot_file, width = 2.5, height = 2.25)
p
dev.off()

# 3) Plot distributions of predictions according to variant classification
var_plot_file <- file.path(results_folder, "figures", "variant_fill_map.svg")
mut_df <- readr::read_tsv(file.path(results_folder, "tables",
                                    "mutation_classification_scores.tsv"))

consider_mutations <- c("3'UTR", "5'UTR", "Intron", "Frame_Shift_Del",
                        "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                        "Missense_Mutation", "Nonsense_Mutation",
                        "Nonstop_Mutation", "RNA", "Splice_Site")

# Generate RAS Gain Variable
mut_df <- mut_df %>% mutate(RAS_gain = max(c(KRAS_gain, NRAS_gain, HRAS_gain)))

silent_df <- mut_df %>% filter(Variant_Classification == "Silent") %>%
  filter(total_status == 0)
delet_df <- mut_df %>% filter(Variant_Classification %in% consider_mutations)

mut_filtered_df <- dplyr::bind_rows(delet_df, silent_df)

# Separate classes of mutations to summarize
copy_num_df <- mut_df %>% filter(RAS_gain == 1) %>%
  filter(TP53 == 0) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Loss")
missense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Missense_Mutation") %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Missense")
nonsense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Nonsense_Mutation") %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Nonsense")
indel_df <- mut_filtered_df %>% filter(Variant_Classification %in%
                                       c("Frame_Shift_Del", "Frame_Shift_Ins",
                                         "In_Frame_Del", "In_Frame_Ins")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Indel")
utr_df <- mut_filtered_df %>%
  filter(Variant_Classification %in% c("3'UTR", "5'UTR", "Intron")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "UTR")
silent_df <- silent_df %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Silent")
splice_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Splice_Site") %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Splice")
wt_df <- mut_df %>% subset(total_status == 0) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "WT")
hyper_df <- mut_df %>%
  filter(hypermutated == 1) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Hyper")

final_df <- dplyr::bind_rows(list(missense_df, nonsense_df, indel_df, utr_df,
                                  splice_df, silent_df, copy_num_df, wt_df,
                                  hyper_df))

colnames(final_df) <- c("ID", "Gene", "Disease", "Weight", "HGVSc", "HGVSp",
                        "Class")

# Plot summary distribution of variant classes prediction scores
ggplot(final_df, aes(Weight, ..count.., fill = Class)) +
  geom_density(position = "fill", size = 0.1) +
  geom_segment(aes(x = 0.5, y = 0, yend = 1, xend = 0.5), colour = "black",
               linetype = "dashed", size = 0.4) +
  labs(list(x = "Probability", y = "Proportion")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + base_theme +
  theme(legend.position = c(1.1, 0.65),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 7),
        plot.margin = unit(c(0.2, 1.5, 0, 0.1),"cm"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 12))
ggsave(var_plot_file, width = 4, height = 3.8)
dev.off()

# 4) Show mutation frequencies and scores
mut_weight_df <- mut_filtered_df %>% filter(!is.na(weight))
mut_weight_df <- mut_weight_df[mut_weight_df$hypermutated != 1, ]

aa_df <- mut_weight_df %>%
  group_by(HGVSp, Variant_Classification, Hugo_Symbol) %>%
  summarise(Mean = mean(weight, na.rm = TRUE),
            SD = sd(weight, na.rm = TRUE),
            count = n(),
            low_CI = get_boot(weight),
            high_CI = get_boot(weight, low = FALSE))

nuc_df <- mut_weight_df %>%
  group_by(HGVSc, Variant_Classification, Hugo_Symbol) %>%
  summarise(Mean = mean(weight),
            SD = sd(weight, na.rm = TRUE),
            count = n(),
            low_CI = get_boot(weight),
            high_CI = get_boot(weight, low = FALSE))

aa_df <- aa_df[order(aa_df$count, decreasing = TRUE),]
nuc_df <- nuc_df[order(nuc_df$count, decreasing = TRUE),]
write.table(aa_df, file = file.path(results_folder, "tables",
                                    "amino_acid_mutation_scores.tsv"),
            sep = "\t", row.names = FALSE)
write.table(nuc_df, file = file.path(results_folder, "tables",
                                     "nucleotide_mutation_scores.tsv"),
            sep = "\t", row.names = FALSE)

# Plot summary distribution of variant classes prediction scores
braf_df <- final_df[complete.cases(final_df), ]
braf_df <- braf_df[braf_df$HGVSp == "p.Val600Glu", ]

braf_df$Disease <- dplyr::recode(braf_df$Disease,
                                 "BLCA" = "Other", "CHOL" = "Other",
                                 "GBM" = "Other", "HNSC" = "Other",
                                 "KIRP" = "Other", "LGG" = "Other",
                                 "READ" = "Other")

braf_plot_file <- file.path(results_folder, "figures",
                            "brafv600e_distribution.svg")
braf_plot <- ggplot(braf_df, aes(Weight, fill = Disease)) +
  geom_density(alpha = 0.4) + theme_bw() +
  ylab("Density") + xlab("BRAFV600E Classifier Score")

svg(braf_plot_file, width = 4, height = 3)
braf_plot
dev.off()

# 5) RAS Summary Counts Distribution
ras_count_file <- file.path(results_folder, "tables",
                            "ras_events_per_sample.tsv")
ras_summary_count_df <- readr::read_tsv(ras_count_file,
                                        col_types = cols(.default = "c",
                                                         "weight" = "d",
                                                         "total_status" = "c"))
ras_summary_count_df$copy_count <- factor(ras_summary_count_df$copy_count,
                                          levels = c("0", "1", "2", "3","4", 
                                                     "5", "6", "7", "8", "9",
                                                     "10"))
ras_summary_count_df$copy_count <-
  dplyr::recode(ras_summary_count_df$copy_count, "6" = ">6", "7" = ">6",
                "8" = ">6", "9" = ">6", "10" = ">6")

# Get summary statistics for each comparison
mut_ras_prop <- ras_summary_count_df %>% group_by(mutation_count) %>%
  summarize(mean_ras = round(mean(as.numeric(total_status)), 2))
cop_ras_prop <- ras_summary_count_df %>% group_by(copy_count) %>%
  summarize(mean_ras = round(mean(as.numeric(total_status)), 2))

mut_ras_count <- ras_summary_count_df %>% group_by(mutation_count) %>% tally()
cop_ras_count <- ras_summary_count_df %>% group_by(copy_count) %>% tally()

# Combine to get summary tables
mut_sum <- dplyr::inner_join(mut_ras_count, mut_ras_prop, by = "mutation_count")
cop_sum <- dplyr::inner_join(cop_ras_count, cop_ras_prop, by = "copy_count")

med_weight <- median(ras_summary_count_df$weight)

classifier_count_theme <- base_theme +
  theme(legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(0.8)),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks = element_line(),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        plot.margin = unit(c(0.2, 2.5, 0.2, 0.2), "cm")) 
  
mut <- ggplot(ras_summary_count_df, aes(x = mutation_count, y = weight)) +
  geom_boxplot(aes(fill = total_status)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(name = "RAS Status", values = c("#3B9AB2", "#F2300F"),
                    labels = c("0" = "Wild-Type", "1" = "Activated")) +
  geom_text(data = mut_sum, aes(x = mutation_count, y = 1.06,
                                label = paste0(n, "\n", mean_ras))) +
  classifier_count_theme +
  labs(list(x = "Number of Other Ras Pathway Mutations",
            y = "RAS Classifier Score"))

cop <- ggplot(ras_summary_count_df, aes(x = copy_count, y = weight)) +
  geom_boxplot(aes(fill = total_status)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(name = "RAS Status", values = c("#3B9AB2", "#F2300F"),
                    labels = c("0" = "Wild-Type", "1" = "Activated")) +
  geom_text(data = cop_sum, aes(x = copy_count, y = 1.06,
                                label = paste0(n, "\n", mean_ras))) +
  classifier_count_theme +
  labs(list(x = "Number of Other Ras Pathway Copy Number Events",
            y = "RAS Classifier Score"))

ras_counts_fig <- file.path(results_folder, "figures", "ras_events_counts.svg")
svg(ras_counts_fig, width = 6, height = 8.6)
plot_grid(mut , cop, align = "v", nrow = 2)
dev.off()

# 6) Performance Metrics Distribution across pathway members
perf_metric_file <- file.path(results_folder, "tables",
                              "all_gene_metric_ranks.tsv")
metric_ranks <- readr::read_tsv(perf_metric_file,
                                col_types = cols(.default = "c",
                                                 "AUROC" = "d",
                                                 "AUPRC" = "d",
                                                 "AUROC Rank" = "i",
                                                 "AUPRC Rank" = "i"))

auprc_violin <- ggplot(metric_ranks, aes(y = AUPRC, x = paste(ras),
                                        fill = paste(ras))) +
  geom_violin() +
  theme(legend.position = "none") +
  xlab("") +
  scale_x_discrete(labels = c("0" = "Other", "1" = "Ras Pathway Genes"))

auroc_violin <- ggplot(metric_ranks, aes(y = AUROC, x = paste(ras),
                                         fill = paste(ras))) +
  geom_violin() +
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  xlab("") +
  scale_x_discrete(labels = c("0" = "Other", "1" = "Ras Pathway Genes"))

auprc_plot <- ggplot(metric_ranks, aes(x = `AUPRC Rank`, y = AUPRC)) +
  geom_point(color = "darkgrey") +
  geom_point(data = metric_ranks[metric_ranks$ras == 1, ], color = "red")

auroc_plot <- ggplot(metric_ranks, aes(x = `AUROC Rank`, y = AUROC)) +
  geom_point(color = "darkgrey") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_point(data = metric_ranks[metric_ranks$ras == 1, ], color = "red")

# Get the top genes by both metrics
top_auprc_genes <- metric_ranks[order(metric_ranks$`AUPRC Rank`), 1:2]
top_auprc_table_grob <- tableGrob(top_auprc_genes[1:20, ])
auprc_plot <- auprc_plot +
  annotation_custom(top_auprc_table_grob, xmin = 10000,
                    xmax = 15000, ymin = 0.1, ymax = 0.45)

top_auroc_genes <- metric_ranks[order(metric_ranks$`AUROC Rank`), c(1, 5)]
top_auroc_table_grob <- tableGrob(top_auroc_genes[1:10, ])
auroc_plot <- auroc_plot +
  annotation_custom(top_auroc_table_grob, xmin = 10000,
                    xmax = 15000, ymin = 0.6, ymax = 0.95)

auprc_distribution_fig <- file.path(results_folder, "figures", 
                                    "auprc_distribution.svg")

svg(auprc_distribution_fig, width = 11.5, height = 7.5)
plot_grid(auprc_plot, auprc_violin, align = "h", ncol = 2)
dev.off()

auroc_distribution_fig <- file.path(results_folder, "figures", 
                                    "auroc_distribution.svg")

svg(auroc_distribution_fig, width = 11, height = 7.5)
plot_grid(auroc_plot, auroc_violin, align = "h", ncol = 2)
dev.off()
