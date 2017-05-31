# Gregory Way 2017
# PanCancer Classifier
# scripts/ddr_summary_figures.R
#
# Visualize summary for DNA Damage Repair/Response (DDR) classification scores
#
# Usage: Run by assigning where the within classifier summary is and where the
#        Pan Cancer classifier summary is
#
#     Rscript --vanilla ddr_summary_figures.R
#
# Output:
# Several figures to summarize DDR findings

library(dplyr)
library(pheatmap)
library(ggplot2)
source(file.path("scripts", "util", "pancancer_util.R"))

results_folder <- file.path("classifiers", "TP53")
results <- parse_summary(file.path(results_folder, "classifier_summary.txt"))

# 1) Heatmap of the distribution of aberrant events across tumors
heatmap_plot_file <- file.path(results_folder, "figures", "tp53_heatmap.pdf")
heat_file <- file.path(results_folder, "summary_counts.csv")
heat_df <- readr::read_csv(heat_file)

prop_matrix <- as.matrix(heat_df[, c('TP53_loss_y', 'TP53_y')])
rownames(prop_matrix) <- heat_df$DISEASE
colnames(prop_matrix) <- c("Loss", "Mutation")

# All diseases that are used in building the classifier
tp53_dis <- results[["Diseases"]]

# Build a vector for heatmap labels
classifier <- c()
for (disease in rownames(prop_matrix)) {
  if (disease %in% tp53_dis) {
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

# 2) Coefficients contributing to the model
coef_plot_file <- file.path(results_folder, "figures", "ddr_coef_plot.svg")
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

p <- add_arrow_label(p = p, x = 1050, y = -0.205, label = "DDB2",
                     offset = c(80, 0.001, -645, -.0002))
p <- add_arrow_label(p = p, x = 1150, y = -0.190, label = "AEN",
                     offset = c(80, 0.001, -450, 0))
p <- add_arrow_label(p = p, x = 1480, y = -0.172, label = "RPS27L",
                     offset = c(80, 0.001, -890, .0003))
p <- add_arrow_label(p = p, x = 950, y = -0.155, label = "MDM2",
                     offset = c(80, -0.001, -680, -.00013))
p <- add_arrow_label(p = p, x = 1700, y = -0.14, label = "BAX",
                     offset = c(80, -0.001, -510, .0002))
p <- add_arrow_label(p = p, x = 2000, y = -0.12, label = "CDKN1A",
                offset = c(90, 0.0001, -950, -.00015))
p <- add_arrow_label(p = p, x = 1800, y = -0.1, label = "XPC",
                     offset = c(80, 0, -550, -.0002))
p <- add_arrow_label(p = p, x = 2100, y = -0.085, label = "MPDU1",
                     offset = c(80, 0, -800, -.0002))
p <- add_arrow_label(p = p, x = 1800, y = -0.07, label = "FDXR",
                     offset = c(80, 0, -650, -.0006))
p <- add_arrow_label(p = p, x = 1900, y = -0.055, label = "PHLDA3",
                     offset = c(80, 0, -890, -.0002))

p <- add_arrow_label(p = p, x = 7000, y = 0.1, label = "KIF1B",
                     offset = c(-80, .001, 550, -.006))
p <- add_arrow_label(p = p, x = 6700, y = 0.085, label = "EEPD1",
                     offset = c(-80, .0006, 460, -.006))
p <- add_arrow_label(p = p, x = 6400, y = 0.07, label = "MIIP",
                     offset = c(-80, .0004, 500, -.0015))
p <- add_arrow_label(p = p, x = 6050, y = 0.054, label = "ID4",
                     offset = c(-80, .0004, 490, -.0015))
p <- add_arrow_label(p = p, x = 5750, y = 0.036, label = "CDC123",
                     offset = c(-80, .0004, 900, -.0015))
p <- add_arrow_label(p = p, x = 5350, y = 0.021, label = "DDX27",
                     offset = c(-80, .0004, 790, -.0015))
p <- add_arrow_label(p = p, x = 6750, y = 0.012, label = "DCAF13",
                     offset = c(-80, -.001, 570, .007))

p <- add_arrow_label(p = p, x = 6500, y = -0.03, label = "log10_mut",
                     offset = c(-50, -0.002, 900, 0.005))
ggsave(coef_plot_file, plot = p, height = 2.5, width = 2.25)

# 3) Plot distributions of predictions according to variant classification
var_plot_file <- file.path(results_folder, "figures", "variant_fill_map.svg")
mut_df <- readr::read_tsv(file.path(results_folder, "tables",
                                    "mutation_classification_scores.tsv"))

consider_mutations <- c("3'UTR", "5'UTR", "Intron", "Frame_Shift_Del",
                        "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                        "Missense_Mutation", "Nonsense_Mutation",
                        "Nonstop_Mutation", "RNA", "Splice_Site")

silent_df <- mut_df %>% filter(Variant_Classification == "Silent") %>%
  filter(total_status == 0)
delet_df <- mut_df %>% filter(Variant_Classification %in% consider_mutations)
mut_filtered_df <- dplyr::bind_rows(delet_df, silent_df)

# Separate classes of mutations to summarize
copy_num_df <- mut_df %>% filter(TP53_loss == 1) %>%
  filter(TP53 == 0) %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Loss")
missense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Missense_Mutation") %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Missense")
nonsense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Nonsense_Mutation") %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Nonsense")
indel_df <- mut_filtered_df %>% filter(Variant_Classification %in%
                                       c("Frame_Shift_Del", "Frame_Shift_Ins",
                                         "In_Frame_Del", "In_Frame_Ins")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Indel")
utr_df <- mut_filtered_df %>%
  filter(Variant_Classification %in% c("3'UTR", "5'UTR", "Intron")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "UTR")
silent_df <- silent_df %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Silent")
splice_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Splice_Site") %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Splice")
wt_df <- mut_df %>% subset(total_status == 0) %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "WT")
hyper_df <- mut_df %>%
  filter(hypermutated == 1) %>%
  select(Variant_Classification, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Hyper")

final_df <- dplyr::bind_rows(list(missense_df, nonsense_df, indel_df, utr_df,
                                  splice_df, silent_df, copy_num_df, wt_df,
                                  hyper_df))

colnames(final_df) <- c("ID", "Disease", "Weight", "HGVSc", "HGVSp", "Class")

# Plot summary distribution of variant classes prediction scores
ggplot(final_df, aes(Weight, ..count.., fill = Class)) +
  geom_density(position = "fill", size = 0.1) +
  geom_segment(aes(x = 0.5, y = 0, yend = 1, xend = 0.5), colour = "black",
               linetype = "dashed", size = 0.4) +
  labs(list(x = "Probability", y = "Proportion")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + base_theme +
  theme(legend.position = c(1.1, 0.65),
        legend.background = element_rect(fill = alpha('white', 0)),
        legend.text = element_text(size = 7),
        plot.margin = unit(c(0.2, 1.5, 0, 0.1),"cm"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 12))
ggsave(var_plot_file, width = 4, height = 3.8)

# 4) Show mutation frequencies and scores
mut_weight_df <- mut_filtered_df %>% filter(!is.na(weight))
mut_weight_df <- mut_weight_df[mut_weight_df$hypermutated != 1, ]

aa_df <- mut_weight_df %>%
  group_by(HGVSp, Variant_Classification) %>%
  summarise(Mean = mean(weight),
            low_CI = quantile(weight, 0.05),
            high_CI = quantile(weight, 0.95),
            count = n())
nuc_df <- mut_weight_df %>%
  group_by(HGVSc, Variant_Classification) %>%
  summarise(Mean = mean(weight),
            low_CI = quantile(weight, 0.05),
            high_CI = quantile(weight, 0.95),
            count = n())

aa_df <- aa_df[order(aa_df$count, decreasing = TRUE),]
nuc_df <- nuc_df[order(nuc_df$count, decreasing = TRUE),]
write.table(aa_df, file = file.path(results_folder, 'tables',
                                    'amino_acid_mutation_scores.tsv'),
            sep = '\t', row.names = FALSE)
write.table(nuc_df, file = file.path(results_folder, 'tables',
                                     'nucleotide_mutation_scores.tsv'),
            sep = '\t', row.names = FALSE)
