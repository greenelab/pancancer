# Gregory Way 2017
# PanCancer Classifier
# scripts/viz/braf_summary_figures.R
#
# Visualize summary for Ras Classifier Scores applied to BRAFV600E mutations
#
# Usage: Run by assigning where the within classifier summary is and where the
#        Pan Cancer classifier summary is
#
#     Rscript --vanilla scripts/viz/braf_summary_figures.R
#
# Output:
# Several figures to summarize BRAF findings

checkpoint::checkpoint("2017-06-01", checkpointLocation = ".")

library(dplyr)
library(pheatmap)
library(ggplot2)
library(readr)
library(cowplot)
library(gridExtra)
source(file.path("scripts", "util", "pancancer_util.R"))

results_folder <- file.path("classifiers", "RAS_noTHCASKCM")
results <- parse_summary(file.path(results_folder, "classifier_summary.txt"))

# 1) Plot distributions of predictions according to variant classification
mut_df <- readr::read_tsv(file.path(results_folder, "tables",
                                    "mutation_classification_scores.tsv"))

consider_mutations <- c("3'UTR", "5'UTR", "Intron", "Frame_Shift_Del",
                        "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                        "Missense_Mutation", "Nonsense_Mutation",
                        "Nonstop_Mutation", "RNA", "Splice_Site")

# Generate RAS Gain Variable
mut_df <- mut_df %>% mutate(RAS_gain = KRAS_gain + NRAS_gain + HRAS_gain)
mut_df$RAS_gain[mut_df$RAS_gain > 1] <- 1

mut_df <- mut_df %>% mutate(RAS_mutation = KRAS + NRAS + HRAS)
mut_df$RAS_mutation[mut_df$RAS_mutation > 1] <- 1

silent_df <- mut_df %>% filter(Variant_Classification == "Silent") %>%
  filter(total_status == 0)
delet_df <- mut_df %>% filter(Variant_Classification %in% consider_mutations)

mut_filtered_df <- dplyr::bind_rows(delet_df, silent_df)

# Separate classes of mutations to summarize
copy_num_df <- mut_df %>% filter(RAS_gain == 1) %>%
  filter(RAS_mutation == 0) %>%
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
var_plot_file <- file.path(results_folder, "figures", "variant_distrib.svg")
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

# 2) Show mutation frequencies and scores
mut_weight_df <- mut_filtered_df %>% filter(!is.na(weight))
mut_weight_df <- mut_weight_df[mut_weight_df$hypermutated != 1, ]

aa_df <- mut_weight_df %>%
  group_by(HGVSp, Variant_Classification, Hugo_Symbol) %>%
  summarise(Mean = mean(weight),
            low_CI = quantile(weight, 0.05),
            high_CI = quantile(weight, 0.95),
            count = n())
nuc_df <- mut_weight_df %>%
  group_by(HGVSc, Variant_Classification, Hugo_Symbol) %>%
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

# 3) Plot summary distribution of BRAFV600E mutations
braf_df <- final_df[complete.cases(final_df), ]
braf_df <- braf_df[braf_df$HGVSp == 'p.Val600Glu', ]

braf_df$Disease <- dplyr::recode(braf_df$Disease,
                                 "BLCA" = "Other", "CHOL" = "Other",
                                 "GBM" = "Other", "HNSC" = "Other",
                                 "KIRP" = "Other", "LGG" = "Other",
                                 "READ" = "Other")

braf_plot_file <- file.path(results_folder, 'figures',
                            'brafv600e_distribution.svg')
braf_plot <- ggplot(braf_df, aes(Weight, fill = Disease)) +
  geom_density(alpha = 0.4) + theme_bw() +
  ylab("Density") + xlab("BRAFV600E Classifier Score")

svg(braf_plot_file, width = 4, height = 3)
braf_plot
dev.off()
