# Gregory Way 2017
# PanCancer Classifier
# scripts/compare_within_models.R
#
# Plots pancancer classifier performance compared to within cancer
#
# Usage: Run by assigning where the within classifier summary is and where the
#        Pan Cancer classifier summary is
#
#     Rscript scripts/compare_within_models.R 
#
#     with the required flags:
#         --pancan_summary      Directory of where classifier summary is
#         --within_dir          Directory of where within cancer-type data is
#
# Output:
# Bar Plots for each comparison

library(ggplot2)
source(file.path("scripts", "util", "pancancer_util.R"))

option_list <- list(optparse::make_option(c("-w", "--within_dir"),
                                          type = "character",
                                          help = "Location of within files"),
                    optparse::make_option(c("-p", "--pancan_summary"),
                                          type = "character",
                                          help = "location of pancan summary"))

opt_parser <- optparse::OptionParser(option_list = option_list);
opt <- optparse::parse_args(opt_parser);

within_folder <- opt$within_dir
pan_summary_dir <- opt$pancan_summary

# Process PanCancer Classifier and summary files
pan_summary <- file.path(pan_summary_dir, "classifier_summary.txt")
pancan_list <- parse_summary(pan_summary)
pancan_df <- process_classifier_summary(pancan_list, "Pan")

# Process Within Cancer Results
within_tissue_files <- list.files(within_folder,
                                  pattern = "classifier_summary.txt",
                                  full.names = TRUE, recursive = TRUE)
within_tissue_data <- data.frame()
for (file in within_tissue_files) {
  file_summary <- parse_summary(file)
  file_frame <- process_classifier_summary(file_summary, "Within")
  within_tissue_data <- rbind(within_tissue_data, file_frame)
}

plot_ready <- plyr::rbind.fill(within_tissue_data, pancan_df)
plot_ready$AUROC <- as.numeric(paste(plot_ready$AUROC))
plot_ready <- plot_ready[complete.cases(plot_ready), ]

ggplot(plot_ready, aes(x = Disease, y = AUROC, fill = Model)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() + within_theme +
  theme(legend.position = c(1.07, 0.65),
        legend.background = element_rect(fill = alpha("white", 0)),
        plot.margin = unit(c(0.2, 1.5, 0, 0.1), "cm")) +
  scale_fill_manual(values = c("brown", "gold")) +
  geom_hline(yintercept = 0.5, linetype = "longdash", size = 0.4,
             color = "black") +
  ylab("CV AUROC") +
  coord_cartesian(ylim = c(0.4, 1)) +
  scale_y_continuous(breaks = seq(0.4, 1, 0.1))

ggsave(file.path(pan_summary_dir, "comparison.svg"), units = "in",
       height = 1.4, width = 4.2)
