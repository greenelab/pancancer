# Gregory Way 2017
# PanCancer Classifier
# scripts/compare_within_models.R
#
# Plots pancancer classifier performance compared to within cancer
#
# Usage: Run by assigning where the within classifier summary is and where the
#        Pan Cancer classifier summary is
#
#     R --no-save --args <within_tissue_folder> <pan_cancer_folder>
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

process_classifier_summary <- function(summary_list, model_type) {
  # Takes in a parsed classifier summary list and outputs a processed dataframe
  #
  # summary_list - a list storing classifier attributes and performance
  # model_type - a string that will indicate the type of model

  disease_perf <- summary_list[["Disease specific performance"]]
  pancan_df <- data.frame(disease_perf[, "disease"])
  pancan_df$Gene <- summary_list[["Genes"]]
  pancan_df$AUROC <- summary_list[["Disease specific performance"]][, "cv"]
  pancan_df$Model <- model_type
  colnames(pancan_df)[1] <- "Disease"
  return(pancan_df)
}

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
  
ggsave(file.path(pan_summary_dir, "figures", "comparison.pdf"), units = "in",
       height = 1.4, width = 4.2, dpi = 600)
