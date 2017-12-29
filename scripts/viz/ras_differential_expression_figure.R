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
# A single plot for supplementary figure S3

library(ggplot2)
library(ggpmisc)
library(ggrepel)

set.seed(123)

# Load Data
diff_exp_file <- file.path('classifiers', 'RAS', 'differential_expression.tsv')
plot_df <- readr::read_tsv(diff_exp_file)

color_logic <- (plot_df$weight > 0.05 | plot_df$weight < -0.05) |
  (plot_df$stat > 28 | plot_df$stat < -16)

ggplot(plot_df, aes(x = weight, y = stat)) +
  geom_point(alpha = 0.5, color = ifelse(color_logic, 'red', 'grey50')) +
  xlab('Ras Classifier Weight') +
  geom_text_repel(data = subset(plot_df,
                                (weight > 0.05 | weight < -0.05) |
                                  (stat > 28 | stat < -16)),
                  arrow = arrow(length = unit(0.01, 'npc')),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  size = 1.5,
                  fontface = 'italic',
            aes(x = weight, y = stat, label = gene)) +
  ylab('Differential Expression Score') +
  theme_bw() +
  theme(axis.text = element_text(size = rel(0.65)),
        axis.title = element_text(size = rel(0.8)),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 3, r = 0, b = 0, l = 0)),
        plot.margin = margin(r = 0.3, l = 0.1, unit = 'cm'))

diff_exp_fig <- file.path('classifiers', 'RAS', 'figures', 'diff_exprs.pdf')
ggplot2::ggsave(diff_exp_fig, dpi = 600, width = 2.3, height = 2)
