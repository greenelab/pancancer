# Gregory Way 2017
# PanCancer Classifier
# pancancer_util.R
#
# Custom ggplot themes, function to extract info from a classifier summary
# file, and several helper functions for plotting.
#
# Usage: sourced by scripts using function

checkpoint::checkpoint("2017-06-01", checkpointLocation = ".")

library(ggplot2)

# ggplot2 themes to be used in various contexts
base_theme <- ggplot2::theme(title = element_text(size = rel(0.6)),
                             axis.title = element_text(size = rel(0.9)),
                             axis.text.x = element_text(size = rel(0.7),
                                                        vjust = 2),
                             axis.text.y = element_text(size = rel(0.7)),
                             axis.line.x = element_blank(),
                             axis.line.y = element_blank(),
                             axis.ticks = element_blank(),
                             legend.text = element_text(size = rel(0.4)),
                             legend.key = element_blank(),
                             legend.key.size = unit(0.4, "lines"),
                             legend.position = c(1.098, 0.7),
                             panel.grid.major = element_line(color = "white",
                                                             size = 0.3),
                             panel.grid.minor = element_line(color = "white",
                                                             size = 0.3),
                             panel.background = element_rect(fill = "white"),
                             plot.margin = unit(c(0.2, 1.5, 0, 0),"cm")) 

within_theme <- theme(title = element_text(size = rel(0.7)),
                      axis.text.x = element_text(angle = 90, size = rel(0.5)),
                      axis.text.y = element_text(size = rel(0.5)),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size = rel(0.7)),
                      legend.position = 'right',
                      plot.margin = unit(rep(0.1, 4), "cm"),
                      legend.text = element_text(size = rel(0.45)),
                      legend.key = element_blank(),
                      legend.key.size = unit(2, "mm"), 
                      strip.text.x = element_text(size = rel(1.0)),
                      panel.grid.major = element_line(color = "gray",
                                                      size = 0.3),
                      panel.background = element_rect(fill = "white"))

parse_summary <- function(summary_info) {
  # Process classifier summary file
  # Arguments:
  #   summary_info: either a dataframe or string storing classifier info
  # Output:
  #    a list of summarized classifier attributes and performance

  if (is.character(summary_info)) {
    summary_info <- readr::read_lines(summary_info)
  }
  summary_list <- list()
  dis_spec_perf <- c()
  alt_spec_perf <- c()
  alt_gene <- FALSE
  for (line in summary_info) {
    line <- unlist(strsplit(line, "\t"))
    if (is.na(line[1]) | (line[1] %in% c("Parameters:", "Results:",
                                         "Disease specific performance:",
                                         "Alternative gene AUROC:"))) {
      if (is.na(line[1])) {
        next
      }
      if (line[1] == "Alternative gene AUROC:") {
        alt_gene <- TRUE
      }
      next
    }
    if (line[1] == "Coefficients:") {
      summary_list[[sub(":", "", line[1])]] <-
        suppressMessages(readr::read_tsv(line[2]))
    } else if (line[1] == "") {
      disease_info <- line[2:length(line)]
      disease <- disease_info[1]
      perf_type <- gsub(":", "", unlist(strsplit(disease_info[2], " "))[2])
      train <- disease_info[3]
      if (alt_gene) {
        test <- disease_info[3]
      } else {
        test <- disease_info[5]
      }
      cv <- disease_info[7]
      disease_summary <- c(disease, train, test, cv, perf_type)
      if (alt_gene) {
        alt_spec_perf <- rbind(alt_spec_perf, disease_summary)
      } else {
        dis_spec_perf <- rbind(dis_spec_perf, disease_summary)
      }
    } else {
      summary_list[[gsub(":", "", line[1])]] <- line[2:length(line)]
    }
  }
  colnames(dis_spec_perf) <- c("disease", "training", "testing", "cv",
                               "performance_type")
  summary_list[["Disease performance"]] <- data.frame(dis_spec_perf)
  if (alt_gene) {
    colnames(alt_spec_perf) <- c("disease", "holdout", "cv", "data_type",
                                 "performance_type")
    summary_list[["Alt gene performance"]] <- data.frame(alt_spec_perf)
  }
  
  return(summary_list)
}

process_classifier_summary <- function(summary_list, model_type,
                                       gene_type = "Disease performance",
                                       gene_class = "Genes",
                                       perf_type = "AUROC") {
  # Takes in a parsed classifier summary list and outputs a processed dataframe
  #
  # summary_list - a list storing classifier attributes and performance
  # model_type - a string that will indicate the type of model
  
  disease_perf <- summary_list[[gene_type]]
  disease_perf <- disease_perf[disease_perf$performance_type == perf_type, ]
  pancan_df <- data.frame(disease_perf[, "disease"])
  pancan_df$Gene <- paste(summary_list[[gene_class]], collapse = "_")
  pancan_df$Performance_type <- disease_perf[, "cv"]
  pancan_df$Model <- model_type
  colnames(pancan_df)[1] <- "Disease"
  return(pancan_df)
}

add_arrow_label <- function(p, x, y, label, offset = c(0, 0, 0, 0)) {
  # Add an arrow and label to a point in a plot
  # Arguments:
  #    p: a ggplot2 object to add annotations to
  #    x: the x coordinate to place the text
  #    y: the y coordinate to place the text
  #    label: the feature of interest (string)
  #    offset: how much to tweak the position of the arrows
  # Output:
  #    an update ggplot2 object

  offset_x_point <- offset[1]
  offset_y_point <- offset[2]
  offset_x_label <- offset[3]
  offset_y_label <- offset[4]
  point <- coef_df[coef_df$feature == label, ]
  x_point <- point$rank
  y_point <- point$weight
  p <- p + annotate("text", x = x, y = y, label = label, size = 2.2,
                    fontface = "bold.italic") +
    annotate("segment", x = x + offset_x_label, y = y + offset_y_label,
             xend = x_point + offset_x_point, yend = y_point + offset_y_point,
             size = 0.18, arrow = arrow(length = unit(0.05, "cm")))
  return(p)
}

get_boot <- function(weight, low = TRUE, B = 1000, min_samples = 5) {
  # Function to input into a `dplyr::summarize()` call
  # The function will compute a bootstrap estimate of the mean as well as
  # the resulting 95% bootstrapped confidence intervals. The call to summarize
  # requires a single returned numeric value, hence the separation of low_ and
  # high_ returned objects.
  #
  # Arguments:
  #    weight: a vector of Ras pathway classifier predictions
  #    low: boolean of whether to return 95% low or high confidence intervals
  #    B: integer of how many bootstrap samples to take, defaults to 1000
  #    min_samples: bootstrapping does not work well with low sample sizes, how
  #                 many samples are required to calculate a bootstrap estimate

  if (length(weight) < min_samples) {
    return(NA)
  }

  boot_obj <- Hmisc::smean.cl.boot(weight, B = B)

  low_ <- as.numeric(boot_obj[2])
  high_ <- as.numeric(boot_obj[3])

  if (low) {
    return(low_)
  } else {
    return(high_)
  }
}
