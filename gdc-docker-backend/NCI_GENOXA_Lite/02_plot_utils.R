# plot_utils.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("03_utils.R")

plot_raw_gene_counts <- function(raw_gene_counts, exp_design, rcs_range, count_type, plot_type, color_palette_func, 
                                 show_x_labels = TRUE, show_y_labels = TRUE, x_axis_label = "Sample Name", 
                                 y_axis_label = "Gene Counts", plot_title = "Raw Gene Count Data Plot", 
                                 legend_title = "Sample_Type", show_plot_title = TRUE, show_legend = TRUE) {
  
  r_cnt <- as.data.frame(t(raw_gene_counts[rcs_range]))
  grp <- exp_design[,2]
  r_cnt$Sample_Type <- grp[rcs_range]
  
  r_cnt <- r_cnt %>%
    rownames_to_column(var = "Sample") %>%
    mutate(Sample = factor(Sample, levels = rownames(r_cnt)))  
  
  r_cnt_long <- r_cnt %>%
    pivot_longer(cols = -c(Sample, Sample_Type), names_to = "Gene", values_to = "Count")
  
  r_cnt_long$log2_Count <- log2(r_cnt_long$Count + 1)
  
  y_var <- if (count_type == "Count") {
    "Count"
  } else {
    "log2_Count"
  }
  
  color_palette <- color_palette_func(n = length(unique(r_cnt_long$Sample_Type)))
  
  p <- ggplot(r_cnt_long, aes(x = Sample, y = !!sym(y_var), fill = Sample_Type)) +
    ggtitle(if (show_plot_title) plot_title else NULL) +
    xlab(if (show_x_labels) x_axis_label else NULL) + 
    ylab(if (show_y_labels) y_axis_label else NULL) +  
    scale_fill_manual(values = color_palette, name = legend_title) +  # Use custom legend title
    extra_theme() +
    theme(
      plot.title = if (show_plot_title) element_text(color = "black", size = 12, face = "bold", hjust = 0.5) else element_blank(),
      axis.title.x = if (show_x_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.title.y = if (show_y_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.text.x = if (show_x_labels) element_text(angle = 90, hjust = 1) else element_blank(),
      axis.text.y = if (show_y_labels) element_text(size = 8) else element_blank(),
      legend.title = if (show_legend) element_text(size = 10, face = "bold") else element_blank(),
      legend.text = if (show_legend) element_text(size = 8) else element_blank(),
      legend.position = if (show_legend) "right" else "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  if (plot_type == "boxplot") {
    p <- p + geom_boxplot()
  } else if (plot_type == "col") {
    p <- p + geom_col()
  } else if (plot_type == "violin") {
    p <- p + geom_violin()
  }
  
  return(p)
}


