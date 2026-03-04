args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}




scatter_burden_effect_vs_non_coding_effect <- function(sub_df) {
  sub_df$neg_non_coding_burden    <- -sub_df$non_coding_burden_effect
  sub_df$neg_non_coding_burden_se <-  sub_df$non_coding_burden_effect_se  # SE unchanged by sign flip

  z <- 1.96  # 95% CI

  pp <- ggplot(sub_df, aes(x = neg_non_coding_burden, y = effect_size)) +
    # vertical 95% CI (y direction)
    geom_errorbar(
      aes(ymin = effect_size - z*effect_size_se,
          ymax = effect_size + z*effect_size_se),
      width = 0, alpha = 0.6
    ) +
    # horizontal 95% CI (x direction) using geom_errorbar + orientation
    geom_errorbar(
      aes(xmin = neg_non_coding_burden - z*neg_non_coding_burden_se,
          xmax = neg_non_coding_burden + z*neg_non_coding_burden_se),
      orientation = "y",
      width = 0, alpha = 0.6
    ) +
    geom_point(color = "dodgerblue3", size = 2, alpha = 0.85) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
    geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
    labs(x = "negative non-coding effect", y = "PTV Burden effect size") +
    figure_theme()

  return(pp)
}

scatter_burden_effect_vs_non_coding_effect_no_errorbar <- function(sub_df) {
  sub_df$neg_non_coding_burden <- -sub_df$non_coding_burden_effect
  sub_df$abs_zed <- abs(sub_df$non_coding_burden_effect /
                        sub_df$non_coding_burden_effect_se)

  pp <- ggplot(
      sub_df,
      aes(x = neg_non_coding_burden,
          y = effect_size,
          size = abs_zed)
    ) +
    geom_point(color = "dodgerblue3", alpha = 0.85) +
    scale_size_continuous(
      name = "|z|",
      range = c(1.5, 6)
    ) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
    geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
    labs(x = "negative non-coding effect",
         y = "PTV Burden effect size") +
    figure_theme()

  return(pp)
}





####################
# Command line args
####################
burden_results_file <- args[1]
visualization_dir <- args[2]


# Load in data
df <- read.table(burden_results_file, header=TRUE, sep="\t")
sub_df <- df[df$borzoi_thresh==0.1,]




# MAke scatter plot
pp <- scatter_burden_effect_vs_non_coding_effect(sub_df)
output_file <- paste0(visualization_dir, "effect_size_scatter.pdf")
ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")

# MAke scatter plot
pp <- scatter_burden_effect_vs_non_coding_effect_no_errorbar(sub_df)
output_file <- paste0(visualization_dir, "effect_size_scatter_no_error_bar.pdf")
ggsave(pp, file=output_file, width=7.2, height=5.0, units="in")
