args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}









load_in_ldsc_like_method_data <- function(ldsc_like_methods_results_dir) {
	n_sims <- 50
	sim_sampling_variances <- c("0.002", "0.004", "0.006", "0.008")

	est_mean_arr <- c()
	est_mean_lb_arr <- c()
	est_mean_ub_arr <- c()
	sim_mean_arr <- c()

	for (samp_var_iter in 1:length(sim_sampling_variances)) {

		sim_sampling_variance <- sim_sampling_variances[samp_var_iter]

		est_vars <- c()
		for (sim_iter in 1:n_sims) {
			filer <- paste0(ldsc_like_methods_results_dir, "sim", sim_iter, "_sim_eqtl_ss_489_simple_mog_", sim_sampling_variance, "_ldsc_like_est_of_sampling_variance.txt")
			tmp_df <- read.table(filer, header=TRUE)

			est_vars <- c(est_vars, tmp_df$borzoi_sampling_variance[1])
		}

		est_mean <- mean(est_vars)
		est_mean_se <- sd(est_vars) / sqrt(length(est_vars))

		est_mean_lb <- est_mean - (1.96*est_mean_se)
		est_mean_ub <- est_mean + (1.96*est_mean_se)


		est_mean_arr <- c(est_mean_arr, est_mean)
		est_mean_lb_arr <- c(est_mean_lb_arr, est_mean_lb)
		est_mean_ub_arr <- c(est_mean_ub_arr, est_mean_ub)
		sim_mean_arr <- c(sim_mean_arr, as.numeric(sim_sampling_variance))

	}

	df <- data.frame(
 		est_mean = est_mean_arr,
 		est_mean_lb = est_mean_lb_arr,
  		est_mean_ub = est_mean_ub_arr,
  		sim_mean = sim_mean_arr
	)

	return(df)
}

make_ldsc_plot <- function(ldsc) {
  p <- ggplot(ldsc, aes(x = sim_mean, y = est_mean)) +
    geom_point(size = 3, color = "#2C7FB8") +
    geom_errorbar(
      aes(ymin = est_mean_lb, ymax = est_mean_ub),
      width = 0,
      size = 1,
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    labs(
      x = "Simulated sampling variance",
      y = "Estimated sampling variance"
    ) +
    figure_theme()

  return(p)
}




####################
# Command line args
####################
ldsc_like_methods_results_dir <- args[1]
viz_dir <- args[2]


# Load in data
ldsc_df <- load_in_ldsc_like_method_data(ldsc_like_methods_results_dir)

# Make ldsc-like plot
pp <- make_ldsc_plot(ldsc_df)
output_file <- paste0(viz_dir, "estimated_vs_simulated_sampling_variance_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")

print(output_file)


