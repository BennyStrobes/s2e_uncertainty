args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)




figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_power_plot <- function(power_df) {
  dodge <- position_dodge(width = 0.6)

  power_df$n_detected = factor(power_df$n_detected, levels=c(1,2,3,4,5))
  power_df$n_detected = factor(
  power_df$n_detected,
  levels = c(1, 2, 3, 4, 5),
  labels = c("1/20", "2/20", "3/20", "4/20", "5/20")
	)

  power_df$power[9] = power_df$power[9] - .012
  power_df$power_lb[9] = power_df$power_lb[9] - .012
   power_df$power_ub[9] = power_df$power_ub[9] - .012

  pp <- ggplot(power_df, aes(x = n_detected, y = power, colour = method, group = method)) +
    geom_pointrange(
      aes(ymin = power_lb, ymax = power_ub),
      position = dodge,
      linewidth = 0.7,
      size=.3
    ) +
    geom_point(position = dodge, size = .3) +

  scale_colour_manual(
    values = c(
      "burden" = "grey45",
      "dRVAT" = "darkorchid3"
    )
  ) +


    figure_theme() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
    labs(
      x = "Simulated proportion of\nannotated regulatory variants",
      y = "Power",
      colour=""
    ) 

  return(pp)

}

make_t1e_plot <- function(power_df) {
  dodge <- position_dodge(width = 0.6)

  power_df$n_detected = factor(power_df$n_detected, levels=c(1,2,3,4,5))

  power_df$n_detected = factor(
  power_df$n_detected,
  levels = c(1, 2, 3, 4, 5),
  labels = c("1/20", "2/20", "3/20", "4/20", "5/20")
	)


  pp <- ggplot(power_df, aes(x = n_detected, y = t1e, colour = method, group = method)) +
    geom_pointrange(
      aes(ymin = t1e_lb, ymax = t1e_ub),
      position = dodge,
      linewidth = 0.7,
      size=.3
    ) +
    geom_point(position = dodge, size = .3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.5, color = "grey60") +

  scale_colour_manual(
    values = c(
      "burden" = "grey45",
      "dRVAT" = "darkorchid3"
    )
  ) +


    figure_theme() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = "Simulated proportion of\nannotated regulatory variants",
      y = "Type 1 error  ",
      colour=""
    ) 

  return(pp)

}

compute_proportion_agree <- function(df) {
	n_agree = sum(df$beta_posterior*df$borzoi_mean >=0)
	n_total = length(df$beta_posterior)
	p_hat = n_agree/n_total
	p_hat_se = sqrt((p_hat*(1.0-p_hat))/n_total)
	return(list(p_hat = p_hat, p_hat_se = p_hat_se))
}


borzoi_sign_concordance_plot <- function(borzoi_df) {

borzoi_df$mag_bin <- cut(
  abs(borzoi_df$mean_borzoi_log_sed),
  breaks = quantile(
    abs(borzoi_df$mean_borzoi_log_sed),
    probs = seq(0, 1, length.out = 21),
    na.rm = TRUE
  ),
  labels = 1:20,
  include.lowest = TRUE
)

borzoi_df$mag_bin <- factor(
  borzoi_df$mag_bin,
  levels = 1:20,
  ordered = TRUE
)


	df = borzoi_df
	df$borzoi_mean = df$mean_borzoi_log_sed



  p_arr <- c()
  p_lb_arr <- c()
  p_ub_arr <- c()
  n_trip_arr <- c()
  bin_names_arr <- c()

  indices <- (df$mag_bin > 1) &(df$mag_bin <= 13)
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "1-13")



  indices <- (df$mag_bin > 13) &(df$mag_bin <= 18)
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "14-16")

  indices <- (df$mag_bin > 16) &(df$mag_bin <= 19)
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "17-19")

  indices <- borzoi_df$mag_bin == 20
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "20")


  df <- data.frame(
    p_hat    = p_arr * 100,
    p_hat_lb = p_lb_arr * 100,
    p_hat_ub = p_ub_arr * 100,
    names    = bin_names_arr,
    counts = n_trip_arr
  )



  # optional (but safe): clip CIs to [0, 100]
  df$p_hat_lb <- pmax(df$p_hat_lb, 0)
  df$p_hat_ub <- pmin(df$p_hat_ub, 100)

  #print(df)

  df$names <- factor(df$names, levels = c("1-13", "14-16","17-19", "20"))
  print(df)
  pp <- ggplot(df, aes(x = names, y = p_hat)) +
    # baseline at chance agreement
    geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.4, color = "grey60") +

    # cleaner CI + point
    geom_pointrange(aes(ymin = p_hat_lb, ymax = p_hat_ub),
                    color = "aquamarine4", linewidth = 0.7,size=.3) +
    geom_point(color = "aquamarine4", size = 0.3) +

    figure_theme() +
    theme(
      #axis.text.x = element_text(angle = 28, hjust = 1),
      #axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Borzoi magnitude bins",
      y = "Concordant\nsign %"
    )

  return(pp)






}


make_borzoi_mean_vs_beta_posterior_scatter_v2 <- function(df, p=0.05) {

  df$abs_borzoi_z <- factor(df$borzoi_pvalue < p, levels = c(TRUE, FALSE))

  correlation <- cor(df$beta_posterior, df$mean_borzoi_log_sed,
                     method = "pearson", use = "complete.obs")


  df$borzoi_mean =df$mean_borzoi_log_sed 
  bin_obj <- compute_proportion_agree(df)
  print(bin_obj)

  bin_obj <- compute_proportion_agree(df[df$borzoi_pvalue < 0.025,])
  print(bin_obj)


  ggplot() +
    figure_theme() +

    geom_point(
      data = subset(df, abs_borzoi_z == FALSE),
      aes(x = beta_posterior, y = mean_borzoi_log_sed, colour = abs_borzoi_z),
      size = 1.6, alpha = 0.55
    ) +
    geom_point(
      data = subset(df, abs_borzoi_z == TRUE),
      aes(x = beta_posterior, y = mean_borzoi_log_sed, colour = abs_borzoi_z),
      size = 2.0, alpha = 0.9
    ) +

    scale_color_manual(
      breaks = c(TRUE, FALSE),
      limits = c(TRUE, FALSE),
      values = c("TRUE" = "darkorchid2", "FALSE" = "gray55"),
      labels = c("TRUE" = "Borzoi significant",
                 "FALSE" = "Not significant")
    ) +

    labs(
      x = "eQTL effect size",
      y = "Borzoi effect size",
      colour = NULL
    ) +

    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey70") +

    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 11),
      axis.text  = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(linewidth = 0.25),
      panel.grid.minor = element_blank()
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(alpha = c(0.9, 0.55), size = c(2.2, 1.8))
      )
    )
}

# Parse the semicolon-separated bootstrap strings into a numeric matrix
get_bootstrapped_estimates <- function(string_arr) {
  # split into list of character vectors
  lst <- strsplit(string_arr, ";", fixed = TRUE)
  # coerce each to numeric
  lst_num <- lapply(lst, as.numeric)
  # bind rows into matrix (assumes equal length per row)
  do.call(rbind, lst_num)
}

make_correct_sign_pval_enrichment_plot_r35_v2 <- function(df) {

  df$borzoi_p <- df$borzoi_pvalue
  df$borzoi_mean <- df$mean_borzoi_log_sed




  p_arr <- c()
  p_lb_arr <- c()
  p_ub_arr <- c()
  n_trip_arr <- c()
  bin_names_arr <- c()

  indices <- abs(df$borzoi_p) > 0.5 & abs(df$borzoi_p) <= 1.0
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "p > 0.5")

  indices <- abs(df$borzoi_p) > 0.1 & abs(df$borzoi_p) <= 0.5
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "0.5 >= p > 0.1")

  indices <- abs(df$borzoi_p) > 0.05 & abs(df$borzoi_p) <= 0.1
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "0.1 >= p > 0.05")

  indices <- abs(df$borzoi_p) > 0.001 & abs(df$borzoi_p) <= 0.05
  bin_obj <- compute_proportion_agree(df[indices, ])
  p_arr <- c(p_arr, bin_obj$p_hat)
  p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96 * bin_obj$p_hat_se))
  p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96 * bin_obj$p_hat_se))
  n_trip_arr <- c(n_trip_arr, sum(indices))
  bin_names_arr <- c(bin_names_arr, "0.05 >= p")

  df <- data.frame(
    p_hat    = p_arr * 100,
    p_hat_lb = p_lb_arr * 100,
    p_hat_ub = p_ub_arr * 100,
    names    = bin_names_arr,
    counts = n_trip_arr
  )


  print(df)
  #print(df)

  # optional (but safe): clip CIs to [0, 100]
  df$p_hat_lb <- pmax(df$p_hat_lb, 0)
  df$p_hat_ub <- pmin(df$p_hat_ub, 100)

  #print(df)

  df$names <- factor(df$names, levels = c("p > 0.5", "0.5 >= p > 0.1", "0.1 >= p > 0.05", "0.05 >= p"))

  pp <- ggplot(df, aes(x = names, y = p_hat)) +
    # baseline at chance agreement
    geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "grey60") +

    # cleaner CI + point
    geom_pointrange(aes(ymin = p_hat_lb, ymax = p_hat_ub),
                    color = "darkorchid2", linewidth = 0.7) +
    geom_point(color = "darkorchid2", size = 2.4) +

    figure_theme() +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Borzoi significance bins",
      y = "Proportion of tests\nwith concordant sign (%)"
    )

  return(pp)
}



organized_sim_results_dir = args[1]
visualization_dir = args[2]
borzoi_eqtl_effects_file = args[3]


df <- read.table(borzoi_eqtl_effects_file, header=TRUE, sep="\t")


# Load in data
#df <- read.table(results_summary_file, header=TRUE, sep="\t")
borzoi_bs <- get_bootstrapped_estimates(df$bs_borzoi_log_sed)  # N x B
B <- ncol(borzoi_bs)
p_plus  <- (1 + rowSums(borzoi_bs >= -0.01, na.rm = TRUE)) / (1 + B)
p_minus <- (1 + rowSums(borzoi_bs <= 0.01, na.rm = TRUE)) / (1 + B)
#p_plus  <- (1 + rowSums(borzoi_bs >= 0.0, na.rm = TRUE)) / (1 + B)
#p_minus <- (1 + rowSums(borzoi_bs <= 0.0, na.rm = TRUE)) / (1 + B)

borzoi_pvalues <- 2 * pmin(p_plus, p_minus)
df$borzoi_pvalue <- borzoi_pvalues
df$bs_borzoi_log_sed <- NULL

nan_indices <- is.na(df$beta_posterior) == FALSE
df <- df[nan_indices,]
borzoi_bs <- borzoi_bs[nan_indices,]




df$borzoi_mean = df$mean_borzoi_log_sed




borzoi_sign_plot = borzoi_sign_concordance_plot(df[df$maf < .05, ])
output_file <- paste0(visualization_dir, "borzoi_sign_concordance.pdf")
ggsave(borzoi_sign_plot, file=output_file,device = cairo_pdf, width=5.2, height=4.12, units="in")



power_file <- paste0(organized_sim_results_dir, "organized_power_sim_results.txt")
t1e_file <- paste0(organized_sim_results_dir, "organized_t1e_sim_results.txt")

power_df <- read.table(power_file, header=TRUE, sep="\t")
t1e_df <- read.table(t1e_file, header=TRUE, sep="\t")


####################
# Make power plot
####################

power_plot <- make_power_plot(power_df)
output_file <- paste0(visualization_dir, "power_plot.pdf")
ggsave(power_plot, file=output_file,device = cairo_pdf, width=7.2, height=4.12, units="in")


t1e_plot <- make_t1e_plot(t1e_df)
output_file <- paste0(visualization_dir, "t1e_plot.pdf")
ggsave(t1e_plot, file=output_file,device = cairo_pdf, width=7.2, height=4.12, units="in")



legender <- get_legend(power_plot + theme(legend.position="top"))
r01_plot <- plot_grid(legender, power_plot + theme(legend.position="none"), t1e_plot + theme(legend.position="none"),borzoi_sign_plot, ncol=1, labels=c(" ", "a", "b", "c"), label_y = 1.13,rel_heights=c(.1, 1, 1, .8))


output_file <- paste0(visualization_dir, "r01_plot.pdf")
ggsave(r01_plot, file=output_file,device = cairo_pdf, width=2.6, height=6.0, units="in")


