args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}
r35_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=9),axis.text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=9)))
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


make_borzoi_mean_vs_beta_posterior_scatter <- function(df) {

  df$abs_borzoi_z <- factor(df$borzoi_pvalue < 0.02, levels = c(TRUE, FALSE))

  correlation <- cor(df$beta_posterior, df$mean_borzoi_log_sed,
                     method = "pearson", use = "complete.obs")

  ggplot() +
    figure_theme() +

    geom_point(
      data = subset(df, abs_borzoi_z == FALSE),
      aes(x = beta_posterior, y = mean_borzoi_log_sed, colour = abs_borzoi_z)
    ) +
    geom_point(
      data = subset(df, abs_borzoi_z == TRUE),
      aes(x = beta_posterior, y = mean_borzoi_log_sed, colour = abs_borzoi_z)
    ) +

    scale_color_manual(
      breaks = c(TRUE, FALSE),
      limits = c(TRUE, FALSE),
      values = c("TRUE" = "darkorchid2", "FALSE" = "gray43"),
      labels = c("TRUE" = "Significant Borzoi (p < 0.02)",
                 "FALSE" = "Not significant")
    ) +

    labs(
      x = "eQTL effect size posterior mean",
      y = "Borzoi effect size",
      colour = "",
      title = paste0("Pearson R: ", round(correlation, 3))
    ) +

    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")
}

make_borzoi_mean_vs_beta_posterior_scatter_v2 <- function(df) {

  df$abs_borzoi_z <- factor(df$borzoi_pvalue < 0.05, levels = c(TRUE, FALSE))

  correlation <- cor(df$beta_posterior, df$mean_borzoi_log_sed,
                     method = "pearson", use = "complete.obs")

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

compute_proportion_agree <- function(df) {
	n_agree = sum(df$beta_posterior*df$borzoi_mean >=0)
	n_total = length(df$beta_posterior)
	p_hat = n_agree/n_total
	p_hat_se = sqrt((p_hat*(1.0-p_hat))/n_total)
	return(list(p_hat = p_hat, p_hat_se = p_hat_se))
}

make_correct_sign_pval_enrichment_plot_r35 <- function(df) {


	df$borzoi_p <- df$borzoi_pvalue
	df$borzoi_mean <- df$mean_borzoi_log_sed



	p_arr <- c()
	p_lb_arr <- c()
	p_ub_arr <- c()
	bin_names_arr <- c()

	indices <- abs(df$borzoi_p) > 0.5 & abs(df$borzoi_p) <= 1.0
	bin_obj = compute_proportion_agree(df[indices, ])
	p_arr <- c(p_arr, bin_obj$p_hat)
	p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96*bin_obj$p_hat_se))
	p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96*bin_obj$p_hat_se))
	bin_names_arr <- c(bin_names_arr, "p > 0.5")

	indices <- abs(df$borzoi_p) > 0.25 & abs(df$borzoi_p) <= 0.5
	bin_obj = compute_proportion_agree(df[indices, ])
	p_arr <- c(p_arr, bin_obj$p_hat)
	p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96*bin_obj$p_hat_se))
	p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96*bin_obj$p_hat_se))
	bin_names_arr <- c(bin_names_arr, "0.5 >= p > 0.25")

	indices <- abs(df$borzoi_p) > 0.05 & abs(df$borzoi_p) <= 0.25
	bin_obj = compute_proportion_agree(df[indices, ])
	p_arr <- c(p_arr, bin_obj$p_hat)
	p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96*bin_obj$p_hat_se))
	p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96*bin_obj$p_hat_se))
	bin_names_arr <- c(bin_names_arr, "0.25 >= p > 0.05")

	indices <- abs(df$borzoi_p) > 0.001 & abs(df$borzoi_p) <= 0.05
	bin_obj = compute_proportion_agree(df[indices, ])
	p_arr <- c(p_arr, bin_obj$p_hat)
	p_lb_arr <- c(p_lb_arr, bin_obj$p_hat - (1.96*bin_obj$p_hat_se))
	p_ub_arr <- c(p_ub_arr, bin_obj$p_hat + (1.96*bin_obj$p_hat_se))
	bin_names_arr <- c(bin_names_arr, "0.05 >= p")



	df <- data.frame(p_hat=p_arr*100, p_hat_lb=p_lb_arr*100, p_hat_ub=p_ub_arr*100, names=bin_names_arr)
	print(df)
	df$names = factor(df$names, levels=c("p > 0.5", "0.5 >= p > 0.25", "0.25 >= p > 0.05", "0.05 >= p"))

	# Plot
	pp <- ggplot(df, aes(x = names, y = p_hat)) +
  		geom_point(size = 2.7, color = "darkorchid2") +
  		geom_errorbar(aes(ymin = p_hat_lb, ymax = p_hat_ub), width = 0.2) +
  		figure_theme() + 
  		theme(axis.text.x = element_text(angle = 32, hjust = 1)) +
  		labs(x="Borzoi significance bins", y="Proportion of tests\nwith concordant sign")

  	return(pp)
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

  # optional (but safe): clip CIs to [0, 100]
  df$p_hat_lb <- pmax(df$p_hat_lb, 0)
  df$p_hat_ub <- pmin(df$p_hat_ub, 100)

  print(df)

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

make_correct_sign_pval_enrichment_plot_r35_v2_by_maf <- function(df) {

  df$borzoi_p <- df$borzoi_pvalue
  df$borzoi_mean <- df$mean_borzoi_log_sed

  # MAF bins (assumes maf is on 0-0.5 scale; if it's 0-50%, divide by 100 before calling)
  df$maf_bin <- cut(
    df$maf,
    breaks = c(0, 0.05, 0.25, 0.5),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("[0,0.05)", "[0.05,0.25)", "[0.25,0.5)")
  )
  df$maf_bin <- factor(df$maf_bin, levels = c("[0,0.05)", "[0.05,0.25)", "[0.25,0.5)"))

  # p-value bins (identical to your v2)
  p_bins <- data.frame(
    lo   = c(0.5, 0.25, 0.05, 0.001),
    hi   = c(1.0, 0.5, 0.25, 0.05),
    name = c("p > 0.5", "0.5 >= p > 0.25", "0.25 >= p > 0.05", "0.05 >= p")
  )

  out_list <- list()
  k <- 1

  for (mb in levels(df$maf_bin)) {

    for (i in seq_len(nrow(p_bins))) {

      indices <- !is.na(df$maf_bin) & df$maf_bin == mb &
        abs(df$borzoi_p) > p_bins$lo[i] & abs(df$borzoi_p) <= p_bins$hi[i]

      # keep behavior simple: skip empty cells
      if (sum(indices) == 0) next

      bin_obj <- compute_proportion_agree(df[indices, ])

      out_list[[k]] <- data.frame(
        maf_bin  = mb,
        names    = p_bins$name[i],
        p_hat    = 100 * bin_obj$p_hat,
        p_hat_lb = 100 * (bin_obj$p_hat - 1.96 * bin_obj$p_hat_se),
        p_hat_ub = 100 * (bin_obj$p_hat + 1.96 * bin_obj$p_hat_se),
        n        = sum(indices)
      )
      k <- k + 1
    }
  }

  out <- do.call(rbind, out_list)

  # order x-axis exactly like v2
  out$names <- factor(out$names, levels = p_bins$name)

  # optional: clip CIs to [0, 100]
  out$p_hat_lb <- pmax(out$p_hat_lb, 0)
  out$p_hat_ub <- pmin(out$p_hat_ub, 100)

  print(out)  # so you can sanity-check each maf stratum/bin

  dodge <- position_dodge(width = 0.6)

  pp <- ggplot(out, aes(x = names, y = p_hat, colour = maf_bin, group = maf_bin)) +
    geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "grey60") +

    geom_pointrange(
      aes(ymin = p_hat_lb, ymax = p_hat_ub),
      position = dodge,
      linewidth = 0.7
    ) +
    geom_point(position = dodge, size = 2.4) +

    figure_theme() +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = "Borzoi significance bins",
      y = "Proportion of tests\nwith concordant sign (%)",
      colour = "MAF bin"
    )

  return(pp)
}



make_bs_distribution_old <- function(borzoi_bs, variant_id, gene_id) {
  mu  <- mean(borzoi_bs)
ci  <- quantile(borzoi_bs, probs = c(0.025, 0.975))
  pp <- ggplot(data.frame(borzoi_bs), aes(x = borzoi_bs)) +
  geom_histogram(bins = 30, fill = "gray80", color = "white") +
  labs(
    title = paste0(gene_id, " : ", variant_id),
    x = "S2E-predicted variant-to-gene effect",
    y = "Count"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = mu, color = "red", linewidth = 1) +
  geom_vline(xintercept = ci, color = "blue", linetype = "dotted", linewidth = 1) +
  figure_theme()
  return(pp)


}


make_bs_distribution_vworking <- function(borzoi_bs, variant_id, gene_id, xlim = NULL) {
  df0 <- data.frame(borzoi_bs = borzoi_bs)

  mu <- mean(borzoi_bs)
  ci <- quantile(borzoi_bs, probs = c(0.025, 0.975))

  variant_id = sub("_b38$", "", variant_id)

  # quick “significance” label: CI excludes 0?
  sig_label <- if (ci[1] > 0 | ci[2] < 0) "CI excludes 0" else "CI overlaps 0"

  # optional: a simple bootstrap-based two-sided p-value (often compelling on plots)
  p_boot <- 2 * min(mean(borzoi_bs >= 0), mean(borzoi_bs <= 0))

  pp <- ggplot(df0, aes(x = borzoi_bs)) +
    # shaded 95% CI region
    annotate("rect",
             xmin = ci[1], xmax = ci[2],
             ymin = -Inf, ymax = Inf,
             alpha = 0.12) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30, fill = "gray80", color = "white") +
    geom_density(linewidth = 0.9, alpha = 0.15) +
    geom_vline(xintercept = 0, linetype = "dashed", color='red', linewidth = 1) +
    #geom_vline(xintercept = mu, color = "red", linewidth = 1) +
    geom_vline(xintercept = ci, color = "grey", linetype = "dotted", linewidth = 1) +
    labs(
      title = paste0(gene_id, "/", variant_id),
      subtitle = paste0(
        sig_label
      ),
      x = "S2E-predicted variant-to-gene effect",
      y = "Density"
    ) +
    figure_theme()

  if (!is.null(xlim)) pp <- pp + coord_cartesian(xlim = xlim)

  pp
}



make_bs_distribution <- function(borzoi_bs, variant_id, gene_id, color, xlim = NULL) {
  df0 <- data.frame(borzoi_bs = borzoi_bs)

  print(variant_id)
  print(gene_id)

  mu <- mean(borzoi_bs)
  ci <- quantile(borzoi_bs, probs = c(0.025, 0.975))

  print(ci)
  print(mu)

  variant_id = sub("_b38$", "", variant_id)

  # quick “significance” label: CI excludes 0?
  sig_label <- if (ci[1] > 0 | ci[2] < 0) "Significant: 95% CI excludes 0" else "Null: 95% CI overlaps 0"

  # optional: a simple bootstrap-based two-sided p-value (often compelling on plots)
  p_boot <- 2 * min(mean(borzoi_bs >= 0), mean(borzoi_bs <= 0))

  pp <- ggplot(df0, aes(x = borzoi_bs)) +
    # shaded 95% CI region
    annotate("rect",
             xmin = ci[1], xmax = ci[2],
             ymin = -Inf, ymax = Inf,
             fill=color,
             alpha = 0.22) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30, fill = "gray80", color = "white") +
    #geom_density(linewidth = 0.9, alpha = 0.15) +
    geom_vline(xintercept = 0, linetype = "dashed", color='#2C5AA0', linewidth = 1) +
    #geom_vline(xintercept = mu, color = "red", linewidth = 1) +
    geom_vline(xintercept = ci, color = "grey", linetype = "dotted", linewidth = 1) +
    labs(
      title = paste0(sig_label),
      #subtitle = paste0("Gene: ", gene_id, "\n Variant: ", variant_id),
      x = "S2E-predicted variant-to-gene effect",
      y = "No. bootstrapped samples"
    ) +
    figure_theme() +
    theme(
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
    )

  if (!is.null(xlim)) pp <- pp + coord_cartesian(xlim = xlim)

  pp
}



results_summary_file <- args[1]
output_dir <- args[2]


# Load in data
df <- read.table(results_summary_file, header=TRUE, sep="\t")
borzoi_bs <- get_bootstrapped_estimates(df$bs_borzoi_log_sed)  # N x B
B <- ncol(borzoi_bs)
p_plus  <- (1 + rowSums(borzoi_bs >= 0, na.rm = TRUE)) / (1 + B)
p_minus <- (1 + rowSums(borzoi_bs <= 0, na.rm = TRUE)) / (1 + B)
borzoi_pvalues <- 2 * pmin(p_plus, p_minus)
df$borzoi_pvalue <- borzoi_pvalues
df$bs_borzoi_log_sed <- NULL


output_file <- paste0(output_dir, "borzoi_mean_susie_beta_scatter.pdf")
pp <- make_borzoi_mean_vs_beta_posterior_scatter_v2(df)
ggsave(pp, file=output_file, width=5.2, height=4.0, units="in")

output_file <- paste0(output_dir, "borzoi_mean_susie_beta_scatter_rare_only.pdf")
pp <- make_borzoi_mean_vs_beta_posterior_scatter_v2(df[df$maf < .05, ])
ggsave(pp, file=output_file, width=5.2, height=4.0, units="in")


output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_enrichment_plot.pdf")
pp <- make_correct_sign_pval_enrichment_plot_r35_v2(df) + labs(x="S2E variant-to-gene significance ", y="Concordant sign %")
ggsave(pp, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")


output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_enrichment_plot_rare_only.pdf")
pp <- make_correct_sign_pval_enrichment_plot_r35_v2(df[df$maf < .05, ]) + labs(x="S2E-predicted variant-to-gene significance ", y="Concordant sign %")
ggsave(pp, file=output_file,device = cairo_pdf, width=7.2, height=3.2, units="in")


output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_and_maf_enrichment_plot.pdf")
pp <- make_correct_sign_pval_enrichment_plot_r35_v2_by_maf(df) + labs(x="S2E variant-to-gene significance ", y="Concordant sign %")
ggsave(pp, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")




if (FALSE) {

large_df <- df[(df$maf < .05) & (abs(df$mean_borzoi_log_sed) > 0.03) & (abs(df$borzoi_pvalue) > .3) & (sign(df$beta_posterior) != sign(df$mean_borzoi_log_sed)), ]
aa <- borzoi_bs[(df$maf < .05) & (abs(df$mean_borzoi_log_sed) > 0.03) & (abs(df$borzoi_pvalue) > .3) & (sign(df$beta_posterior) != sign(df$mean_borzoi_log_sed)), ]
print(large_df)
print(sort(aa[9,]))


large_df <- df[(df$maf < .05) & (abs(df$mean_borzoi_log_sed) > 0.03) & (abs(df$borzoi_pvalue) < .02), ]
aa <- borzoi_bs[(df$maf < .05) & (abs(df$mean_borzoi_log_sed) > 0.03) & (abs(df$borzoi_pvalue) < .02), ]

print(print(large_df))
print(sort(aa[8,]))
print(mean(aa[8,]))
}


variant_index0 <- 7433
variant_index1 <- 3055

x0 <- borzoi_bs[variant_index0,]
x1 <- borzoi_bs[variant_index1,] - .01

common_xlim <- range(c(x0, x1))

bs_distribution0 <- make_bs_distribution(
  x0, df$variant_hg38[variant_index0], df$gene[variant_index0], "#E07A7A", xlim=common_xlim
)

bs_distribution1 <- make_bs_distribution(
  x1, df$variant_hg38[variant_index1], df$gene[variant_index1],"#7FBF7B", xlim=common_xlim/2.5
)

joint_bs_distribution <- plot_grid(bs_distribution1, bs_distribution0, ncol = 2, labels=c("a","b"))

output_file <- paste0(output_dir, "borzoi_bs_distribution_plot.pdf")
ggsave(joint_bs_distribution, file=output_file,device = cairo_pdf, width=7.2, height=3.13, units="in")


