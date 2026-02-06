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

compute_proportion_agree <- function(df) {
	n_agree = sum(df$beta_posterior*df$borzoi_mean >=0)
	n_total = length(df$beta_posterior)
	p_hat = n_agree/n_total
	p_hat_se = sqrt((p_hat*(1.0-p_hat))/n_total)
	return(list(p_hat = p_hat, p_hat_se = p_hat_se))
}



compute_power <- function(df, p_thresh) {
  n_agree = sum(df$borzoi_p <= p_thresh)
  n_total = length(df$borzoi_p)
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
	#print(df)
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


make_correct_sign_pval_enrichment_plot_cmp_pvalue_with_magnitude <- function(df) {

print(sum(df$borzoi_pvalue < .2))

# dat = your original data.frame of variant-gene(-tissue) pairs
dat <- df   # only if df is truly your data.frame; otherwise set dat <- <your_data>

# sanity check
stopifnot(is.data.frame(dat))

dat$borzoi_p    <- dat$borzoi_pvalue
dat$borzoi_mean <- dat$mean_borzoi_log_sed

cutoffs <- c(100, 250, 500, 1000)

.compute_topN_row <- function(dat, ord_idx, N, label, group_name) {
  N_use <- min(N, length(ord_idx))
  idx   <- ord_idx[seq_len(N_use)]
  bin_obj <- compute_proportion_agree(dat[idx, , drop = FALSE])

  p_hat <- bin_obj$p_hat
  lb    <- p_hat - 1.96 * bin_obj$p_hat_se
  ub    <- p_hat + 1.96 * bin_obj$p_hat_se

  data.frame(
    group    = group_name,
    bin      = label,
    N        = N_use,
    p_hat    = 100 * p_hat,
    p_hat_lb = 100 * pmax(lb, 0),
    p_hat_ub = 100 * pmin(ub, 1),
    stringsAsFactors = FALSE
  )
}

ord_sig <- order(dat$borzoi_p, decreasing = FALSE, na.last = NA)
ord_mag <- order(abs(dat$borzoi_mean), decreasing = TRUE, na.last = NA)

sumdf <- do.call(
  rbind,
  lapply(cutoffs, function(N) {
    rbind(
      .compute_topN_row(dat, ord_sig, N, paste0("Top ", N), "Most significant (lowest p)"),
      .compute_topN_row(dat, ord_mag, N, paste0("Top ", N), "Largest |borzoi effect|")
    )
  })
)

sumdf$bin <- factor(sumdf$bin, levels = paste0("Top ", cutoffs))
sumdf$group <- factor(sumdf$group, levels = c("Most significant (lowest p)", "Largest |borzoi effect|"))

pp <- ggplot(sumdf, aes(x = bin, y = p_hat, color = group)) +
  #geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "grey60") +
  geom_pointrange(aes(ymin = p_hat_lb, ymax = p_hat_ub),
                  position = position_dodge(width = 0.45),
                  linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.45), size = 2.4) +
  scale_color_manual(
    values = c(
      "Most significant (lowest p)" = "darkorchid2",
      "Largest |borzoi effect|"     = "#d62728"
    ),
    breaks = c("Most significant (lowest p)", "Largest |borzoi effect|")  # keep order
  ) +
  figure_theme() +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom"
  ) +
  labs(
    x = "Ranking bins (top N)",
    y = "Concordant sign (%)",
    color = NULL
  )

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

  #print(out)  # so you can sanity-check each maf stratum/bin

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



make_power_plot_controlled_for_abs_distance <- function(df) {

  df$borzoi_p <- df$borzoi_pvalue
  df$borzoi_mean <- df$mean_borzoi_log_sed

  df$abs_borzoi_mean <- abs(df$mean_borzoi_log_sed)

  df$abs_dist = abs(df$dist_to_tss)





  # MAF bins (assumes maf is on 0-0.5 scale; if it's 0-50%, divide by 100 before calling)
  df$distance_bin <- cut(
    df$abs_dist,
    breaks = c(0.0, 1000.0, 5000.0, 20000.0, 100000.0),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("[0,1)", "[1,5)", "[5,20)", "[20,100)")
  )
  df$distance_bin <- factor(df$distance_bin, levels = c("[0,1)", "[1,5)", "[5,20)", "[20,100)"))

  # p-value bins (identical to your v2)
  p_bins <- data.frame(
    lo   = c(0.0000001, 0.000001),
    hi   = c(0.2, 0.05),
    name = c( "p <= 0.2", "p <= 0.05")
  )




  out_list <- list()
  k <- 1

  for (mb in levels(df$distance_bin)) {

    for (i in seq_len(nrow(p_bins))) {



      indices <- !is.na(df$distance_bin) & df$distance_bin == mb




      # keep behavior simple: skip empty cells
      if (sum(indices) == 0) next

      bin_obj <- compute_power(df[indices, ], p_bins$hi[i])

      out_list[[k]] <- data.frame(
        distance_bin  = mb,
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

  #print(out)  # so you can sanity-check each maf stratum/bin

  dodge <- position_dodge(width = 0.6)

  out$distance_bin <- factor(out$distance_bin, levels = c("[0,1)", "[1,5)", "[5,20)", "[20,100)"))


  pp <- ggplot(out, aes(x = distance_bin, y = p_hat, colour = names, group = names)) +

    geom_pointrange(
      aes(ymin = p_hat_lb, ymax = p_hat_ub),
      position = dodge,
      linewidth = 0.7
    ) +
    geom_point(position = dodge, size = 2.4) +


  scale_colour_manual(
    values = c(
      "p <= 0.2" = "darkorchid4",
      "p <= 0.05"   = "darkorchid1"
    )
  ) +


    figure_theme() +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = "Absolute distance to TSS (kb)",
      y = "Power (%)",
      colour = "Significance"
    ) 

  return(pp)
}




make_power_plot_controlled_for_distance <- function(df) {

  df$borzoi_p <- df$borzoi_pvalue
  df$borzoi_mean <- df$mean_borzoi_log_sed

  df$abs_borzoi_mean <- abs(df$mean_borzoi_log_sed)






  # MAF bins (assumes maf is on 0-0.5 scale; if it's 0-50%, divide by 100 before calling)
  df$distance_bin <- cut(
    df$dist_to_tss,
    breaks = c(-100000.0, -20000.0, -5000.0, -1000.0, 0.0, 1000.0, 5000.0, 20000.0, 100000.0),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("[-100,-20)", "[-20,-5)", "[-5,-1)", "[-1,0)", "[0,1)", "[1,5)", "[5,20)", "[20,100)")
  )
  df$distance_bin <- factor(df$distance_bin, levels = c("[-100,-20)", "[-20,-5)", "[-5,-1)", "[-1,0)", "[0,1)", "[1,5)", "[5,20)", "[20,100)"))

  # p-value bins (identical to your v2)
  p_bins <- data.frame(
    lo   = c(0.0000001, 0.000001),
    hi   = c(0.2, 0.05),
    name = c( "p <= 0.2", "p <= 0.05")
  )




  out_list <- list()
  k <- 1

  for (mb in levels(df$distance_bin)) {

    for (i in seq_len(nrow(p_bins))) {



      indices <- !is.na(df$distance_bin) & df$distance_bin == mb




      # keep behavior simple: skip empty cells
      if (sum(indices) == 0) next

      bin_obj <- compute_power(df[indices, ], p_bins$hi[i])

      out_list[[k]] <- data.frame(
        distance_bin  = mb,
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

  #print(out)  # so you can sanity-check each maf stratum/bin

  dodge <- position_dodge(width = 0.6)

  out$distance_bin <- factor(out$distance_bin, levels = c("[-100,-20)", "[-20,-5)", "[-5,-1)", "[-1,0)", "[0,1)", "[1,5)", "[5,20)", "[20,100)"))


  pp <- ggplot(out, aes(x = distance_bin, y = p_hat, colour = names, group = names)) +

    geom_pointrange(
      aes(ymin = p_hat_lb, ymax = p_hat_ub),
      position = dodge,
      linewidth = 0.7
    ) +
    geom_point(position = dodge, size = 2.4) +


  scale_colour_manual(
    values = c(
      "p <= 0.2" = "darkorchid4",
      "p <= 0.05"   = "darkorchid1"
    )
  ) +


    figure_theme() +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = "Distance to TSS (kb)",
      y = "Power (%)",
      colour = "Significance"
    ) + 
    geom_vline(
    xintercept = 4.5,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  )

  return(pp)
}





make_correct_sign_pval_enrichment_plot_controlled_for_distance <- function(df) {

  df$borzoi_p <- df$borzoi_pvalue
  df$borzoi_mean <- df$mean_borzoi_log_sed

  df$abs_borzoi_mean <- abs(df$mean_borzoi_log_sed)


  # MAF bins (assumes maf is on 0-0.5 scale; if it's 0-50%, divide by 100 before calling)
  df$distance_bin <- cut(
    df$dist_to_tss,
    breaks = c(-100000.0, -20000.0, -5000.0, -1000.0, 0.0, 1000.0, 5000.0, 20000.0, 100000.0),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("[-100,-20)", "[-20,-5)", "[-5,-1)", "[-1,0)", "[0,1)", "[1,5)", "[5,20)", "[20,100)")
  )
  df$distance_bin <- factor(df$distance_bin, levels = c("[-100,-20)", "[-20,-5)", "[-5,-1)", "[-1,0)", "[0,1)", "[1,5)", "[5,20)", "[20,100)"))

  # p-value bins (identical to your v2)
  p_bins <- data.frame(
    lo   = c(0.5, 0.05, 0.001),
    hi   = c(1.0,0.5, 0.05),
    name = c("p > 0.5", "0.5 >= p > 0.05", "0.05 >= p")
  )




  out_list <- list()
  k <- 1

  for (mb in levels(df$distance_bin)) {

    for (i in seq_len(nrow(p_bins))) {



      indices <- !is.na(df$distance_bin) & df$distance_bin == mb &
        abs(df$borzoi_p) > p_bins$lo[i] & abs(df$borzoi_p) <= p_bins$hi[i]




      # keep behavior simple: skip empty cells
      if (sum(indices) == 0) next

      bin_obj <- compute_proportion_agree(df[indices, ])

      out_list[[k]] <- data.frame(
        distance_bin  = mb,
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

  out$distance_bin <- factor(out$distance_bin, levels = c("[-100,-20)", "[-20,-5)", "[-5,-1)", "[-1,0)", "[0,1)", "[1,5)", "[5,20)", "[20,100)"))



  # order x-axis exactly like v2
  out$names <- factor(out$names, levels = p_bins$name)

  # optional: clip CIs to [0, 100]
  out$p_hat_lb <- pmax(out$p_hat_lb, 0)
  out$p_hat_ub <- pmin(out$p_hat_ub, 100)

  #print(out)  # so you can sanity-check each maf stratum/bin

  dodge <- position_dodge(width = 0.6)


  pp <- ggplot(out, aes(x = distance_bin, y = p_hat, colour = names, group = names)) +
    geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "grey60") +

    geom_pointrange(
      aes(ymin = p_hat_lb, ymax = p_hat_ub),
      position = dodge,
      linewidth = 0.7
    ) +
    geom_point(position = dodge, size = 2.4) +

  scale_colour_manual(
    values = c(
      "p > 0.5" = "grey45",
      "0.5 >= p > 0.05" = "darkorchid4",
      "0.05 >= p"   = "darkorchid1"
    )
  ) +


    figure_theme() +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = "Absolute distance to TSS (kb)",
      y = "Concordant sign %",
      colour = "Significance"
    ) +
    coord_cartesian(ylim = c(48.0, NA)) +
    geom_vline(
    xintercept = 4.5,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  )
  return(pp)
}




make_correct_sign_pval_enrichment_plot_controlled_for_abs_distance <- function(df) {

  df$borzoi_p <- df$borzoi_pvalue
  df$borzoi_mean <- df$mean_borzoi_log_sed

  df$abs_borzoi_mean <- abs(df$mean_borzoi_log_sed)

  df$abs_dist = abs(df$dist_to_tss)

  # MAF bins (assumes maf is on 0-0.5 scale; if it's 0-50%, divide by 100 before calling)
  df$distance_bin <- cut(
    df$abs_dist,
    breaks = c(0.0, 1000.0, 5000.0, 20000.0, 100000.0),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("[0,1)", "[1,5)", "[5,20)", "[20,100)")
  )
  df$distance_bin <- factor(df$distance_bin, levels = c("[0,1)", "[1,5)", "[5,20)", "[20,100)"))

  # p-value bins (identical to your v2)
  p_bins <- data.frame(
    lo   = c(0.5, 0.05, 0.001),
    hi   = c(1.0,0.5, 0.05),
    name = c("p > 0.5", "0.5 >= p > 0.05", "0.05 >= p")
  )




  out_list <- list()
  k <- 1

  for (mb in levels(df$distance_bin)) {

    for (i in seq_len(nrow(p_bins))) {



      indices <- !is.na(df$distance_bin) & df$distance_bin == mb &
        abs(df$borzoi_p) > p_bins$lo[i] & abs(df$borzoi_p) <= p_bins$hi[i]




      # keep behavior simple: skip empty cells
      if (sum(indices) == 0) next

      bin_obj <- compute_proportion_agree(df[indices, ])

      out_list[[k]] <- data.frame(
        distance_bin  = mb,
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

  out$distance_bin <- factor(out$distance_bin, levels = c("[0,1)", "[1,5)", "[5,20)", "[20,100)"))



  # order x-axis exactly like v2
  out$names <- factor(out$names, levels = p_bins$name)

  # optional: clip CIs to [0, 100]
  out$p_hat_lb <- pmax(out$p_hat_lb, 0)
  out$p_hat_ub <- pmin(out$p_hat_ub, 100)

  #print(out)  # so you can sanity-check each maf stratum/bin

  dodge <- position_dodge(width = 0.6)


  pp <- ggplot(out, aes(x = distance_bin, y = p_hat, colour = names, group = names)) +
    geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "grey60") +

    geom_pointrange(
      aes(ymin = p_hat_lb, ymax = p_hat_ub),
      position = dodge,
      linewidth = 0.7
    ) +
    geom_point(position = dodge, size = 2.4) +

  scale_colour_manual(
    values = c(
      "p > 0.5" = "grey45",
      "0.5 >= p > 0.05" = "darkorchid4",
      "0.05 >= p"   = "darkorchid1"
    )
  ) +


    figure_theme() +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = "Absolute distance to TSS (kb)",
      y = "Concordant sign %",
      colour = "Significance"
    ) +
    coord_cartesian(ylim = c(48.0, NA))

  return(pp)
}



make_correct_sign_pval_enrichment_plot_controlled_for_magnitude <- function(df) {

  df$borzoi_p <- df$borzoi_pvalue
  df$borzoi_mean <- df$mean_borzoi_log_sed

  df$abs_borzoi_mean <- abs(df$mean_borzoi_log_sed)

  # MAF bins (assumes maf is on 0-0.5 scale; if it's 0-50%, divide by 100 before calling)
  df$magnitude_bin <- cut(
    df$abs_borzoi_mean,
    breaks = c(.01, 0.025, .05, .1, .2, .3, .4, .5),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("[0.01,0.025)", "[0.025,0.05)", "[0.05,0.2)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)", "[0.4,0.5)")
  )
  df$magnitude_bin <- factor(df$magnitude_bin, levels = c("[0.01,0.025)", "[0.025,0.05)", "[0.05,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)", "[0.4,0.5)"))

  # p-value bins (identical to your v2)
  p_bins <- data.frame(
    lo   = c(0.5, 0.05, 0.001),
    hi   = c(1.0,0.5, 0.05),
    name = c("p > 0.5", "0.5 >= p > 0.05", "0.05 >= p")
  )




  out_list <- list()
  k <- 1

  for (mb in levels(df$magnitude_bin)) {

    for (i in seq_len(nrow(p_bins))) {



      indices <- !is.na(df$magnitude_bin) & df$magnitude_bin == mb &
        abs(df$borzoi_p) > p_bins$lo[i] & abs(df$borzoi_p) <= p_bins$hi[i]




      # keep behavior simple: skip empty cells
      if (sum(indices) == 0) next

      bin_obj <- compute_proportion_agree(df[indices, ])

      out_list[[k]] <- data.frame(
        magnitude_bin  = mb,
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

  #print(out)  # so you can sanity-check each maf stratum/bin

  dodge <- position_dodge(width = 0.6)


  pp <- ggplot(out, aes(x = magnitude_bin, y = p_hat, colour = names, group = names)) +
    geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "grey60") +

    geom_pointrange(
      aes(ymin = p_hat_lb, ymax = p_hat_ub),
      position = dodge,
      linewidth = 0.7
    ) +
    geom_point(position = dodge, size = 2.4) +

  scale_colour_manual(
    values = c(
      "p > 0.5" = "grey45",
      "0.5 >= p > 0.05" = "darkorchid4",
      "0.05 >= p"   = "darkorchid1"
    )
  ) +


    figure_theme() +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = "Borzoi effect magnitude",
      y = "Concordant sign %",
      colour = "Significance"
    ) +
    coord_cartesian(ylim = c(48.0, NA))

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

make_mean_vs_stdev_borzoi_effect_scatter <- function(df) {
  df$abs_mean_borzoi = abs(df$mean_borzoi_log_sed)
  pp <-ggplot(data = df, aes(x = abs_mean_borzoi, y = sdev_borzoi_log_sed)) +
    geom_point(color="slateblue", alpha=.5,size=1.0) + 
    labs(x="|Borzoi effect|", y="SD(Borzoi effect)") +
    geom_abline(intercept = 0, slope = 1) +    # y = x
    #geom_abline(intercept = 0, slope = -1) +   # y = -x 
    figure_theme()

  return(pp)
}

make_mean_vs_neg_log10_borzoi_effect_scatter <- function(df) {
  df$abs_mean_borzoi = abs(df$mean_borzoi_log_sed)
  df$neglog10_borzoi_pvalues = -log10(df$borzoi_pvalue)
  pp <-ggplot(data = df, aes(x = abs_mean_borzoi, y = neglog10_borzoi_pvalues)) +
    geom_point(color="slateblue", alpha=.5,size=1.0) + 
    labs(x="|Borzoi effect|", y="-log10(Borzoi pvalue)") +
    #geom_abline(intercept = 0, slope = 1) +    # y = x
    #geom_abline(intercept = 0, slope = -1) +   # y = -x 
    figure_theme()

  return(pp)
}

sign_consistency_regression <- function(df) {
  df$same_sign <- as.integer(
    sign(df$mean_borzoi_log_sed) == sign(df$beta_posterior)
  )

  fit <- glm(
    #same_sign ~ abs(mean_borzoi_log_sed) + sqrt(abs(mean_borzoi_log_sed)) + I(mean_borzoi_log_sed^2) + I(-log10(borzoi_pvalue)),
    same_sign ~  abs(mean_borzoi_log_sed) + I(-log10(borzoi_pvalue)),
    data = df,
    family = binomial(link = "logit")
  )

  summary(fit)

}

rank_compare_perm_test <- function(
  dat,
  cutoffs = c(100, 250, 500, 1000),
  B = 2000,
  seed = 1,
  keep_perm_dists = FALSE
) {
  stopifnot(is.data.frame(dat))
  stopifnot(all(c("borzoi_pvalue", "mean_borzoi_log_sed") %in% colnames(dat)))

  set.seed(seed)

  # local copies (avoid name collisions)
  dat <- dat
  dat$borzoi_p    <- dat$borzoi_pvalue
  dat$borzoi_mean <- dat$mean_borzoi_log_sed

  # helper: compute p_hat on Top-N by two orderings, and their diff
  compute_topN_stats <- function(d) {
    ord_sig <- order(d$borzoi_p, decreasing = FALSE, na.last = NA)                # smallest p first
    ord_mag <- order(abs(d$borzoi_mean), decreasing = TRUE, na.last = NA)         # largest |effect| first

    out_list <- vector("list", length(cutoffs))

    for (i in seq_along(cutoffs)) {
      N <- cutoffs[i]

      Nsig <- min(N, length(ord_sig))
      Nmag <- min(N, length(ord_mag))

      idx_sig <- ord_sig[seq_len(Nsig)]
      idx_mag <- ord_mag[seq_len(Nmag)]

      p_sig <- compute_proportion_agree(d[idx_sig, , drop = FALSE])$p_hat
      p_mag <- compute_proportion_agree(d[idx_mag, , drop = FALSE])$p_hat

      out_list[[i]] <- list(
        N = N,
        Nsig = Nsig,
        Nmag = Nmag,
        p_sig = p_sig,
        p_mag = p_mag,
        diff = p_sig - p_mag,
        overlap = length(intersect(idx_sig, idx_mag))
      )
    }

    out_list
  }

  # observed stats
  obs_list <- compute_topN_stats(dat)

  # storage for permutation distributions (one vector per cutoff)
  perm_mat <- matrix(NA_real_, nrow = B, ncol = length(cutoffs))

  # permutation loop: shuffle rows to break association between (p, effect) and concordance structure
  for (b in seq_len(B)) {
    perm_dat <- dat[sample.int(nrow(dat)), , drop = FALSE]
    perm_list <- compute_topN_stats(perm_dat)
    perm_mat[b, ] <- vapply(perm_list, `[[`, numeric(1), "diff")
  }

  # summarize results
  results <- do.call(
    rbind,
    lapply(seq_along(cutoffs), function(i) {
      N <- cutoffs[i]
      obs <- obs_list[[i]]

      perm_dist <- perm_mat[, i]
      pval <- mean(abs(perm_dist) >= abs(obs$diff))  # two-sided

      # permutation CI (central 95%)
      ci <- as.numeric(stats::quantile(perm_dist, probs = c(0.025, 0.975), na.rm = TRUE))

      data.frame(
        N = N,
        Nsig = obs$Nsig,
        Nmag = obs$Nmag,
        p_sig = obs$p_sig,
        p_mag = obs$p_mag,
        diff_sig_minus_mag = obs$diff,
        diff_sig_minus_mag_pct = 100 * obs$diff,
        p_value = pval,
        perm_ci_low = ci[1],
        perm_ci_high = ci[2],
        perm_ci_low_pct = 100 * ci[1],
        perm_ci_high_pct = 100 * ci[2],
        stringsAsFactors = FALSE
      )
    })
  )

  overlap <- do.call(
    rbind,
    lapply(obs_list, function(obs) {
      data.frame(
        N = obs$N,
        overlap = obs$overlap,
        overlap_frac_of_N = obs$overlap / min(obs$Nsig, obs$Nmag),
        stringsAsFactors = FALSE
      )
    })
  )

  out <- list(
    results = results,
    overlap = overlap,
    B = B,
    seed = seed,
    cutoffs = cutoffs
  )

  if (keep_perm_dists) {
    out$perm_dists <- perm_mat
  }

  out
}



results_summary_file <- args[1]
output_dir <- args[2]


# Load in data
df <- read.table(results_summary_file, header=TRUE, sep="\t")
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

indices = (df$borzoi_pvalue < .025) & (((df$beta_posterior)*(df$mean_borzoi_log_sed)) < 0.0)

print(df[indices,])



sub_bs <- borzoi_bs[indices,]

cbind(
  min = apply(sub_bs, 1, min, na.rm = TRUE),
  max = apply(sub_bs, 1, max, na.rm = TRUE)
)


sig_indices = df$borzoi_pvalue < .2
obj=compute_proportion_agree(df[sig_indices,])
print(obj$p_hat)

sig_indices = df$borzoi_pvalue < .1
obj=compute_proportion_agree(df[sig_indices,])
print(obj$p_hat)

sig_indices = df$borzoi_pvalue < .05
obj=compute_proportion_agree(df[sig_indices,])
print(obj$p_hat)

sig_indices = df$borzoi_pvalue < .025
print(sum(sig_indices))
obj=compute_proportion_agree(df[sig_indices,])
print(obj$p_hat)

if (FALSE) {
sig_indices = abs(df$mean_borzoi_log_sed) > .24
print(sum(sig_indices))
obj=compute_proportion_agree(df[sig_indices,])
print(obj$p_hat)
}
print("DONE")







print("DONE")

output_file <- paste0(output_dir, "borzoi_power_to_detect_by_abs_dist_to_tss_enrichment_plot.pdf")
pp <- make_power_plot_controlled_for_abs_distance(df) 
ggsave(pp, file=output_file,device = cairo_pdf, width=7.2, height=4.12, units="in")


output_file <- paste0(output_dir, "borzoi_power_to_detect_by_dist_to_tss_enrichment_plot.pdf")
pp <- make_power_plot_controlled_for_distance(df) 
ggsave(pp, file=output_file,device = cairo_pdf, width=7.2, height=4.12, units="in")



output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_control_for__absdist_to_tss_enrichment_plot.pdf")
pp <- make_correct_sign_pval_enrichment_plot_controlled_for_abs_distance(df) 
ggsave(pp, file=output_file,device = cairo_pdf, width=7.2, height=4.12, units="in")


output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_control_for_dist_to_tss_enrichment_plot.pdf")
pp <- make_correct_sign_pval_enrichment_plot_controlled_for_distance(df) 
ggsave(pp, file=output_file,device = cairo_pdf, width=7.2, height=4.12, units="in")



output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_control_for_magnitude_enrichment_plot.pdf")
pp <- make_correct_sign_pval_enrichment_plot_controlled_for_magnitude(df)
ggsave(pp, file=output_file,device = cairo_pdf, width=7.2, height=4.12, units="in")


output_file <- paste0(output_dir, "borzoi_mean_susie_beta_scatter.pdf")
pp <- make_borzoi_mean_vs_beta_posterior_scatter_v2(df, p=.05)
ggsave(pp, file=output_file, width=5.2, height=4.0, units="in")

output_file <- paste0(output_dir, "borzoi_mean_susie_beta_all_snps_scatter.pdf")
pp <- make_borzoi_mean_vs_beta_posterior_scatter_v2(df,p=0.0)
ggsave(pp, file=output_file, width=5.2, height=4.0, units="in")

if (FALSE) {
output_file <- paste0(output_dir, "borzoi_mean_susie_beta_scatter_rare_only.pdf")
pp <- make_borzoi_mean_vs_beta_posterior_scatter_v2(df[df$maf < .05, ])
ggsave(pp, file=output_file, width=5.2, height=4.0, units="in")
}

output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_enrichment_plot.pdf")
pp <- make_correct_sign_pval_enrichment_plot_r35_v2(df) + labs(x="S2E variant-to-gene significance ", y="Concordant sign %")
ggsave(pp, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")


if (FALSE) {
output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_enrichment_plot_rare_only.pdf")
pp <- make_correct_sign_pval_enrichment_plot_r35_v2(df[df$maf < .05, ]) + labs(x="S2E-predicted variant-to-gene significance ", y="Concordant sign %")
ggsave(pp, file=output_file,device = cairo_pdf, width=7.2, height=3.2, units="in")
}


output_file <- paste0(output_dir, "borzoi_correct_sign_pval_significance_and_maf_enrichment_plot.pdf")
pp <- make_correct_sign_pval_enrichment_plot_r35_v2_by_maf(df) + labs(x="S2E variant-to-gene significance ", y="Concordant sign %")
ggsave(pp, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")



output_file <- paste0(output_dir, "abs_borzoi_mean_borzoi_vs_stdev_borzoi.pdf")
pp <- make_mean_vs_stdev_borzoi_effect_scatter(df) 
ggsave(pp, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")

output_file <- paste0(output_dir, "abs_borzoi_mean_borzoi_vs_neglog10p.pdf")
pp <- make_mean_vs_neg_log10_borzoi_effect_scatter(df) 
ggsave(pp, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")


if (FALSE) {
output_file <- paste0(output_dir, "borzoi_correct_sign_enrichment_plot_cmp_pval_with_magnitude.pdf")
pp <- make_correct_sign_pval_enrichment_plot_cmp_pvalue_with_magnitude(df) 
ggsave(pp, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")

# Get significance on this
rank_sig_obj = rank_compare_perm_test(df)
}




if (FALSE) {
variant_index0 <- 259
x0 <- borzoi_bs[variant_index0,]


bs_distribution0 <- make_bs_distribution(
  x0, df$variant_hg38[variant_index0], df$gene[variant_index0], "White"
)

output_file <- paste0(output_dir, "borzoi_bs_distribution_plot2.pdf")
ggsave(bs_distribution0 + labs(title=paste0(df$variant_hg38[variant_index0], " : ", df$gene[variant_index0], " : Whole Blood")) , file=output_file,device = cairo_pdf, width=7.2, height=2.75, units="in")
}

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


if (FALSE) {
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
}

