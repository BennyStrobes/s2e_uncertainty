args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}







make_expression_correlation_histogram <- function(df, color, gene_subset_string) {
	pp <- ggplot(df, aes(x = expression_correlation)) +
  geom_histogram(bins = 30,fill = color, color = "white") +
  labs(x = "Expression correlation",y = "No. genes", title=gene_subset_string) +
  figure_theme() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "steelblue4") +
   theme(plot.title = element_text(hjust = 0.5))

  return(pp)
}


make_mean_correlation_standard_error_plot <- function(df) {

  sorted_p <- sort(df$corr_distribution_p)
  #thresh = sorted_p[100]
  #print(sum(df$corr_distribution_p <= thresh))

  mean_arr <- c()
  mean_lb_arr <- c()
  mean_ub_arr <- c()
  n_gene_arr <- c()
  bin_names_arr <- c()

  indices <- abs(df$corr_distribution_p) <= 1.0
  corrs = df$raw_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)

  NN = 5000
  thresh = sorted_p[NN]
  indices <- abs(df$corr_distribution_p) <= thresh
  corrs = df$raw_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("N=", NN)
  bin_names_arr <- c(bin_names_arr, namer)


  NN = 2000
  thresh = sorted_p[NN]
  indices <- abs(df$corr_distribution_p) <= thresh
  corrs = df$raw_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("N=", NN)
  bin_names_arr <- c(bin_names_arr, namer)


  NN = 1000
  thresh = sorted_p[NN]
  indices <- abs(df$corr_distribution_p) <= thresh
  corrs = df$raw_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("N=", NN)
  bin_names_arr <- c(bin_names_arr, namer)

  NN = 500
  thresh = sorted_p[NN]
  indices <- abs(df$corr_distribution_p) <= thresh
  corrs = df$raw_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("N=", NN)
  bin_names_arr <- c(bin_names_arr, namer)

  NN = 200
  thresh = sorted_p[NN]
  indices <- abs(df$corr_distribution_p) <= thresh
  corrs = df$raw_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("N=", NN)
  bin_names_arr <- c(bin_names_arr, namer)

  NN = 100
  thresh = sorted_p[NN]
  indices <- abs(df$corr_distribution_p) <= thresh
  corrs = df$raw_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("N=", NN)
  bin_names_arr <- c(bin_names_arr, namer)



  df2 <- data.frame(
    mean_correlation    = mean_arr,
    mean_correlation_lb = mean_lb_arr,
    mean_correlation_ub = mean_ub_arr,
    names    = bin_names_arr,
    n_genes = n_gene_arr
  )


  df2$names = factor(df2$names, levels=as.character(df2$names))

  print(df2)
  pp <- ggplot(df2, aes(x = names, y = mean_correlation)) +
    # baseline at chance agreement
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey60") +

    # cleaner CI + point
    geom_pointrange(aes(ymin = mean_correlation_lb, ymax = mean_correlation_ub),
                    color = "darkorchid2", linewidth = 0.7) +
    geom_point(color = "darkorchid2", size = 2.4) +

    figure_theme() +
    labs(
      x = "Top N genes",
      y = "Average expression correlation"
    )

    return(pp)
}


make_directional_fsr_calibration_plot <- function(df, pvalue_threshold=1.0) {
  point_color <- "#007C89"
  reference_color <- "#3F3F46"

  df <- df[!is.na(df$directional_FSR) &
             !is.na(df$raw_correlation) &
             !is.na(df$raw_correlation_pvalue) &
             df$raw_correlation_pvalue <= pvalue_threshold, ]

  bin_breaks <- c(0, 0.1, 0.2, 0.3, 0.4, Inf)
  bin_labels <- c("[0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", ">0.4")

  df$directional_FSR_bin <- cut(
    df$directional_FSR,
    breaks=bin_breaks,
    labels=bin_labels,
    include.lowest=TRUE,
    right=TRUE
  )

  observed_fsr_arr <- c()
  observed_fsr_lb_arr <- c()
  observed_fsr_ub_arr <- c()
  expected_fsr_arr <- c()
  n_gene_arr <- c()
  bin_arr <- c()
  x_arr <- c()

  for (bin_iter in 1:length(bin_labels)) {
    bin_name <- bin_labels[bin_iter]
    bin_df <- df[df$directional_FSR_bin == bin_name, ]
    n_genes <- nrow(bin_df)

    if (n_genes == 0) {
      observed_fsr <- NA
      observed_fsr_lb <- NA
      observed_fsr_ub <- NA
      expected_fsr <- NA
    } else {
      n_false_sign <- sum(bin_df$raw_correlation < 0)
      observed_fsr <- n_false_sign/n_genes
      observed_fsr_se <- sqrt(observed_fsr*(1.0 - observed_fsr)/n_genes)
      observed_fsr_lb <- max(0.0, observed_fsr - 1.96*observed_fsr_se)
      observed_fsr_ub <- min(1.0, observed_fsr + 1.96*observed_fsr_se)
      expected_fsr <- mean(bin_df$directional_FSR)
    }

    observed_fsr_arr <- c(observed_fsr_arr, observed_fsr)
    observed_fsr_lb_arr <- c(observed_fsr_lb_arr, observed_fsr_lb)
    observed_fsr_ub_arr <- c(observed_fsr_ub_arr, observed_fsr_ub)
    expected_fsr_arr <- c(expected_fsr_arr, expected_fsr)
    n_gene_arr <- c(n_gene_arr, n_genes)
    bin_arr <- c(bin_arr, bin_name)
    x_arr <- c(x_arr, bin_iter)
  }

  calibration_df <- data.frame(
    directional_FSR_bin=bin_arr,
    x=x_arr,
    observed_fsr=observed_fsr_arr,
    observed_fsr_lb=observed_fsr_lb_arr,
    observed_fsr_ub=observed_fsr_ub_arr,
    expected_fsr=expected_fsr_arr,
    n_genes=n_gene_arr
  )



  calibration_df$directional_FSR_bin <- factor(calibration_df$directional_FSR_bin, levels=bin_labels)
  calibration_df$bin_label <- paste0(calibration_df$directional_FSR_bin, "\nN=", calibration_df$n_genes)
  print(calibration_df)

  pp <- ggplot(calibration_df, aes(x=x, y=observed_fsr)) +
    geom_segment(
      aes(x=x - 0.35, xend=x + 0.35, y=expected_fsr, yend=expected_fsr),
      color=reference_color,
      linetype="dashed",
      linewidth=0.6,
      na.rm=TRUE
    ) +
    geom_line(color=point_color, linewidth=0.7, na.rm=TRUE) +
    geom_pointrange(
      aes(ymin=observed_fsr_lb, ymax=observed_fsr_ub),
      color=point_color,
      linewidth=0.7,
      na.rm=TRUE
    ) +
    geom_point(color=point_color, size=2.4, na.rm=TRUE) +
    scale_x_continuous(breaks=calibration_df$x, labels=calibration_df$bin_label) +
    scale_y_continuous(limits=c(0, 0.5), breaks=seq(0, 0.5, by=0.1)) +
    figure_theme() +
    theme(axis.text.x=element_text(size=9), plot.title=element_text(hjust=0.5)) +
    labs(
      title="Calibration",
      x="Expression-FSR",
      y="Observed FSR"
    )

  return(pp)
}

make_correlation_by_directional_fsr_bin_plot <- function(df, pvalue_threshold=1.0) {
  point_color <- "#007C89"

  df <- df[!is.na(df$directional_FSR) &
             !is.na(df$raw_correlation) &
             !is.na(df$raw_correlation_pvalue) &
             df$raw_correlation_pvalue <= pvalue_threshold, ]

  bin_breaks <- c(0, 0.1, 0.2, 0.3, 0.4, Inf)
  bin_labels <- c("[0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", ">0.4")

  df$directional_FSR_bin <- cut(
    df$directional_FSR,
    breaks=bin_breaks,
    labels=bin_labels,
    include.lowest=TRUE,
    right=TRUE
  )

  mean_correlation_arr <- c()
  mean_correlation_lb_arr <- c()
  mean_correlation_ub_arr <- c()
  n_gene_arr <- c()
  bin_arr <- c()
  x_arr <- c()

  for (bin_iter in 1:length(bin_labels)) {
    bin_name <- bin_labels[bin_iter]
    bin_df <- df[df$directional_FSR_bin == bin_name, ]
    n_genes <- nrow(bin_df)

    if (n_genes == 0) {
      mean_correlation <- NA
      mean_correlation_lb <- NA
      mean_correlation_ub <- NA
    } else {
      corrs <- bin_df$raw_correlation
      mean_correlation <- mean(corrs)

      if (n_genes == 1) {
        mean_correlation_lb <- mean_correlation
        mean_correlation_ub <- mean_correlation
      } else {
        se_correlation <- sd(corrs)/sqrt(n_genes)
        ci_width <- qt(0.975, df=n_genes - 1)*se_correlation
        mean_correlation_lb <- mean_correlation - ci_width
        mean_correlation_ub <- mean_correlation + ci_width
      }
    }

    mean_correlation_arr <- c(mean_correlation_arr, mean_correlation)
    mean_correlation_lb_arr <- c(mean_correlation_lb_arr, mean_correlation_lb)
    mean_correlation_ub_arr <- c(mean_correlation_ub_arr, mean_correlation_ub)
    n_gene_arr <- c(n_gene_arr, n_genes)
    bin_arr <- c(bin_arr, bin_name)
    x_arr <- c(x_arr, bin_iter)
  }

  correlation_df <- data.frame(
    directional_FSR_bin=bin_arr,
    x=x_arr,
    mean_correlation=mean_correlation_arr,
    mean_correlation_lb=mean_correlation_lb_arr,
    mean_correlation_ub=mean_correlation_ub_arr,
    n_genes=n_gene_arr
  )
  correlation_df$directional_FSR_bin <- factor(correlation_df$directional_FSR_bin, levels=bin_labels)
  correlation_df$bin_label <- paste0(correlation_df$directional_FSR_bin, "\nN=", correlation_df$n_genes)

  print(correlation_df)
  y_max <- max(correlation_df$mean_correlation_ub, na.rm=TRUE)*1.08
  if (!is.finite(y_max)) {
    y_max <- 0.5
  }
  y_max <- max(y_max, 0.05)

  pp <- ggplot(correlation_df, aes(x=x, y=mean_correlation)) +
    geom_hline(yintercept=0, linetype="dashed", linewidth=0.5, color="grey60") +
    geom_line(color=point_color, linewidth=0.7, na.rm=TRUE) +
    geom_pointrange(
      aes(ymin=mean_correlation_lb, ymax=mean_correlation_ub),
      color=point_color,
      linewidth=0.7,
      na.rm=TRUE
    ) +
    geom_point(color=point_color, size=2.4, na.rm=TRUE) +
    scale_x_continuous(breaks=correlation_df$x, labels=correlation_df$bin_label) +
    coord_cartesian(ylim=c(0, y_max)) +
    figure_theme() +
    theme(axis.text.x=element_text(size=9), plot.title=element_text(hjust=0.5)) +
    labs(
      title="Accuracy",
      x="Expression-FSR",
      y="Average\nexpression correlation"
    )

  return(pp)
}

make_cumulative_correlation_by_directional_fsr_plot <- function(df, max_fsr=0.5, grid_size=101) {
  df <- df[!is.na(df$directional_FSR) &
             !is.na(df$raw_correlation), ]

  fsr_thresholds <- seq(0, max_fsr, length.out=grid_size)
  mean_correlation_arr <- c()
  mean_correlation_lb_arr <- c()
  mean_correlation_ub_arr <- c()
  n_gene_arr <- c()

  for (fsr_threshold in fsr_thresholds) {
    threshold_df <- df[df$directional_FSR <= fsr_threshold, ]
    n_genes <- nrow(threshold_df)

    if (n_genes == 0) {
      mean_correlation <- NA
      mean_correlation_lb <- NA
      mean_correlation_ub <- NA
    } else {
      corrs <- threshold_df$raw_correlation
      mean_correlation <- mean(corrs)

      if (n_genes == 1) {
        mean_correlation_lb <- mean_correlation
        mean_correlation_ub <- mean_correlation
      } else {
        se_correlation <- sd(corrs)/sqrt(n_genes)
        ci_width <- qt(0.975, df=n_genes - 1)*se_correlation
        mean_correlation_lb <- mean_correlation - ci_width
        mean_correlation_ub <- mean_correlation + ci_width
      }
    }

    mean_correlation_arr <- c(mean_correlation_arr, mean_correlation)
    mean_correlation_lb_arr <- c(mean_correlation_lb_arr, mean_correlation_lb)
    mean_correlation_ub_arr <- c(mean_correlation_ub_arr, mean_correlation_ub)
    n_gene_arr <- c(n_gene_arr, n_genes)
  }

  cumulative_df <- data.frame(
    fsr_threshold=fsr_thresholds,
    mean_correlation=mean_correlation_arr,
    mean_correlation_lb=mean_correlation_lb_arr,
    mean_correlation_ub=mean_correlation_ub_arr,
    n_genes=n_gene_arr
  )

  print(cumulative_df)

  pp <- ggplot(cumulative_df, aes(x=fsr_threshold, y=mean_correlation)) +
    geom_hline(yintercept=0, linetype="dashed", linewidth=0.5, color="grey60") +
    geom_ribbon(
      aes(ymin=mean_correlation_lb, ymax=mean_correlation_ub),
      fill="darkorchid2",
      alpha=0.18,
      color=NA,
      na.rm=TRUE
    ) +
    geom_line(color="darkorchid2", linewidth=0.8, na.rm=TRUE) +
    scale_x_continuous(limits=c(0, max_fsr), breaks=seq(0, max_fsr, by=0.1)) +
    coord_cartesian(ylim=c(0, 0.5)) +
    figure_theme() +
    labs(
      x="Predicted FSR threshold",
      y="Average expression correlation"
    )

  return(pp)
}

make_marginal_effect_ci_coverage_plot <- function(df, confidence_levels=c(seq(0.5, 0.95, by=0.05), 0.99), include_observed_se=TRUE) {
  point_color <- "#007C89"
  reference_color <- "#3F3F46"

  interval_variance <- df$borzoi_ld_predicted_variance_std
  y_axis_label <- "Empirical\ncoverage"
  if (include_observed_se) {
    interval_variance <- interval_variance + df$observed_marginal_effect_se_sq_std
    y_axis_label <- "Empirical\ncoverage"
  }

  df$interval_variance <- interval_variance
  df <- df[is.finite(df$observed_marginal_effect_std) &
             is.finite(df$borzoi_ld_predicted_mean_std) &
             is.finite(df$interval_variance) &
             df$interval_variance > 0, ]

  empirical_coverage_arr <- c()
  empirical_coverage_lb_arr <- c()
  empirical_coverage_ub_arr <- c()
  n_variant_arr <- c()

  for (confidence_level in confidence_levels) {
    z_score <- qnorm((1.0 + confidence_level)/2.0)
    interval_sdev <- sqrt(df$interval_variance)
    interval_lb <- df$borzoi_ld_predicted_mean_std - (z_score*interval_sdev)
    interval_ub <- df$borzoi_ld_predicted_mean_std + (z_score*interval_sdev)
    covered <- df$observed_marginal_effect_std >= interval_lb &
      df$observed_marginal_effect_std <= interval_ub

    n_covered <- sum(covered)
    n_variants <- length(covered)
    empirical_coverage <- n_covered/n_variants
    binom_ci <- binom.test(n_covered, n_variants)$conf.int

    empirical_coverage_arr <- c(empirical_coverage_arr, empirical_coverage)
    empirical_coverage_lb_arr <- c(empirical_coverage_lb_arr, binom_ci[1])
    empirical_coverage_ub_arr <- c(empirical_coverage_ub_arr, binom_ci[2])
    n_variant_arr <- c(n_variant_arr, n_variants)
  }

  coverage_df <- data.frame(
    nominal_coverage=confidence_levels,
    empirical_coverage=empirical_coverage_arr,
    empirical_coverage_lb=empirical_coverage_lb_arr,
    empirical_coverage_ub=empirical_coverage_ub_arr,
    n_variants=n_variant_arr
  )

  print(coverage_df)

  pp <- ggplot(coverage_df, aes(x=nominal_coverage, y=empirical_coverage)) +
    geom_abline(intercept=0, slope=1, color=reference_color, linetype="dashed", linewidth=0.6) +
    geom_line(color=point_color, linewidth=0.7) +
    geom_pointrange(
      aes(ymin=empirical_coverage_lb, ymax=empirical_coverage_ub),
      color=point_color,
      linewidth=0.7
    ) +
    geom_point(color=point_color, size=2.4) +
    scale_x_continuous(
      limits=c(0.5, 1.0),
      breaks=seq(0.5, 1.0, by=0.1),
      labels=function(x) paste0(round(100*x), "%")
    ) +
    scale_y_continuous(
      limits=c(0.5, 1.0),
      breaks=seq(0.5, 1.0, by=0.1),
      labels=function(x) paste0(round(100*x), "%")
    ) +
    figure_theme() +
    labs(
      x="Nominal CI level",
      y=y_axis_label
    )

  return(pp)
}

make_marginal_effect_confident_subset_ci_coverage_plot <- function(df, confidence_levels=c(seq(0.5, 0.95, by=0.05), 0.99), confident_level=0.9, include_observed_se=TRUE) {
  point_color <- "#C44E52"
  reference_color <- "#3F3F46"

  interval_variance <- df$borzoi_ld_predicted_variance_std
  y_axis_label <- "Empirical coverage in confident subset"
  if (include_observed_se) {
    interval_variance <- interval_variance + df$observed_marginal_effect_se_sq_std
  }

  df$interval_variance <- interval_variance
  df <- df[is.finite(df$observed_marginal_effect_std) &
             is.finite(df$borzoi_ld_predicted_mean_std) &
             is.finite(df$interval_variance) &
             df$interval_variance > 0, ]

  confident_z_score <- qnorm((1.0 + confident_level)/2.0)
  interval_sdev <- sqrt(df$interval_variance)
  confident_lb <- df$borzoi_ld_predicted_mean_std - (confident_z_score*interval_sdev)
  confident_ub <- df$borzoi_ld_predicted_mean_std + (confident_z_score*interval_sdev)
  df <- df[confident_lb > 0 | confident_ub < 0, ]
  if (nrow(df) == 0) {
    stop("No variant-gene pairs have a predicted confident CI excluding zero.")
  }

  empirical_coverage_arr <- c()
  empirical_coverage_lb_arr <- c()
  empirical_coverage_ub_arr <- c()
  n_variant_arr <- c()

  for (confidence_level in confidence_levels) {
    z_score <- qnorm((1.0 + confidence_level)/2.0)
    interval_sdev <- sqrt(df$interval_variance)
    interval_lb <- df$borzoi_ld_predicted_mean_std - (z_score*interval_sdev)
    interval_ub <- df$borzoi_ld_predicted_mean_std + (z_score*interval_sdev)
    covered <- df$observed_marginal_effect_std >= interval_lb &
      df$observed_marginal_effect_std <= interval_ub

    n_covered <- sum(covered)
    n_variants <- length(covered)
    empirical_coverage <- n_covered/n_variants
    binom_ci <- binom.test(n_covered, n_variants)$conf.int

    empirical_coverage_arr <- c(empirical_coverage_arr, empirical_coverage)
    empirical_coverage_lb_arr <- c(empirical_coverage_lb_arr, binom_ci[1])
    empirical_coverage_ub_arr <- c(empirical_coverage_ub_arr, binom_ci[2])
    n_variant_arr <- c(n_variant_arr, n_variants)
  }

  coverage_df <- data.frame(
    nominal_coverage=confidence_levels,
    empirical_coverage=empirical_coverage_arr,
    empirical_coverage_lb=empirical_coverage_lb_arr,
    empirical_coverage_ub=empirical_coverage_ub_arr,
    n_variants=n_variant_arr
  )

  print(coverage_df)

  pp <- ggplot(coverage_df, aes(x=nominal_coverage, y=empirical_coverage)) +
    geom_abline(intercept=0, slope=1, color=reference_color, linetype="dashed", linewidth=0.6) +
    geom_line(color=point_color, linewidth=0.7) +
    geom_pointrange(
      aes(ymin=empirical_coverage_lb, ymax=empirical_coverage_ub),
      color=point_color,
      linewidth=0.7
    ) +
    geom_point(color=point_color, size=2.4) +
    scale_x_continuous(
      limits=c(0.5, 1.0),
      breaks=seq(0.5, 1.0, by=0.1),
      labels=function(x) paste0(round(100*x), "%")
    ) +
    scale_y_continuous(
      limits=c(0.5, 1.0),
      breaks=seq(0.5, 1.0, by=0.1),
      labels=function(x) paste0(round(100*x), "%")
    ) +
    figure_theme() +
    labs(
      x="Nominal CI level",
      y=y_axis_label
    )

  return(pp)
}

load_ldscore_grid_squared_prior <- function(borzoi_based_prior_output_stem) {
  summary_file <- paste0(borzoi_based_prior_output_stem, "_ldscore_grid_squared_summary.txt")
  param_file <- paste0(borzoi_based_prior_output_stem, "_ldscore_grid_squared_params.txt")
  curve_file <- paste0(borzoi_based_prior_output_stem, "_ldscore_grid_squared_curve.txt")

  prior_summary <- read.table(summary_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  prior_params <- read.table(param_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  prior_curve <- read.table(curve_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  return(list(
    summary=prior_summary,
    params=prior_params,
    curve=prior_curve
  ))
}

make_prior_variance_curve_plot <- function(prior_curve) {
  prior_curve$pred_allelic_variance <- pmax(prior_curve$pred_allelic_variance, 0.0)
  pp <- ggplot(prior_curve, aes(x=abs_delta_allelic, y=sqrt(pred_allelic_variance))) +
    geom_line(color="#007C89", linewidth=0.8) +
    figure_theme() +
    labs(
      x="Absolute Borzoi allelic effect",
      y="Fitted causal-effect prior SD"
    )

  return(pp)
}

make_prior_density_by_borzoi_magnitude_plot <- function(df, n_bins=4, n_grid=600) {
  df <- df[is.finite(df$borzoi_effect_size_allelic) &
             is.finite(df$borzoi_prior_mean_allelic) &
             is.finite(df$borzoi_prior_variance_allelic) &
             df$borzoi_prior_variance_allelic > 0, ]
  if (nrow(df) == 0) {
    stop("No rows with finite Borzoi prior mean and variance.")
  }

  df$abs_borzoi_effect <- abs(df$borzoi_effect_size_allelic)
  quantile_breaks <- unique(as.numeric(quantile(df$abs_borzoi_effect, probs=seq(0, 1, length.out=n_bins + 1), na.rm=TRUE)))
  if (length(quantile_breaks) < 3) {
    stop("Not enough distinct Borzoi effect values to create magnitude bins.")
  }

  df$borzoi_magnitude_bin <- cut(
    df$abs_borzoi_effect,
    breaks=quantile_breaks,
    include.lowest=TRUE,
    right=TRUE
  )

  x_limit <- quantile(abs(df$borzoi_prior_mean_allelic) + 4.0*sqrt(df$borzoi_prior_variance_allelic), 0.99, na.rm=TRUE)
  x_limit <- max(x_limit, quantile(abs(df$borzoi_prior_mean_allelic), 0.99, na.rm=TRUE))
  if (!is.finite(x_limit) || x_limit == 0) {
    x_limit <- 0.05
  }
  x_grid <- seq(-x_limit, x_limit, length.out=n_grid)

  density_df <- data.frame()
  bin_levels <- levels(df$borzoi_magnitude_bin)
  max_rows_per_bin <- 5000
  set.seed(1)
  for (bin_name in bin_levels) {
    bin_df <- df[df$borzoi_magnitude_bin == bin_name, ]
    if (nrow(bin_df) == 0) {
      next
    }
    if (nrow(bin_df) > max_rows_per_bin) {
      bin_df <- bin_df[sample(1:nrow(bin_df), max_rows_per_bin), ]
    }
    density_vals <- rep(0, length(x_grid))
    for (row_iter in 1:nrow(bin_df)) {
      density_vals <- density_vals + dnorm(
        x_grid,
        mean=bin_df$borzoi_prior_mean_allelic[row_iter],
        sd=sqrt(bin_df$borzoi_prior_variance_allelic[row_iter])
      )
    }
    density_vals <- density_vals/nrow(bin_df)
    bin_label <- paste0(bin_name, "\nN=", nrow(bin_df))
    density_df <- rbind(
      density_df,
      data.frame(
        causal_effect=x_grid,
        density=density_vals,
        borzoi_magnitude_bin=bin_label
      )
    )
  }

  pp <- ggplot(density_df, aes(x=causal_effect, y=density, color=borzoi_magnitude_bin)) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.5, color="grey60") +
    geom_line(linewidth=0.75) +
    figure_theme() +
    labs(
      x="Causal allelic effect",
      y="Fitted prior density",
      color="|Borzoi effect| bin"
    )

  return(pp)
}

make_prior_density_at_selected_borzoi_effects_plot <- function(prior_info, borzoi_effects=c(-0.1, 0.8), n_grid=800) {
  prior_params <- prior_info$params
  prior_curve <- prior_info$curve

  a_prior <- prior_params$value[prior_params$parameter == "a_prior"][1]
  if (is.na(a_prior)) {
    stop("Could not find a_prior in prior parameter file.")
  }

  prior_curve$pred_allelic_variance <- pmax(prior_curve$pred_allelic_variance, 0.0)
  prior_variance <- approx(
    x=prior_curve$abs_delta_allelic,
    y=prior_curve$pred_allelic_variance,
    xout=abs(borzoi_effects),
    rule=2
  )$y
  prior_variance <- pmax(prior_variance, 0.0)
  prior_mean <- a_prior*borzoi_effects
  prior_sd <- sqrt(prior_variance)

  x_limit <- max(abs(prior_mean) + 4.0*prior_sd)
  if (!is.finite(x_limit) || x_limit == 0) {
    x_limit <- 0.05
  }
  x_grid <- seq(-x_limit, x_limit, length.out=n_grid)

  density_df <- data.frame()
  for (effect_iter in 1:length(borzoi_effects)) {
    effect_val <- borzoi_effects[effect_iter]
    effect_label <- paste0("Borzoi effect = ", effect_val)
    density_df <- rbind(
      density_df,
      data.frame(
        causal_effect=x_grid,
        density=dnorm(x_grid, mean=prior_mean[effect_iter], sd=prior_sd[effect_iter]),
        borzoi_effect=effect_label
      )
    )
  }

  density_df$borzoi_effect <- factor(
    density_df$borzoi_effect,
    levels=paste0("Borzoi effect = ", borzoi_effects)
  )

  pp <- ggplot(density_df, aes(x=causal_effect, y=density, color=borzoi_effect)) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.5, color="grey60") +
    geom_line(linewidth=0.8) +
    figure_theme() +
    theme(
      legend.position=c(0.02, 1.51),
      legend.justification=c(0, 1),
      legend.background=element_rect(fill="white", color=NA),
      legend.key=element_blank()
    ) +
    labs(
      x="Causal allelic effect",
      y="\nGDL-Bridge    \nDensity   ",
      color=""
    )

  return(pp)
}

make_marginal_effect_directional_power_plot <- function(df, confidence_level=0.9, include_observed_se=TRUE) {
  negative_color <- "#C44E52"
  positive_color <- "#007C89"

  interval_variance <- df$borzoi_ld_predicted_variance_std
  if (include_observed_se) {
    interval_variance <- interval_variance + df$observed_marginal_effect_se_sq_std
  }

  df$interval_variance <- interval_variance
  df <- df[is.finite(df$observed_marginal_effect_std) &
             is.finite(df$borzoi_ld_predicted_mean_std) &
             is.finite(df$interval_variance) &
             df$interval_variance > 0, ]

  z_score <- qnorm((1.0 + confidence_level)/2.0)
  interval_sdev <- sqrt(df$interval_variance)
  interval_lb <- df$borzoi_ld_predicted_mean_std - (z_score*interval_sdev)
  interval_ub <- df$borzoi_ld_predicted_mean_std + (z_score*interval_sdev)

  df$predicted_direction <- NA
  df$predicted_direction[interval_ub < 0] <- "Predicted negative"
  df$predicted_direction[interval_lb > 0] <- "Predicted positive"
  df$predicted_sign <- NA
  df$predicted_sign[interval_ub < 0] <- -1.0
  df$predicted_sign[interval_lb > 0] <- 1.0
  df <- df[!is.na(df$predicted_direction), ]
  if (nrow(df) == 0) {
    stop("No variant-gene pairs have a predicted CI excluding zero.")
  }

  summary_df <- data.frame(
    predicted_direction=c("Predicted negative", "Predicted positive"),
    n_variant_pairs=c(
      sum(df$predicted_direction == "Predicted negative"),
      sum(df$predicted_direction == "Predicted positive")
    ),
    mean_observed_effect=c(
      mean(df$observed_marginal_effect_std[df$predicted_direction == "Predicted negative"]),
      mean(df$observed_marginal_effect_std[df$predicted_direction == "Predicted positive"])
    ),
    fraction_observed_in_predicted_direction=c(
      mean(df$observed_marginal_effect_std[df$predicted_direction == "Predicted negative"] < 0),
      mean(df$observed_marginal_effect_std[df$predicted_direction == "Predicted positive"] > 0)
    )
  )
  print(summary_df)

  direction_levels <- c("Predicted negative", "Predicted positive")
  direction_labels <- paste0(
    direction_levels,
    "\nN=",
    summary_df$n_variant_pairs[match(direction_levels, summary_df$predicted_direction)]
  )
  names(direction_labels) <- direction_levels
  df$predicted_direction <- factor(df$predicted_direction, levels=direction_levels)

	  pp <- ggplot(df, aes(x=predicted_direction, y=observed_marginal_effect_std, fill=predicted_direction)) +
	    geom_hline(yintercept=0, linetype="dashed", linewidth=0.5, color="grey60") +
	    geom_violin(width=0.72, alpha=0.35, color=NA, trim=FALSE) +
	    geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.85, linewidth=0.45) +
	    geom_jitter(width=0.12, height=0, alpha=0.35, size=0.75, shape=16, color="black") +
	    annotate("text", x=1.5, y=Inf, label="Predicted 90% CI excludes zero", vjust=1.25, size=3.3) +
	    scale_x_discrete(labels=direction_labels) +
	    scale_fill_manual(values=c("Predicted negative"=negative_color, "Predicted positive"=positive_color)) +
	    figure_theme() +
	    theme(
	      legend.position="none",
	      plot.margin=margin(5.5, 5.5, 0, 5.5)
	    ) +
	    labs(
	      x=NULL,
	      y="Observed\neffect"
	    )

  return(pp)
}

load_hybrid_expression_prediction_summary <- function(hybrid_expression_prediction_dir,
                                                      target_tissue="Whole_Blood",
                                                      target_sample="GTEX-1LB8K-0005-SM-DIPED.1",
                                                      anno_method="borzoi_effect_sizes",
                                                      distribution="ldscore_grid_squared",
                                                      standardize_geno="False",
                                                      training_sample_sizes=c(25, 50, 100, 200, 250),
                                                      max_training_sample_size=300) {
  method_columns <- c(
    "borzoi_correlation",
    "no_borzoi_correlation",
    "hybrid_correlation"
  )
  method_labels <- c(
    "borzoi_correlation"="Borzoi",
    "no_borzoi_correlation"="Elastic Net",
    "hybrid_correlation"="Hybrid"
  )

  summary_df <- data.frame()
  for (training_sample_size in training_sample_sizes) {
    hybrid_prediction_file <- paste0(
      hybrid_expression_prediction_dir,
      "expr_pred_summary_",
      target_tissue,
      "_",
      target_sample,
      "_",
      anno_method,
      "_",
      distribution,
      "_",
      standardize_geno,
      "_",
      training_sample_size,
      "_",
      max_training_sample_size,
      ".txt"
    )
    if (!file.exists(hybrid_prediction_file)) {
      warning(paste0("Skipping missing hybrid prediction file: ", hybrid_prediction_file))
      next
    }

    df <- read.table(hybrid_prediction_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    for (method_column in method_columns) {
      corrs <- df[[method_column]]
      corrs <- corrs[is.finite(corrs)]
      n_genes <- length(corrs)
      mean_correlation <- mean(corrs)
      se_correlation <- sd(corrs)/sqrt(n_genes)
      summary_df <- rbind(
        summary_df,
        data.frame(
          training_sample_size=training_sample_size,
          method=method_labels[method_column],
          mean_correlation=mean_correlation,
          mean_correlation_lb=mean_correlation - (1.96*se_correlation),
          mean_correlation_ub=mean_correlation + (1.96*se_correlation),
          n_genes=n_genes
        )
      )
    }
  }

  if (nrow(summary_df) == 0) {
    stop("No hybrid expression prediction files were loaded.")
  }

  summary_df$training_sample_size <- factor(
    summary_df$training_sample_size,
    levels=as.character(training_sample_sizes)
  )
  summary_df$method <- factor(summary_df$method, levels=c("Borzoi", "Elastic Net", "Hybrid"))

  return(summary_df)
}

make_hybrid_expression_prediction_barplot <- function(summary_df) {
  method_colors <- c(
    "Borzoi"="#007C89",
    "Elastic Net"="#7A7A7A",
    "Hybrid"="#C44E52"
  )
  dodge_width <- 0.78

  pp <- ggplot(summary_df, aes(x=training_sample_size, y=mean_correlation, fill=method)) +
    geom_col(position=position_dodge(width=dodge_width), width=0.68, color="white", linewidth=0.25) +
    geom_errorbar(
      aes(ymin=mean_correlation_lb, ymax=mean_correlation_ub),
      position=position_dodge(width=dodge_width),
      width=0.18,
      linewidth=0.45
    ) +
    scale_fill_manual(values=method_colors) +
    figure_theme() +
    labs(
      x="Training sample size",
      y="Mean expression correlation",
      fill=""
    )

  return(pp)
}

make_hybrid_expression_prediction_lineplot <- function(summary_df) {
  summary_df <- summary_df[summary_df$method != "Borzoi", ]
  summary_df$method <- factor(summary_df$method, levels=c("Hybrid", "Elastic Net"))
  method_colors <- c(
    "Hybrid"="#C44E52",
    "Elastic Net"="#7A7A7A"
  )

  pp <- ggplot(summary_df, aes(x=as.numeric(as.character(training_sample_size)), y=mean_correlation, color=method, group=method)) +
    geom_line(linewidth=0.75) +
    geom_pointrange(
      aes(ymin=mean_correlation_lb, ymax=mean_correlation_ub),
      linewidth=0.55,
      size=0.55
    ) +
    scale_color_manual(values=method_colors) +
    scale_x_continuous(breaks=sort(unique(as.numeric(as.character(summary_df$training_sample_size))))) +
    figure_theme() +
    theme(
      legend.position=c(0.47, 0.55),
      legend.justification=c(0, 1),
      legend.background=element_rect(fill="white", color=NA),
      legend.key=element_blank()
    ) +
    labs(
      x="Training sample size",
      y="Mean expression\ncorrelation",
      color=""
    )

  return(pp)
}

load_hybrid_expression_prediction_delta_by_fsr <- function(hybrid_expression_prediction_dir,
                                                           target_tissue="Whole_Blood",
                                                           target_sample="GTEX-1LB8K-0005-SM-DIPED.1",
                                                           anno_method="borzoi_effect_sizes",
                                                           distribution="ldscore_grid_squared",
                                                           standardize_geno="False",
                                                           training_sample_sizes=c(25, 50, 100, 250),
                                                           max_training_sample_size=300,
                                                           fsr_thresholds=c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  delta_df <- data.frame()

  for (training_sample_size in training_sample_sizes) {
    hybrid_prediction_file <- paste0(
      hybrid_expression_prediction_dir,
      "expr_pred_summary_",
      target_tissue,
      "_",
      target_sample,
      "_",
      anno_method,
      "_",
      distribution,
      "_",
      standardize_geno,
      "_",
      training_sample_size,
      "_",
      max_training_sample_size,
      ".txt"
    )
    if (!file.exists(hybrid_prediction_file)) {
      warning(paste0("Skipping missing hybrid prediction file: ", hybrid_prediction_file))
      next
    }

    df <- read.table(hybrid_prediction_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    df$hybrid_no_borzoi_delta <- df$hybrid_correlation - df$no_borzoi_correlation

    for (fsr_threshold in fsr_thresholds) {
      threshold_delta <- df$hybrid_no_borzoi_delta[
        is.finite(df$hybrid_no_borzoi_delta) &
          is.finite(df$borzoi_directional_FSR) &
          df$borzoi_directional_FSR < fsr_threshold
      ]
      n_genes <- length(threshold_delta)
      if (n_genes == 0) {
        mean_delta <- NA
        mean_delta_lb <- NA
        mean_delta_ub <- NA
      } else {
        mean_delta <- mean(threshold_delta)
        se_delta <- sd(threshold_delta)/sqrt(n_genes)
        mean_delta_lb <- mean_delta - (1.96*se_delta)
        mean_delta_ub <- mean_delta + (1.96*se_delta)
      }

      delta_df <- rbind(
        delta_df,
        data.frame(
          training_sample_size=training_sample_size,
          fsr_threshold=fsr_threshold,
          fsr_threshold_label=paste0("Borzoi FSR < ", fsr_threshold),
          mean_delta=mean_delta,
          mean_delta_lb=mean_delta_lb,
          mean_delta_ub=mean_delta_ub,
          n_genes=n_genes
        )
      )
    }
  }

  if (nrow(delta_df) == 0) {
    stop("No hybrid expression prediction files were loaded.")
  }

  delta_df$training_sample_size <- factor(
    delta_df$training_sample_size,
    levels=as.character(training_sample_sizes)
  )
  delta_df$fsr_threshold_label <- factor(
    delta_df$fsr_threshold_label,
    levels=paste0("Borzoi FSR < ", fsr_thresholds)
  )

  return(delta_df)
}

load_hybrid_expression_prediction_delta_by_fsr_bins <- function(hybrid_expression_prediction_dir,
                                                                target_tissue="Whole_Blood",
                                                                target_sample="GTEX-1LB8K-0005-SM-DIPED.1",
                                                                anno_method="borzoi_effect_sizes",
                                                                distribution="ldscore_grid_squared",
                                                                standardize_geno="False",
                                                                training_sample_sizes=c(25, 50, 100, 200, 250),
                                                                max_training_sample_size=300,
                                                                fsr_bin_lower_bounds=c(-Inf, 0.1, 0.3),
                                                                fsr_bin_upper_bounds=c(0.1, 0.3, Inf),
                                                                fsr_bin_labels=c("< .1", "[.1, .3)", ">= .3")) {
  delta_df <- data.frame()

  for (training_sample_size in training_sample_sizes) {
    hybrid_prediction_file <- paste0(
      hybrid_expression_prediction_dir,
      "expr_pred_summary_",
      target_tissue,
      "_",
      target_sample,
      "_",
      anno_method,
      "_",
      distribution,
      "_",
      standardize_geno,
      "_",
      training_sample_size,
      "_",
      max_training_sample_size,
      ".txt"
    )
    if (!file.exists(hybrid_prediction_file)) {
      warning(paste0("Skipping missing hybrid prediction file: ", hybrid_prediction_file))
      next
    }

    df <- read.table(hybrid_prediction_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    df$hybrid_no_borzoi_delta <- df$hybrid_correlation - df$no_borzoi_correlation

    for (bin_iter in seq_along(fsr_bin_labels)) {
      lower_bound <- fsr_bin_lower_bounds[bin_iter]
      upper_bound <- fsr_bin_upper_bounds[bin_iter]
      bin_indices <- is.finite(df$hybrid_no_borzoi_delta) & is.finite(df$borzoi_directional_FSR)
      if (is.finite(lower_bound)) {
        bin_indices <- bin_indices & df$borzoi_directional_FSR >= lower_bound
      }
      if (is.finite(upper_bound)) {
        bin_indices <- bin_indices & df$borzoi_directional_FSR < upper_bound
      }

      bin_delta <- df$hybrid_no_borzoi_delta[bin_indices]
      n_genes <- length(bin_delta)
      if (n_genes == 0) {
        mean_delta <- NA
        mean_delta_lb <- NA
        mean_delta_ub <- NA
      } else {
        mean_delta <- mean(bin_delta)
        se_delta <- sd(bin_delta)/sqrt(n_genes)
        mean_delta_lb <- mean_delta - (1.96*se_delta)
        mean_delta_ub <- mean_delta + (1.96*se_delta)
      }

      delta_df <- rbind(
        delta_df,
        data.frame(
          training_sample_size=training_sample_size,
          fsr_bin_label=fsr_bin_labels[bin_iter],
          mean_delta=mean_delta,
          mean_delta_lb=mean_delta_lb,
          mean_delta_ub=mean_delta_ub,
          n_genes=n_genes
        )
      )
    }
  }

  if (nrow(delta_df) == 0) {
    stop("No hybrid expression prediction files were loaded.")
  }

  delta_df$training_sample_size <- factor(
    delta_df$training_sample_size,
    levels=as.character(training_sample_sizes)
  )
  delta_df$fsr_bin_label <- factor(delta_df$fsr_bin_label, levels=fsr_bin_labels)

  return(delta_df)
}

make_hybrid_expression_delta_by_fsr_plot <- function(delta_df) {
  pp <- ggplot(delta_df, aes(x=training_sample_size, y=mean_delta, color=fsr_threshold_label, group=fsr_threshold_label)) +
    geom_hline(yintercept=0, linetype="dashed", linewidth=0.5, color="grey60") +
    geom_line(linewidth=0.75, na.rm=TRUE) +
    geom_pointrange(
      aes(ymin=mean_delta_lb, ymax=mean_delta_ub),
      linewidth=0.55,
      size=0.55,
      na.rm=TRUE
    ) +
    figure_theme() +
    labs(
      x="Training sample size",
      y="Mean hybrid - Elastic Net correlation",
      color=""
    )

  return(pp)
}

make_hybrid_expression_delta_by_fsr_simple_plot <- function(delta_df, fsr_bin_labels=c("< .1", "[.1, .3)", ">= .3")) {
  plot_df <- delta_df
  plot_df$fsr_bin_label <- factor(plot_df$fsr_bin_label, levels=fsr_bin_labels)
  fsr_colors <- colorRampPalette(c("#F2C94C", "#007C89"))(length(fsr_bin_labels))
  names(fsr_colors) <- fsr_bin_labels

  pp <- ggplot(plot_df, aes(x=as.numeric(as.character(training_sample_size)), y=mean_delta, color=fsr_bin_label, group=fsr_bin_label)) +
    geom_hline(yintercept=0, linetype="dashed", linewidth=0.5, color="grey60") +
    geom_line(linewidth=0.75, na.rm=TRUE) +
    geom_pointrange(
      aes(ymin=mean_delta_lb, ymax=mean_delta_ub),
      linewidth=0.55,
      size=0.55,
      na.rm=TRUE
    ) +
    scale_x_continuous(breaks=sort(unique(as.numeric(as.character(plot_df$training_sample_size))))) +
    scale_color_manual(
      values=fsr_colors,
      breaks=fsr_bin_labels
    ) +
    figure_theme() +
    theme(
      legend.position=c(0.995, 1.078),
      legend.justification=c(1, 1),
      legend.background=element_rect(fill="white", color=NA),
      legend.key=element_blank()
    ) +
    labs(
      x="Training sample size",
      y="Hybrid - Elastic Net\ncorrelation",
      color="Expression FSR"
    )

  return(pp)
}

#######################
# Command line args
#######################
expr_pred_output_dir = args[1]
visualization_dir = args[2]
if (length(args) >= 3 && !is.na(args[3])) {
  hybrid_expression_prediction_dir=args[3]
} else {
  hybrid_expression_prediction_dir <- sub("expression_prediction/?$", "new_hybrid_expression_prediction/", expr_pred_output_dir)
}

borzoi_based_prior_output_dir <- sub("expression_prediction/?$", "borzoi_based_prior/", expr_pred_output_dir)
ref_tissue <- "Whole_Blood"
ref_sample <- "GTEX-1LB8K-0005-SM-DIPED.1"
anno_method <- "borzoi_effect_sizes"
borzoi_based_prior_output_stem <- paste0(
  borzoi_based_prior_output_dir,
  "borzoi_based_prior_ldscore_grid_squared_based_",
  ref_tissue,
  "_",
  ref_sample,
  "_",
  anno_method
)


######################
# Expression prediction
#####################
if (FALSE) {
prediction_summary_file <- paste0(expr_pred_output_dir, "expr_pred_summary_Skin_Not_Sun_Exposed_Suprapubic_GTEX-133LE-2326-SM-5K7W3.1_Whole_Blood_borzoi_effect_sizes_ldscore_grid_squared_False.txt")
prediction_summary_file <- paste0(expr_pred_output_dir, "expr_pred_summary_Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1_Whole_Blood_borzoi_effect_sizes_ldscore_grid_squared_False.txt")
prediction_summary_file <- paste0(expr_pred_output_dir, "expr_pred_summary_Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1_Muscle_Skeletal_borzoi_effect_sizes_ldscore_grid_squared_False.txt")

df <- read.table(prediction_summary_file, header=TRUE, sep="\t")



#directional_fsr_calibration_plot <- make_directional_fsr_calibration_plot(df[df$cis_snp_h2_pvalue < .01,])
directional_fsr_calibration_plot <- make_directional_fsr_calibration_plot(df[df$cis_snp_h2_pvalue < .05,])
directional_fsr_calibration_output_file <- paste0(visualization_dir, "directional_fsr_calibration_observed_negative_raw_correlation.pdf")
ggsave(directional_fsr_calibration_plot, file=directional_fsr_calibration_output_file, device="pdf", width=7.2, height=3.63, units="in")

correlation_by_directional_fsr_bin_plot <- make_correlation_by_directional_fsr_bin_plot(df[df$cis_snp_h2_pvalue < .05,])
correlation_by_directional_fsr_bin_output_file <- paste0(visualization_dir, "average_expression_correlation_by_predicted_fsr_bin.pdf")
ggsave(correlation_by_directional_fsr_bin_plot, file=correlation_by_directional_fsr_bin_output_file, device="pdf", width=7.2, height=3.63, units="in")

prediction_accuracy_fsr_calibration_plot <- plot_grid(
  directional_fsr_calibration_plot,
  correlation_by_directional_fsr_bin_plot,
  labels=c("a", "b"),
  nrow=1,
  align="h",
  axis="tb",
  rel_widths=c(1, 1.1)
)
prediction_accuracy_fsr_calibration_output_file <- paste0(visualization_dir, "prediction_accuracy_and_fsr_calibration_by_predicted_fsr.pdf")
ggsave(prediction_accuracy_fsr_calibration_plot, file=prediction_accuracy_fsr_calibration_output_file, device="pdf", width=7.3, height=2.1, units="in")


######################
# Marginal effect prediction
#####################
marginal_effect_pred_file <- paste0(expr_pred_output_dir, "marginal_eqtl_effect_prediction_Muscle_Skeletal_GTEX-13QJ3-0726-SM-5SI68.1_Whole_Blood_borzoi_effect_sizes_ldscore_grid_squared_False_ld_pruned_r2_0.2_seed_1.txt.gz")
marginal_effect_pred_file <- paste0(expr_pred_output_dir, "marginal_eqtl_effect_prediction_Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1_borzoi_effect_sizes_ldscore_grid_squared_False_ld_pruned_r2_0.2_seed_1.txt.gz")
marginal_effect_pred_file <- paste0(expr_pred_output_dir, "marginal_eqtl_effect_prediction_Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1_Muscle_Skeletal_borzoi_effect_sizes_ldscore_grid_squared_False_ld_pruned_r2_0.2_seed_1.txt.gz")

marginal_effect_df <- read.table(
  gzfile(marginal_effect_pred_file),
  header=TRUE,
  sep="\t",
  stringsAsFactors=FALSE
)

ldscore_grid_squared_prior <- load_ldscore_grid_squared_prior(borzoi_based_prior_output_stem)

prior_variance_curve_plot <- make_prior_variance_curve_plot(ldscore_grid_squared_prior$curve)
prior_variance_curve_output_file <- paste0(visualization_dir, "borzoi_prior_variance_curve.pdf")
ggsave(prior_variance_curve_plot, file=prior_variance_curve_output_file, device="pdf", width=4.2, height=3.6, units="in")
print(prior_variance_curve_output_file)

prior_density_by_borzoi_magnitude_plot <- make_prior_density_by_borzoi_magnitude_plot(marginal_effect_df)
prior_density_by_borzoi_magnitude_output_file <- paste0(visualization_dir, "borzoi_prior_density_by_borzoi_magnitude.pdf")
ggsave(prior_density_by_borzoi_magnitude_plot, file=prior_density_by_borzoi_magnitude_output_file, device="pdf", width=5.2, height=3.6, units="in")
print(prior_density_by_borzoi_magnitude_output_file)

prior_density_selected_borzoi_effects_plot <- make_prior_density_at_selected_borzoi_effects_plot(ldscore_grid_squared_prior)
prior_density_selected_borzoi_effects_output_file <- paste0(visualization_dir, "borzoi_prior_density_selected_borzoi_effects.pdf")
ggsave(prior_density_selected_borzoi_effects_plot, file=prior_density_selected_borzoi_effects_output_file, device="pdf", width=7.4, height=2.0, units="in")
print(prior_density_selected_borzoi_effects_output_file)

prior_distribution_shape_plot <- plot_grid(
  prior_variance_curve_plot,
  prior_density_by_borzoi_magnitude_plot,
  labels=c("A", "B"),
  nrow=1,
  align="h",
  axis="tb",
  rel_widths=c(1.0, 1.2)
)
prior_distribution_shape_output_file <- paste0(visualization_dir, "borzoi_prior_distribution_shape.pdf")
ggsave(prior_distribution_shape_plot, file=prior_distribution_shape_output_file, device="pdf", width=8.0, height=3.6, units="in")
print(prior_distribution_shape_output_file)

marginal_effect_ci_coverage_plot <- make_marginal_effect_ci_coverage_plot(marginal_effect_df)
marginal_effect_ci_coverage_output_file <- paste0(visualization_dir, "marginal_effect_ci_empirical_coverage.pdf")
ggsave(marginal_effect_ci_coverage_plot, file=marginal_effect_ci_coverage_output_file, device="pdf", width=4.2, height=3.6, units="in")
print(marginal_effect_ci_coverage_output_file)

marginal_effect_confident_subset_ci_coverage_plot <- make_marginal_effect_confident_subset_ci_coverage_plot(marginal_effect_df)
marginal_effect_confident_subset_ci_coverage_output_file <- paste0(visualization_dir, "marginal_effect_confident_subset_ci_empirical_coverage.pdf")
ggsave(marginal_effect_confident_subset_ci_coverage_plot, file=marginal_effect_confident_subset_ci_coverage_output_file, device="pdf", width=4.2, height=3.6, units="in")
print(marginal_effect_confident_subset_ci_coverage_output_file)

marginal_effect_directional_power_plot <- make_marginal_effect_directional_power_plot(marginal_effect_df)
marginal_effect_directional_power_output_file <- paste0(visualization_dir, "marginal_effect_predicted_direction_observed_effects.pdf")
ggsave(marginal_effect_directional_power_plot, file=marginal_effect_directional_power_output_file, device="pdf", width=4.2, height=3.6, units="in")
print(marginal_effect_directional_power_output_file)

marginal_effect_coverage_power_plot <- plot_grid(
  marginal_effect_ci_coverage_plot,
  marginal_effect_directional_power_plot,
  labels=c("b", "c"),
  label_y=1.07,
  nrow=1,
  align="h",
  axis="tb"
)
joint <- plot_grid(prior_density_selected_borzoi_effects_plot, marginal_effect_coverage_power_plot, labels=c("a",""), ncol=1, rel_heights=c(1, 1.2))

marginal_effect_coverage_power_output_file <- paste0(visualization_dir, "marginal_effect_coverage_and_directional_power.pdf")
ggsave(joint, file=marginal_effect_coverage_power_output_file, device="pdf", width=7.4, height=2.5, units="in")
print(marginal_effect_coverage_power_output_file)

}

######################
# Hybrid expression prediction
#####################
hybrid_expression_prediction_summary_df <- load_hybrid_expression_prediction_summary(
  hybrid_expression_prediction_dir=hybrid_expression_prediction_dir
)
print(hybrid_expression_prediction_summary_df)

hybrid_expression_prediction_barplot <- make_hybrid_expression_prediction_barplot(hybrid_expression_prediction_summary_df)
hybrid_expression_prediction_barplot_output_file <- paste0(visualization_dir, "hybrid_expression_prediction_mean_correlation_by_sample_size.pdf")
ggsave(hybrid_expression_prediction_barplot, file=hybrid_expression_prediction_barplot_output_file, device="pdf", width=7.4, height=3.6, units="in")
print(hybrid_expression_prediction_barplot_output_file)

hybrid_expression_prediction_delta_by_fsr_df <- load_hybrid_expression_prediction_delta_by_fsr(
  hybrid_expression_prediction_dir=hybrid_expression_prediction_dir
)
print(hybrid_expression_prediction_delta_by_fsr_df)

hybrid_expression_prediction_delta_by_fsr_plot <- make_hybrid_expression_delta_by_fsr_plot(hybrid_expression_prediction_delta_by_fsr_df)
hybrid_expression_prediction_delta_by_fsr_output_file <- paste0(visualization_dir, "hybrid_expression_prediction_delta_by_borzoi_fsr.pdf")
ggsave(hybrid_expression_prediction_delta_by_fsr_plot, file=hybrid_expression_prediction_delta_by_fsr_output_file, device="pdf", width=7.4, height=3.6, units="in")
print(hybrid_expression_prediction_delta_by_fsr_output_file)

hybrid_expression_prediction_delta_by_fsr_bin_df <- load_hybrid_expression_prediction_delta_by_fsr_bins(
  hybrid_expression_prediction_dir=hybrid_expression_prediction_dir
)
print(hybrid_expression_prediction_delta_by_fsr_bin_df)

hybrid_expression_prediction_lineplot <- make_hybrid_expression_prediction_lineplot(hybrid_expression_prediction_summary_df)
hybrid_expression_prediction_delta_by_fsr_simple_plot <- make_hybrid_expression_delta_by_fsr_simple_plot(hybrid_expression_prediction_delta_by_fsr_bin_df)
hybrid_expression_prediction_joint_plot <- plot_grid(
  hybrid_expression_prediction_lineplot,
  hybrid_expression_prediction_delta_by_fsr_simple_plot,
  labels=c("a", "b"),
  nrow=1,
  align="h",
  axis="tb"
)
hybrid_expression_prediction_joint_output_file <- paste0(visualization_dir, "hybrid_expression_prediction_summary_and_delta.pdf")
ggsave(hybrid_expression_prediction_joint_plot, file=hybrid_expression_prediction_joint_output_file, device="pdf", width=7.4, height=2.43, units="in")
print(hybrid_expression_prediction_joint_output_file)
