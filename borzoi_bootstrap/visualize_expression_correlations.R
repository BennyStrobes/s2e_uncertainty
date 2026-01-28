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

  mean_arr <- c()
  mean_lb_arr <- c()
  mean_ub_arr <- c()
  n_gene_arr <- c()
  bin_names_arr <- c()

  indices <- abs(df$correlation_bs_pvalue) <= 1.0
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 1.0", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)

  indices <- abs(df$correlation_bs_pvalue) <= 0.5
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 0.5", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)


  indices <- abs(df$correlation_bs_pvalue) <= 0.4
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 0.4", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)

  indices <- abs(df$correlation_bs_pvalue) <= 0.3
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 0.3", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)

  indices <- abs(df$correlation_bs_pvalue) <= 0.2
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 0.2", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)

  indices <- abs(df$correlation_bs_pvalue) <= 0.15
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 0.15", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)

  indices <- abs(df$correlation_bs_pvalue) <= 0.1
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 0.1", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)


  indices <- abs(df$correlation_bs_pvalue) <= 0.05
  corrs = df$expression_correlation[indices]
  mean_corr = mean(corrs)
  se_corr = sd(corrs)/sqrt(length(corrs))
  mean_arr <- c(mean_arr, mean_corr)
  mean_lb_arr <- c(mean_lb_arr, mean_corr - (1.96*se_corr))
  mean_ub_arr <- c(mean_ub_arr, mean_corr + (1.96*se_corr))
  n_gene_arr <- c(n_gene_arr, sum(indices))
  namer <- paste0("p <= 0.05", "\n", "N=", sum(indices))
  bin_names_arr <- c(bin_names_arr, namer)


  df2 <- data.frame(
    mean_correlation    = mean_arr,
    mean_correlation_lb = mean_lb_arr,
    mean_correlation_ub = mean_ub_arr,
    names    = bin_names_arr,
    n_genes = n_gene_arr
  )


  df2$names = factor(df2$names, levels=as.character(df2$names))
  pp <- ggplot(df2, aes(x = names, y = mean_correlation)) +
    # baseline at chance agreement
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey60") +

    # cleaner CI + point
    geom_pointrange(aes(ymin = mean_correlation_lb, ymax = mean_correlation_ub),
                    color = "darkorchid2", linewidth = 0.7) +
    geom_point(color = "darkorchid2", size = 2.4) +

    figure_theme() +
    labs(
      x = "Borzoi significance bins",
      y = "Average expression correlation"
    )
  print(pp)

}

make_max_z_correlation_scatterplot <- function(df, titler, color) {

	df$neglog10_p <- -log10(2 * pnorm(-df$max_abs_eqtl_z))

	pp <- ggplot(df, aes(x = max_abs_eqtl_z, y = expression_correlation)) +
  		geom_point(alpha = 0.7, size = 1,color=color) +
  		labs(x="Max(abs(eQTL Z))", y="Expression correlation", title=titler) +
  		figure_theme()+ 
  		geom_hline(yintercept = 0, color = "black") +
  		coord_cartesian(ylim = c(-.55, .55)) +
  		 theme(plot.title = element_text(hjust = 0.5))
  	return(pp)
}

make_neg_log10_p_correlation_scatterplot <- function(df, titler, color) {

	df$neglog10_p <- -log10(2 * pnorm(-df$max_abs_eqtl_z))

	pp <- ggplot(df, aes(x = neglog10_p, y = expression_correlation)) +
  		geom_point(alpha = 0.7, size = 1, color=color) +
  		labs(x="-log10(p)", y="Expression correlation", title=titler) +
  		figure_theme()+ 
  		geom_hline(yintercept = 0, color = "black") +
  		coord_cartesian(ylim = c(-.55, .55)) +
  		 theme(plot.title = element_text(hjust = 0.5))
  	return(pp)
}


make_bs_distribution <- function(borzoi_bs, gene_id, color, xlim = NULL) {
  df0 <- data.frame(borzoi_bs = borzoi_bs)

  print(gene_id)

  mu <- mean(borzoi_bs)
  ci <- quantile(borzoi_bs, probs = c(0.025, 0.975))

  print(ci)
  print(mu)


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
      title = paste0(sig_label, "\n", gene_id),
      #subtitle = paste0("Gene: ", gene_id, "\n Variant: ", variant_id),
      x = "Cross-individual expression correlation",
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

#####################
# Command line args
#####################
tissue_name = args[1]
expr_correlation_summary_file = args[2]
visualize_expression_correlation_dir = args[3]


# Load in data
df = read.table(expr_correlation_summary_file, header=TRUE, sep="\t")



geneid ="ENSG00000235098.9"
x0 <- read.table(paste0("/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_bootstrap/expression_correlation/Muscle_Skeletal_", geneid, ".txt"), header=FALSE)$V1
bs_distribution0 <- make_bs_distribution(
  x0, geneid, "#E07A7A"
)


geneid="ENSG00000169598.17"
x1 <- read.table(paste0("/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_bootstrap/expression_correlation/Muscle_Skeletal_", geneid, ".txt"), header=FALSE)$V1
bs_distribution1 <- make_bs_distribution(
  x1, geneid,"#7FBF7B"
)

joint_bs_distribution <- plot_grid(bs_distribution1, bs_distribution0, ncol = 2, labels=c("a","b"))

output_file <- paste0(visualize_expression_correlation_dir, "example_bootstrapped_sampling_distribution_plots.pdf")
ggsave(joint_bs_distribution, file=output_file,device = cairo_pdf, width=7.2, height=3.13, units="in")



#####################################
# Histogram of expresssion correlations
######################################
if (FALSE) {
# Plot histogram of expression correlations with all genes
all_genes_expression_correlation_histogram = make_expression_correlation_histogram(df, "grey50", "All genes")
# Plot histogram of expression correlations with sig genes
sig_genes_expression_correlation_histogram = make_expression_correlation_histogram(df[df$correlation_bs_pvalue <= 0.1,], "darkorchid2", "Significant genes (p <= 0.1)")

# Join to gether with cowplot
joint_hist <- plot_grid(all_genes_expression_correlation_histogram, sig_genes_expression_correlation_histogram, ncol=2)
# save to output
output_file <- paste0(visualize_expression_correlation_dir, "expression_correlation_histogram.pdf")
ggsave(joint_hist, file=output_file,device = cairo_pdf, width=5.2, height=3.63, units="in")
}



#####################################
# Average correlation mean-standard error plot
#####################################
if (FALSE) {
mean_correlation_standard_error_dot_plot = make_mean_correlation_standard_error_plot(df)
output_file <- paste0(visualize_expression_correlation_dir, "average_expression_correlation_mean_se.pdf")
ggsave(mean_correlation_standard_error_dot_plot, file=output_file,device = cairo_pdf, width=7.2, height=3.63, units="in")
}


if (FALSE) {
#####################################
# Correlation vs max Z scatter
#####################################
corr_max_z_scatter_all_genes = make_max_z_correlation_scatterplot(df, "All genes", "grey50")
corr_max_z_scatter_sig_genes = make_max_z_correlation_scatterplot(df[df$correlation_bs_pvalue <= 0.1, ], "Significant genes (p <= 0.1)", "darkorchid2")

joint_scatter <- plot_grid(corr_max_z_scatter_all_genes, corr_max_z_scatter_sig_genes, ncol=2)

output_file <- paste0(visualize_expression_correlation_dir, "max_z_vs_correlation_scatter.pdf")
ggsave(joint_scatter, file=output_file,device = cairo_pdf, width=7.2, height=3.43, units="in")


#####################################
# Correlation vs -log10(p) scatter
#####################################
corr_neg_log10p_scatter_all_genes = make_neg_log10_p_correlation_scatterplot(df, "All genes", "grey50")
corr_neg_log10p_scatter_sig_genes = make_neg_log10_p_correlation_scatterplot(df[df$correlation_bs_pvalue <= 0.1, ], "Significant genes (p <= 0.1)", "darkorchid2")

joint_scatter <- plot_grid(corr_neg_log10p_scatter_all_genes, corr_neg_log10p_scatter_sig_genes, ncol=2)

output_file <- paste0(visualize_expression_correlation_dir, "neg_log10_p_vs_correlation_scatter.pdf")
ggsave(joint_scatter, file=output_file,device = cairo_pdf, width=7.2, height=3.43, units="in")
}





