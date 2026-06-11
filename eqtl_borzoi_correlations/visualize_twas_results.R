args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title=element_text(face="plain", size=11),
		text=element_text(size=11),
		axis.text=element_text(size=11),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		panel.background=element_blank(),
		axis.line=element_line(colour="black"),
		legend.text=element_text(size=11),
		legend.title=element_text(size=11)))
}

normalize_dir <- function(input_dir) {
	if (grepl("/$", input_dir)) {
		return(input_dir)
	}
	return(paste0(input_dir, "/"))
}

sample_proportion_ci <- function(n_successes, n_trials, confidence=0.95) {
	if (n_trials == 0) {
		return(c(NA, NA))
	}
	z <- qnorm(1.0 - (1.0 - confidence)/2.0)
	p_hat <- n_successes/n_trials
	se <- sqrt(p_hat*(1.0 - p_hat)/n_trials)
	return(c(max(0.0, p_hat - z*se), min(1.0, p_hat + z*se)))
}

make_discordant_sign_rate_by_expression_fsr_plot <- function(df) {
	df <- df[is.finite(df$twas_zed) & is.finite(df$borzoi_zed) & is.finite(df$borzoi_fsr), ]
	df <- df[df$twas_zed != 0.0 & df$borzoi_zed != 0.0, ]
	df <- df[abs(df$twas_zed) > 2.0 & abs(df$borzoi_zed) > 2.0, ]
	df$discordant_sign <- sign(df$twas_zed) != sign(df$borzoi_zed)
	df$expression_fsr_bin <- cut(df$borzoi_fsr,
		breaks=c(0.0, 0.1, 0.2, 0.3, 0.4, Inf),
		labels=c("[0,.1]", "(.1,.2]", "(.2,.3]", "(.3,.4]", ">.4"),
		include.lowest=TRUE,
		right=TRUE)

	bin_labels <- levels(df$expression_fsr_bin)
	plot_df <- data.frame()
	for (bin_label in bin_labels) {
		bin_df <- df[df$expression_fsr_bin == bin_label, ]
		n_genes <- nrow(bin_df)
			n_discordant <- sum(bin_df$discordant_sign)
			discordant_rate <- n_discordant/n_genes
			ci <- sample_proportion_ci(n_discordant, n_genes)
		plot_df <- rbind(plot_df, data.frame(
			expression_fsr_bin=bin_label,
			n_genes=n_genes,
			n_discordant=n_discordant,
			discordant_rate=discordant_rate,
			discordant_rate_lb=ci[1],
			discordant_rate_ub=ci[2]))
	}
	plot_df$expression_fsr_bin <- factor(plot_df$expression_fsr_bin, levels=bin_labels)

	print(plot_df)
	pp <- ggplot(plot_df, aes(x=expression_fsr_bin, y=discordant_rate)) +
		geom_errorbar(aes(ymin=discordant_rate_lb, ymax=discordant_rate_ub), width=0.18, color="#3F3F46", linewidth=0.5) +
		geom_point(color="#007C89", size=2.2) +
		labs(x="Expression FSR bin", y="SEWAS-TWAS\nSign Discordance Rate") +
		figure_theme()
	return(pp)
}

if (length(args) != 1) {
	stop("Usage: Rscript visualize_twas_results.R <twas_dir>")
}

twas_dir <- normalize_dir(args[1])
matched_tissue_fsr_summary_file <- paste0(twas_dir, "matched_tissue_fsr_summary.txt")
scatterplot_output_file <- paste0(twas_dir, "matched_tissue_twas_vs_borzoi_z_scatter.pdf")
discordant_sign_rate_output_file <- paste0(twas_dir, "matched_tissue_discordant_sign_rate_by_expression_fsr.pdf")

df <- read.table(matched_tissue_fsr_summary_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
df$borzoi_fsr_lt_0.1 <- df$borzoi_fsr < 0.1

scatterplot <- ggplot(df, aes(x=twas_zed, y=borzoi_zed)) +
	geom_hline(yintercept=0.0, linetype="dashed", color="gray60", linewidth=0.3) +
	geom_vline(xintercept=0.0, linetype="dashed", color="gray60", linewidth=0.3) +
	geom_point(data=df[!df$borzoi_fsr_lt_0.1, ], color="gray60", alpha=0.55, size=1.2) +
	geom_point(data=df[df$borzoi_fsr_lt_0.1, ], aes(color=borzoi_fsr_lt_0.1), alpha=0.8, size=1.2) +
	scale_color_manual(values=c("FALSE"="gray60", "TRUE"="firebrick3"), name="Borzoi FSR < 0.1") +
	labs(x="TWAS z-score", y="Borzoi TWAS z-score") +
	figure_theme()

ggsave(scatterplot_output_file, scatterplot, device="pdf", width=4.8, height=4.2, units="in")
print(scatterplot_output_file)

discordant_sign_rate_plot <- make_discordant_sign_rate_by_expression_fsr_plot(df)
ggsave(discordant_sign_rate_output_file, discordant_sign_rate_plot, device="pdf", width=4.9, height=2.0, units="in")
print(discordant_sign_rate_output_file)

