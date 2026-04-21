args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(
		theme_minimal(base_size=11) +
		theme(
			plot.title=element_text(face="bold", size=12, hjust=.5),
			axis.title=element_text(size=11, color="#25313B"),
			axis.text=element_text(size=10, color="#33424F"),
			panel.grid.minor=element_blank(),
			panel.grid.major.x=element_blank(),
			panel.grid.major.y=element_blank(),
			panel.background=element_rect(fill="#FBFCFD", color=NA),
			plot.background=element_rect(fill="#FBFCFD", color=NA),
			axis.line.x=element_line(colour="#7E8A97", linewidth=.35),
			axis.line.y=element_line(colour="#7E8A97", linewidth=.35),
			legend.position="none"
		)
	)
}

make_bias_barplot <- function(summary_df, plot_title) {
	summary_df$annotation_name <- factor(summary_df$annotation_name, levels=summary_df$annotation_name)
	summary_df$xpos <- seq_len(nrow(summary_df))
	return(
		ggplot(summary_df, aes(x=annotation_name, y=mean)) +
			geom_col(fill="#688B84", color="#365854", width=.68, linewidth=.35) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.16, linewidth=.45, color="#22313B") +
			geom_segment(aes(x=xpos-.30, xend=xpos+.30, y=true_simulated_value, yend=true_simulated_value), color="#D06A4B", linetype="dashed", linewidth=.75, inherit.aes=FALSE) +
			geom_point(aes(y=true_simulated_value), color="#D06A4B", size=1.8, inherit.aes=TRUE) +
			xlab("") +
			ylab("Value") +
			ggtitle(plot_title) +
			figure_theme() +
			theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1))
	)
}

make_coverage_plot <- function(coverage_df, plot_title) {
	coverage_df$annotation_name <- factor(coverage_df$annotation_name, levels=unique(coverage_df$annotation_name))
	coverage_df$coverage_label <- factor(coverage_df$coverage_label, levels=c("99%", "95%", "90%"))
	y_min <- min(c(coverage_df$ci_lower, coverage_df$nominal_coverage), na.rm=TRUE)
	y_max <- max(c(coverage_df$ci_upper, coverage_df$nominal_coverage), na.rm=TRUE)
	y_padding <- 0.02
	y_min <- max(0, y_min - y_padding)
	y_max <- min(1, y_max + y_padding)
	dodge_width <- 0.35
	return(
		ggplot(coverage_df, aes(x=annotation_name, y=observed_coverage, color=coverage_label, group=coverage_label)) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.14, linewidth=.45, position=position_dodge(width=dodge_width)) +
			geom_line(aes(group=coverage_label), linewidth=.7, position=position_dodge(width=dodge_width)) +
			geom_point(size=2, position=position_dodge(width=dodge_width)) +
			geom_hline(aes(yintercept=nominal_coverage, color=coverage_label), linetype="dotted", linewidth=.6, show.legend=FALSE) +
			scale_color_manual(values=c("90%"="#4C78A8", "95%"="#D06A4B", "99%"="#7A5EA6"), name="") +
			xlab("") +
			ylab("Coverage") +
			ggtitle(plot_title) +
			coord_cartesian(ylim=c(y_min, y_max)) +
			figure_theme() +
			theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1), legend.position="right")
	)
}

extract_simulation_iter_from_file <- function(file_name) {
	matches <- regmatches(file_name, regexpr("sim[0-9]+", file_name))
	if (length(matches) == 0 || matches[1] == "") {
		return(NA)
	}
	return(matches[1])
}

load_estimated_corr_summary_across_sims <- function(correlation_inference_results_dir) {
	result_files <- list.files(correlation_inference_results_dir, pattern="_bootstrap_stats\\.txt$", full.names=TRUE)
	if (length(result_files) == 0) {
		stop(paste("No estimated result files found in", correlation_inference_results_dir))
	}

	all_est <- data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp <- read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp$simulation_name <- basename(result_file)
		all_est <- rbind(all_est, tmp)
	}
	if (nrow(all_est) == 0) {
		stop(paste("Estimated result files were found in", correlation_inference_results_dir, "but all were empty"))
	}

	output_names <- c("correlation", "calibration_slope")
	summary_list <- list()
	for (output_name in output_names) {
		subset_df <- all_est[all_est$output_name == output_name, ]
		annotation_names <- unique(subset_df$annotation_name)
		annotation_names <- annotation_names[order(annotation_names)]

		output_summary <- data.frame(
			annotation_name=annotation_names,
			mean=rep(NA, length(annotation_names)),
			ci_lower=rep(NA, length(annotation_names)),
			ci_upper=rep(NA, length(annotation_names)),
			n_sims=rep(NA, length(annotation_names)),
			stringsAsFactors=FALSE
		)

		for (anno_iter in seq_along(annotation_names)) {
			anno_name <- annotation_names[anno_iter]
			anno_df <- subset_df[subset_df$annotation_name == anno_name, ]
			anno_means <- anno_df$mean
			n_sims <- length(anno_means)
			mean_val <- mean(anno_means)
			sem_val <- sd(anno_means)/sqrt(n_sims)

			output_summary$mean[anno_iter] <- mean_val
			output_summary$ci_lower[anno_iter] <- mean_val - 1.96*sem_val
			output_summary$ci_upper[anno_iter] <- mean_val + 1.96*sem_val
			output_summary$n_sims[anno_iter] <- n_sims
		}
		summary_list[[output_name]] <- output_summary
	}
	return(summary_list)
}

load_estimated_corr_per_sim <- function(correlation_inference_results_dir) {
	result_files <- list.files(correlation_inference_results_dir, pattern="_bootstrap_stats\\.txt$", full.names=TRUE)
	if (length(result_files) == 0) {
		stop(paste("No estimated result files found in", correlation_inference_results_dir))
	}

	all_est <- data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp <- read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp$simulation_iter <- extract_simulation_iter_from_file(basename(result_file))
		all_est <- rbind(all_est, tmp)
	}
	if (nrow(all_est) == 0) {
		stop(paste("Estimated result files were found in", correlation_inference_results_dir, "but all were empty"))
	}
	return(all_est)
}

load_true_corr_summary_across_sims <- function(correlation_borzoi_est_effect_dir) {
	result_files <- list.files(correlation_borzoi_est_effect_dir, pattern="_true_sim_effect_summary\\.txt$", full.names=TRUE)
	if (length(result_files) == 0) {
		stop(paste("No true simulated summary files found in", correlation_borzoi_est_effect_dir))
	}

	all_true <- data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp <- read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp$simulation_name <- basename(result_file)
		all_true <- rbind(all_true, tmp)
	}
	if (nrow(all_true) == 0) {
		stop(paste("True simulated summary files were found in", correlation_borzoi_est_effect_dir, "but all were empty"))
	}

	annotation_names <- unique(all_true$annotation_name)
	annotation_names <- annotation_names[order(annotation_names)]

	true_correlation_df <- data.frame(
		annotation_name=annotation_names,
		true_simulated_value=rep(NA, length(annotation_names)),
		stringsAsFactors=FALSE
	)
	true_calibration_slope_df <- data.frame(
		annotation_name=annotation_names,
		true_simulated_value=rep(NA, length(annotation_names)),
		stringsAsFactors=FALSE
	)

	for (anno_iter in seq_along(annotation_names)) {
		anno_name <- annotation_names[anno_iter]
		anno_df <- all_true[all_true$annotation_name == anno_name, ]
		true_correlation_df$true_simulated_value[anno_iter] <- mean(anno_df$pearson_correlation)
		true_calibration_slope_df$true_simulated_value[anno_iter] <- mean(anno_df$regression_slope)
	}

	return(list(correlation=true_correlation_df, calibration_slope=true_calibration_slope_df))
}

load_true_corr_per_sim <- function(correlation_borzoi_est_effect_dir) {
	result_files <- list.files(correlation_borzoi_est_effect_dir, pattern="_true_sim_effect_summary\\.txt$", full.names=TRUE)
	if (length(result_files) == 0) {
		stop(paste("No true simulated summary files found in", correlation_borzoi_est_effect_dir))
	}

	all_true <- data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp <- read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp$simulation_iter <- extract_simulation_iter_from_file(basename(result_file))
		all_true <- rbind(all_true, tmp)
	}
	if (nrow(all_true) == 0) {
		stop(paste("True simulated summary files were found in", correlation_borzoi_est_effect_dir, "but all were empty"))
	}
	return(all_true)
}

make_coverage_summary <- function(all_est, all_true, output_name, true_column_name) {
	est_df <- all_est[all_est$output_name == output_name, c("simulation_iter", "annotation_name", "mean", "bootstrap_se")]
	true_df <- all_true[, c("simulation_iter", "annotation_name", true_column_name)]
	colnames(true_df)[3] <- "true_value"
	merged_df <- merge(est_df, true_df, by=c("simulation_iter", "annotation_name"))
	coverage_levels <- c(.90, .95, .99)
	coverage_summary <- data.frame()
		for (coverage_level in coverage_levels) {
			z_value <- qnorm((1.0 + coverage_level)/2.0)
			tmp_df <- merged_df
			tmp_df$ci_lower <- tmp_df$mean - z_value*tmp_df$bootstrap_se
			tmp_df$ci_upper <- tmp_df$mean + z_value*tmp_df$bootstrap_se
			tmp_df$covered <- (tmp_df$true_value >= tmp_df$ci_lower) & (tmp_df$true_value <= tmp_df$ci_upper)
		annotation_names <- unique(tmp_df$annotation_name)
		annotation_names <- annotation_names[order(annotation_names)]
			for (anno_name in annotation_names) {
				anno_df <- tmp_df[tmp_df$annotation_name == anno_name, ]
				observed_coverage <- mean(anno_df$covered)
				n_sims <- nrow(anno_df)
				coverage_sem <- sqrt(observed_coverage*(1.0-observed_coverage)/n_sims)
				coverage_ci_lower <- max(0, observed_coverage - 1.96*coverage_sem)
				coverage_ci_upper <- min(1, observed_coverage + 1.96*coverage_sem)
				coverage_summary <- rbind(
					coverage_summary,
					data.frame(
						annotation_name=anno_name,
						coverage_label=paste0(as.integer(coverage_level*100), "%"),
						nominal_coverage=coverage_level,
						observed_coverage=observed_coverage,
						coverage_sem=coverage_sem,
						ci_lower=coverage_ci_lower,
						ci_upper=coverage_ci_upper,
						stringsAsFactors=FALSE
					)
				)
			}
	}
	return(coverage_summary)
}


####################
# Command line args
####################
correlation_inference_results_dir <- args[1]
correlation_borzoi_est_effect_dir <- args[2]
viz_dir <- args[3]


#################
# Load in data
################
estimated_result_summary <- load_estimated_corr_summary_across_sims(correlation_inference_results_dir)
true_result_summary <- load_true_corr_summary_across_sims(correlation_borzoi_est_effect_dir)
all_estimated_results <- load_estimated_corr_per_sim(correlation_inference_results_dir)
all_true_results <- load_true_corr_per_sim(correlation_borzoi_est_effect_dir)
estimated_correlation_summary_df <- estimated_result_summary[["correlation"]]
estimated_calibration_slope_summary_df <- estimated_result_summary[["calibration_slope"]]
estimated_correlation_summary_df <- merge(estimated_correlation_summary_df, true_result_summary[["correlation"]], by="annotation_name", all.x=TRUE, sort=FALSE)
estimated_calibration_slope_summary_df <- merge(estimated_calibration_slope_summary_df, true_result_summary[["calibration_slope"]], by="annotation_name", all.x=TRUE, sort=FALSE)


#################
# Bias plots
#################
write.table(estimated_correlation_summary_df, file=file.path(viz_dir, "estimated_correlation_summary_across_sims.txt"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(estimated_calibration_slope_summary_df, file=file.path(viz_dir, "estimated_calibration_slope_summary_across_sims.txt"), sep="\t", quote=FALSE, row.names=FALSE)

correlation_bias_plot <- make_bias_barplot(estimated_correlation_summary_df, "Correlation")
calibration_bias_plot <- make_bias_barplot(estimated_calibration_slope_summary_df, "Calibration Slope")

ggsave(filename=file.path(viz_dir, "estimated_correlation_bias_barplot.pdf"), plot=correlation_bias_plot, width=7.2, height=3.5)
ggsave(filename=file.path(viz_dir, "estimated_calibration_slope_bias_barplot.pdf"), plot=calibration_bias_plot, width=7.2, height=3.5)


#################
# Coverage plots
#################
correlation_coverage_df <- make_coverage_summary(all_estimated_results, all_true_results, "correlation", "pearson_correlation")
calibration_coverage_df <- make_coverage_summary(all_estimated_results, all_true_results, "calibration_slope", "regression_slope")

write.table(correlation_coverage_df, file=file.path(viz_dir, "estimated_correlation_coverage_summary.txt"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(calibration_coverage_df, file=file.path(viz_dir, "estimated_calibration_slope_coverage_summary.txt"), sep="\t", quote=FALSE, row.names=FALSE)

correlation_coverage_plot <- make_coverage_plot(correlation_coverage_df, "Correlation Coverage")
calibration_coverage_plot <- make_coverage_plot(calibration_coverage_df, "Calibration Slope Coverage")

ggsave(filename=file.path(viz_dir, "estimated_correlation_coverage_plot.pdf"), plot=correlation_coverage_plot, width=7.2, height=3.6)
ggsave(filename=file.path(viz_dir, "estimated_calibration_slope_coverage_plot.pdf"), plot=calibration_coverage_plot, width=7.2, height=3.6)

print(file.path(viz_dir, "estimated_calibration_slope_bias_barplot.pdf"))
print(file.path(viz_dir, "estimated_calibration_slope_coverage_plot.pdf"))
