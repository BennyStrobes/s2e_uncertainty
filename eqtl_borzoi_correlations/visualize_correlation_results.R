args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
library(scales)
options(warn=1)

figure_theme <- function() {
	return(
		theme_cowplot(font_size=11) +
		theme(
			plot.title = element_text(face="bold", size=14, hjust=0),
			text = element_text(size=11),
			axis.text = element_text(size=10, color="#1F2937"),
			axis.title = element_text(size=11, face="bold"),
			panel.grid.major.y = element_line(color="#D1D5DB", linewidth=.35),
			panel.grid.major.x = element_blank(),
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line.x = element_line(colour="#111827"),
			axis.line.y = element_blank(),
			legend.position = "none"
		)
	)
}

make_compete_calibration_heatmap <- function(df) {
	all_tissues = sort(unique(c(df$output_tissue, df$input_tissue)))
	df$output_tissue = factor(df$output_tissue, levels=all_tissues)
	df$input_tissue = factor(df$input_tissue, levels=all_tissues)
	return(
		ggplot(df, aes(x=input_tissue, y=output_tissue, fill=calibration_slope)) +
			geom_tile(color="white", linewidth=.35) +
			geom_text(aes(label=sprintf("%.2f", calibration_slope)), size=2.6, color="#111827") +
			scale_fill_gradient2(
				low="#2F5D8C",
				mid="#F8F5F0",
				high="#B85C38",
				midpoint=0.0,
				labels=number_format(accuracy=.01)
			) +
			labs(
				title="Cross-Tissue Calibration Heatmap",
				x="Input tissue",
				y="Output tissue",
				fill="Calibration\nslope"
			) +
			figure_theme() +
			theme(
				legend.position="right",
				legend.title=element_text(face="bold"),
				axis.text.x=element_text(angle=45, hjust=1, vjust=1),
				axis.text.y=element_text(hjust=1),
				panel.grid=element_blank(),
				axis.line=element_blank(),
				plot.margin=margin(8, 14, 8, 8)
			)
	)
}

extract_cross_tissue_pair_from_file <- function(result_file) {
	file_name = basename(result_file)
	pair_name = sub("^ld_corr_results_", "", file_name)
	pair_name = sub("_intercept_bootstrap_stats\\.txt$", "", pair_name)
	matches = regmatches(
		pair_name,
		regexec("^(.+)_GTEX-[^_]+_(.+)_GTEX-[^_]+$", pair_name)
	)[[1]]
	if (length(matches) != 3) {
		stop(paste("Could not parse cross-tissue result file name:", file_name))
	}
	return(matches[2:3])
}

load_whole_blood_cross_tissue_correlations <- function(cross_tissue_ld_corr_results_output_dir) {
	result_files = list.files(
		cross_tissue_ld_corr_results_output_dir,
		pattern="^ld_corr_results_.*_intercept_bootstrap_stats\\.txt$",
		full.names=TRUE
	)
	if (length(result_files) == 0) {
		stop(paste("No cross-tissue intercept bootstrap stats files found in", cross_tissue_ld_corr_results_output_dir))
	}

	whole_blood_df = data.frame()
	for (result_file in result_files) {
		tissue_pair = extract_cross_tissue_pair_from_file(result_file)
		if ("Whole_Blood" %in% tissue_pair == FALSE) {
			next
		}
		other_tissue = tissue_pair[tissue_pair != "Whole_Blood"]
		if (length(other_tissue) != 1) {
			next
		}
		tmp_df = read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp_df = tmp_df[tmp_df$output_name == "correlation", ]
		if (nrow(tmp_df) == 0) {
			next
		}
		tmp_df$other_tissue = other_tissue
		tmp_df$result_file = basename(result_file)
		whole_blood_df = rbind(whole_blood_df, tmp_df)
	}
	if (nrow(whole_blood_df) == 0) {
		stop(paste("No Whole_Blood cross-tissue correlation rows were loaded from", cross_tissue_ld_corr_results_output_dir))
	}
	if (any(duplicated(whole_blood_df$other_tissue))) {
		duplicate_tissues = paste(unique(whole_blood_df$other_tissue[duplicated(whole_blood_df$other_tissue)]), collapse=", ")
		stop(paste("Multiple Whole_Blood cross-tissue correlation files found for:", duplicate_tissues))
	}
	return(whole_blood_df)
}


load_whole_blood_cross_tissue_high_mag_correlations <- function(cross_tissue_ld_corr_results_output_dir) {
	result_files = list.files(
		cross_tissue_ld_corr_results_output_dir,
		pattern="^ld_corr_results_.*_borzoi_magnitude_bins_bootstrap_stats\\.txt$",
		full.names=TRUE
	)
	if (length(result_files) == 0) {
		stop(paste("No cross-tissue intercept bootstrap stats files found in", cross_tissue_ld_corr_results_output_dir))
	}

	whole_blood_df = data.frame()
	for (result_file in result_files) {
		tissue_pair = extract_cross_tissue_pair_from_file(result_file)
		if ("Whole_Blood" %in% tissue_pair == FALSE) {
			next
		}
		other_tissue = tissue_pair[tissue_pair != "Whole_Blood"]
		if (length(other_tissue) != 1) {
			next
		}
		tmp_df = read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp_df = tmp_df[tmp_df$output_name == "correlation", ]
		if (nrow(tmp_df) == 0) {
			next
		}
		tmp_df$other_tissue = other_tissue
		tmp_df$result_file = basename(result_file)
		whole_blood_df = rbind(whole_blood_df, tmp_df)
	}
	if (nrow(whole_blood_df) == 0) {
		stop(paste("No Whole_Blood cross-tissue correlation rows were loaded from", cross_tissue_ld_corr_results_output_dir))
	}
	if (any(duplicated(whole_blood_df$other_tissue))) {
		duplicate_tissues = paste(unique(whole_blood_df$other_tissue[duplicated(whole_blood_df$other_tissue)]), collapse=", ")
		stop(paste("Multiple Whole_Blood cross-tissue correlation files found for:", duplicate_tissues))
	}
	return(whole_blood_df)
}



make_readable_tissue_labels <- function(tissue_names, wrap_width=16) {
	tissue_labels = gsub("_", " ", tissue_names)
	tissue_labels[tissue_names == "Adipose_Subcutaneous"] = "Adipose Sub."
	tissue_labels[tissue_names == "Skin_Not_Sun_Exposed_Suprapubic"] = "Skin (no sun)"
	return(sapply(tissue_labels, function(tissue_label) {
		paste(strwrap(tissue_label, width=wrap_width), collapse="\n")
	}))
}

make_whole_blood_cross_tissue_correlation_plot <- function(df) {
	plot_df = df[df$output_name == "correlation", ]
	plot_df$ci_lower = plot_df$mean - 1.96*plot_df$bootstrap_se
	plot_df$ci_upper = plot_df$mean + 1.96*plot_df$bootstrap_se
	plot_df = plot_df[order(plot_df$mean), ]
	plot_df$other_tissue = factor(plot_df$other_tissue, levels=plot_df$other_tissue)
	return(
		ggplot(plot_df, aes(x=other_tissue, y=mean)) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.18, linewidth=.5, color="#111827") +
			geom_point(size=2.6, shape=21, stroke=.45, fill="#2F7D5B", color="#111827") +
			geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
			scale_y_continuous(labels=number_format(accuracy=.01)) +
			scale_x_discrete(labels=function(x) make_readable_tissue_labels(x)) +
			labs(title="Whole Blood vs other") +
			xlab(NULL) +
					ylab("Cross-context\nCorrelation") +
					figure_theme() +
					theme(
						plot.title=element_text(size=11, hjust=.5),
						axis.text.x=element_text(size=7.5, angle=35, hjust=1, vjust=1, lineheight=.85),
						panel.grid.major.x=element_blank(),
						plot.margin=margin(5, 14, 0, 8)
					)
	)
}

make_bar_plot <- function(df, output_name, ylab) {
	plot_df = df[df$output_name == output_name, ]
	plot_df$annotation_name = factor(plot_df$annotation_name, levels=plot_df$annotation_name)
	plot_df$annotation_label = paste0("bin", seq_len(nrow(plot_df)))
	if (output_name == "calibration_slope") {
		plot_title = "Muscle Skeletal Calibration"
		bar_color = "#B85C38"
	} else {
		plot_title = "Muscle Skeletal Correlation"
		bar_color = "#5C7C4F"
	}
	return(
		ggplot(plot_df, aes(x=annotation_name, y=mean)) +
		geom_col(fill=bar_color, color="#111827", width=.72, linewidth=.25) +
		geom_errorbar(aes(ymin=empirical_ci_lower, ymax=empirical_ci_upper), width=.18, linewidth=.45, color="#111827") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		scale_x_discrete(labels=plot_df$annotation_label) +
		labs(title=plot_title) +
		xlab("Borzoi magnitude bin") +
		ylab(ylab) +
		figure_theme() +
		theme(
			axis.text.x=element_text(angle=35, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}

make_magnitude_correlation_plot <- function(df) {
	plot_df = df[df$output_name == "correlation", ]
	ordered_annotations = unique(plot_df$annotation_name)
	plot_df$annotation_name = factor(plot_df$annotation_name, levels=ordered_annotations)
	plot_df$annotation_label = c("<0.001", "0.001-0.01", "0.01-0.075", "0.075-0.2", ">=0.2")
	plot_df$ci_lower = plot_df$mean - 1.96*plot_df$bootstrap_se
	plot_df$ci_upper = plot_df$mean + 1.96*plot_df$bootstrap_se
	return(
		ggplot(plot_df, aes(x=annotation_name, y=mean)) +
			geom_col(fill="#3F88C5", color="#111827", width=.72, linewidth=.25) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.18, linewidth=.45, color="#111827") +
			geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
			scale_y_continuous(labels=number_format(accuracy=.01)) +
			scale_x_discrete(labels=plot_df$annotation_label) +
			labs(title="Borzoi magnitude strata") +
			xlab(NULL) +
			ylab("S-LDMC\nCorrelation") +
			figure_theme() +
			theme(
				plot.title=element_text(size=10, hjust=.5, margin=margin(b=1)),
				axis.text.x=element_text(size=8, angle=28, hjust=1, vjust=1),
				axis.title.y=element_text(size=9),
				axis.title.x=element_blank(),
				plot.margin=margin(1, 6, 1, 4)
			)
	)
}

make_distance_correlation_plot <- function(df) {
	plot_df = df[df$output_name == "correlation", ]
	ordered_annotations = unique(plot_df$annotation_name)
	plot_df$annotation_name = factor(plot_df$annotation_name, levels=ordered_annotations)
	plot_df$annotation_label = make_distance_bin_labels(ordered_annotations)
	plot_df$ci_lower = plot_df$mean - 1.96*plot_df$bootstrap_se
	plot_df$ci_upper = plot_df$mean + 1.96*plot_df$bootstrap_se
	return(
		ggplot(plot_df, aes(x=annotation_name, y=mean)) +
			geom_col(fill="#3F88C5", color="#111827", width=.72, linewidth=.25) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.18, linewidth=.45, color="#111827") +
			geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
			scale_y_continuous(labels=number_format(accuracy=.01)) +
			scale_x_discrete(labels=plot_df$annotation_label) +
			labs(title="Distance to TSS strata") +
			xlab(NULL) +
			ylab("S-LDMC\nCorrelation") +
			figure_theme() +
			theme(
				plot.title=element_text(size=10, hjust=.5, margin=margin(b=1)),
				axis.text.x=element_text(size=8, angle=25, hjust=1, vjust=1),
				axis.title.y=element_text(size=9),
				axis.title.x=element_blank(),
				plot.margin=margin(1, 6, 1, 4)
			)
	)
}

order_annotation_names <- function(annotation_names) {
	annotation_names = unique(annotation_names)
	parsed_ids = suppressWarnings(as.numeric(sub(".*[^0-9]([0-9]+)$", "\\1", annotation_names)))
	if (any(is.na(parsed_ids))) {
		return(sort(annotation_names))
	}
	return(annotation_names[order(parsed_ids)])
}

make_simulation_annotation_labels <- function(annotation_names) {
	parsed_ids = suppressWarnings(as.numeric(sub(".*[^0-9]([0-9]+)$", "\\1", annotation_names)))
	if (any(is.na(parsed_ids))) {
		return(annotation_names)
	}
	return(paste0("bin ", parsed_ids + 1))
}

extract_simulation_iter_from_file <- function(file_name) {
	matches = regmatches(file_name, regexpr("sim[0-9]+", file_name))
	if (length(matches) == 0 || matches[1] == "") {
		return(NA)
	}
	return(matches[1])
}

get_simulation_sem <- function(values) {
	values = values[is.finite(values)]
	if (length(values) <= 1) {
		return(0.0)
	}
	return(sd(values)/sqrt(length(values)))
}

get_true_simulation_results_dir <- function(simulation_results_dir) {
	clean_dir = sub("/+$", "", simulation_results_dir)
	return(file.path(dirname(clean_dir), "corr_sim_borzoi_eqtl_effects"))
}

get_xc_simulation_results_dir <- function(simulation_results_dir) {
	clean_dir = sub("/+$", "", simulation_results_dir)
	return(file.path(dirname(clean_dir), "xc_corr_sim_inference_results"))
}

get_xc_true_simulation_results_dir <- function(simulation_results_dir) {
	clean_dir = sub("/+$", "", simulation_results_dir)
	return(file.path(dirname(clean_dir), "xc_corr_sim_borzoi_eqtl_effects"))
}

extract_xc_simulation_rho_from_file <- function(file_name) {
	matches = regmatches(file_name, regexec("_rho_([^_]+)_", file_name))[[1]]
	if (length(matches) != 2) {
		return(NA)
	}
	return(matches[2])
}

load_simulation_correlation_summary <- function(simulation_results_dir) {
	result_files = list.files(
		simulation_results_dir,
		pattern="_(ld_corr|fm_corr)_results_bootstrap_stats\\.txt$",
		full.names=TRUE
	)
	if (length(result_files) == 0) {
		stop(paste("No LD-corr or fine-mapped simulation result files found in", simulation_results_dir))
	}

	all_est = data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp = read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp = tmp[tmp$output_name == "correlation", ]
		if (nrow(tmp) == 0) {
			next
		}
		file_name = basename(result_file)
		tmp$method = ifelse(grepl("_fm_corr_results_bootstrap_stats\\.txt$", file_name), "Fine-map", "S-LDMC")
		tmp$simulation_iter = extract_simulation_iter_from_file(file_name)
		all_est = rbind(all_est, tmp[, c("simulation_iter", "annotation_name", "method", "mean")])
	}
	if (nrow(all_est) == 0) {
		stop(paste("Simulation result files were found in", simulation_results_dir, "but no correlation rows were loaded"))
	}

	true_results_dir = get_true_simulation_results_dir(simulation_results_dir)
	true_files = list.files(true_results_dir, pattern="_true_sim_effect_summary\\.txt$", full.names=TRUE)
	if (length(true_files) == 0) {
		stop(paste("No true simulated summary files found in", true_results_dir))
	}

	all_true = data.frame()
	for (true_file in true_files) {
		if (is.na(file.info(true_file)$size) || file.info(true_file)$size == 0) {
			next
		}
		tmp = read.table(true_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp$simulation_iter = extract_simulation_iter_from_file(basename(true_file))
		tmp$method = "True simulated"
		tmp$mean = tmp$pearson_correlation
		all_true = rbind(all_true, tmp[, c("simulation_iter", "annotation_name", "method", "mean")])
	}
	if (nrow(all_true) == 0) {
		stop(paste("True simulated summary files were found in", true_results_dir, "but no rows were loaded"))
	}

	all_results = rbind(all_est, all_true)
	ordered_annotations = order_annotation_names(all_results$annotation_name)
	method_levels = c("Fine-map", "S-LDMC", "True simulated")
	summary_df = data.frame()
	for (method_name in method_levels) {
		for (annotation_name in ordered_annotations) {
			values = all_results$mean[all_results$method == method_name & all_results$annotation_name == annotation_name]
			values = values[is.finite(values)]
			if (length(values) == 0) {
				next
			}
			mean_value = mean(values)
			sem_value = get_simulation_sem(values)
			summary_df = rbind(
				summary_df,
				data.frame(
					annotation_name=annotation_name,
					method=method_name,
					mean=mean_value,
					ci_lower=mean_value - 1.96*sem_value,
					ci_upper=mean_value + 1.96*sem_value,
					n_sims=length(values),
					stringsAsFactors=FALSE
				)
			)
		}
	}
	summary_df$annotation_name = factor(summary_df$annotation_name, levels=ordered_annotations)
	summary_df$annotation_label = make_simulation_annotation_labels(as.character(summary_df$annotation_name))
	summary_df$method = factor(summary_df$method, levels=method_levels)
	rightmost_annotation = ordered_annotations[length(ordered_annotations)]
	rightmost_sldmc = summary_df$annotation_name == rightmost_annotation & summary_df$method == "S-LDMC"
	return(summary_df)
}

make_simulation_correlation_plot <- function(simulation_summary_df) {
	bar_df = simulation_summary_df[simulation_summary_df$method %in% c("Fine-map", "S-LDMC"), ]
	truth_df = simulation_summary_df[simulation_summary_df$method == "True simulated", ]
	return(
			ggplot(bar_df, aes(x=annotation_name, y=mean, fill=method)) +
			geom_col(position=position_dodge(width=.78), width=.68, color="#111827", linewidth=.25) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position=position_dodge(width=.78), width=.16, linewidth=.4, color="#111827") +
			geom_errorbar(data=truth_df, aes(x=annotation_name, ymin=mean, ymax=mean), inherit.aes=FALSE, width=.58, linewidth=.9, color="#6D3A9C") +
			geom_point(data=truth_df, aes(x=annotation_name, y=mean), inherit.aes=FALSE, shape=21, size=1.8, stroke=.55, fill="white", color="#6D3A9C") +
			geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
			scale_fill_manual(values=c("Fine-map"="#9CA3AF", "S-LDMC"="#3F88C5")) +
			scale_x_discrete(labels=make_simulation_annotation_labels(levels(simulation_summary_df$annotation_name))) +
			scale_y_continuous(labels=number_format(accuracy=.01)) +
			labs(title="Simulation correlation estimates", fill=NULL) +
			xlab("Annotation bin") +
			ylab("Correlation") +
			figure_theme() +
			theme(
				legend.position="top",
				plot.title=element_text(size=11, hjust=.5),
				axis.text.x=element_text(angle=0, hjust=.5),
				axis.title=element_text(size=10, face="bold"),
				plot.margin=margin(5, 8, 5, 5)
			)
		)
	}

load_xc_simulation_correlation_summary <- function(simulation_results_dir) {
	xc_simulation_results_dir = get_xc_simulation_results_dir(simulation_results_dir)
	result_files = list.files(
		xc_simulation_results_dir,
		pattern="_xc_(ld_corr|fm_corr)_results_bootstrap_stats\\.txt$",
		full.names=TRUE
	)
	if (length(result_files) == 0) {
		stop(paste("No cross-context LD-corr or fine-mapped simulation result files found in", xc_simulation_results_dir))
	}

	all_est = data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp = read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp = tmp[tmp$output_name == "correlation", ]
		file_name = basename(result_file)
		if (grepl("_xc_fm_corr_results_bootstrap_stats\\.txt$", file_name)) {
			tmp = tmp[tmp$estimator_name == "fine_mapped_posterior_mean", ]
			tmp$method = "Fine-map"
		} else {
			tmp$method = "XC-LDMC"
		}
		if (nrow(tmp) == 0) {
			next
		}
		tmp$simulation_iter = extract_simulation_iter_from_file(file_name)
		tmp$rho_delta = extract_xc_simulation_rho_from_file(file_name)
		all_est = rbind(all_est, tmp[, c("simulation_iter", "rho_delta", "method", "mean")])
	}
	if (nrow(all_est) == 0) {
		stop(paste("Cross-context simulation result files were found in", xc_simulation_results_dir, "but no correlation rows were loaded"))
	}

	xc_true_results_dir = get_xc_true_simulation_results_dir(simulation_results_dir)
	true_files = list.files(xc_true_results_dir, pattern="_true_sim_effect_summary\\.txt$", full.names=TRUE)
	if (length(true_files) == 0) {
		stop(paste("No cross-context true simulated summary files found in", xc_true_results_dir))
	}

	all_true = data.frame()
	for (true_file in true_files) {
		if (is.na(file.info(true_file)$size) || file.info(true_file)$size == 0) {
			next
		}
		tmp = read.table(true_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp = tmp[tmp$effect_type == "delta", ]
		if (nrow(tmp) == 0) {
			next
		}
		tmp$simulation_iter = extract_simulation_iter_from_file(basename(true_file))
		tmp$rho_delta = extract_xc_simulation_rho_from_file(basename(true_file))
		tmp$method = "True simulated"
		tmp$mean = tmp$pearson_correlation
		all_true = rbind(all_true, tmp[, c("simulation_iter", "rho_delta", "method", "mean")])
	}
	if (nrow(all_true) == 0) {
		stop(paste("Cross-context true simulated summary files were found in", xc_true_results_dir, "but no delta rows were loaded"))
	}

	all_results = rbind(all_est, all_true)
	rho_levels = unique(all_results$rho_delta[order(as.numeric(all_results$rho_delta))])
	method_levels = c("Fine-map", "XC-LDMC", "True simulated")
	summary_df = data.frame()
	for (method_name in method_levels) {
		for (rho_delta in rho_levels) {
			values = all_results$mean[all_results$method == method_name & all_results$rho_delta == rho_delta]
			values = values[is.finite(values)]
			if (length(values) == 0) {
				next
			}
			mean_value = mean(values)
			sem_value = get_simulation_sem(values)
			summary_df = rbind(
				summary_df,
				data.frame(
					rho_delta=rho_delta,
					method=method_name,
					mean=mean_value,
					ci_lower=mean_value - 1.96*sem_value,
					ci_upper=mean_value + 1.96*sem_value,
					n_sims=length(values),
					stringsAsFactors=FALSE
				)
			)
		}
	}
	summary_df$rho_delta = factor(summary_df$rho_delta, levels=rho_levels)
	summary_df$method = factor(summary_df$method, levels=method_levels)
	return(summary_df)
}

make_xc_simulation_correlation_plot <- function(simulation_summary_df) {
	bar_df = simulation_summary_df[simulation_summary_df$method %in% c("Fine-map", "XC-LDMC"), ]
	truth_df = simulation_summary_df[simulation_summary_df$method == "True simulated", ]
	return(
			ggplot(bar_df, aes(x=rho_delta, y=mean, fill=method)) +
			geom_col(position=position_dodge(width=.78), width=.68, color="#111827", linewidth=.25) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position=position_dodge(width=.78), width=.16, linewidth=.4, color="#111827") +
			geom_errorbar(data=truth_df, aes(x=rho_delta, ymin=mean, ymax=mean), inherit.aes=FALSE, width=.58, linewidth=.9, color="#6D3A9C") +
			geom_point(data=truth_df, aes(x=rho_delta, y=mean), inherit.aes=FALSE, shape=21, size=1.8, stroke=.55, fill="white", color="#6D3A9C") +
			geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
			scale_fill_manual(values=c("Fine-map"="#9CA3AF", "XC-LDMC"="#2F7D5B")) +
			scale_y_continuous(labels=number_format(accuracy=.01)) +
			labs(title="Cross-tissue simulation", fill=NULL) +
			xlab("Simulated cross-context correlation") +
			ylab("Cross-context\nCorrelation") +
			figure_theme() +
			theme(
				legend.position="top",
				plot.title=element_text(size=11, hjust=.5),
				axis.text.x=element_text(angle=0, hjust=.5),
				axis.title=element_text(size=10, face="bold"),
				plot.margin=margin(5, 8, 5, 5)
			)
	)
}

make_stratified_bar_plot <- function(df, output_name, ylab, plot_title, xlab_title, bar_color) {
	plot_df = df[df$output_name == output_name, ]
	plot_df$annotation_name = factor(plot_df$annotation_name, levels=plot_df$annotation_name)
	plot_df$annotation_label = paste0("bin", seq_len(nrow(plot_df)))
	return(
		ggplot(plot_df, aes(x=annotation_name, y=mean)) +
		geom_col(fill=bar_color, color="#111827", width=.72, linewidth=.25) +
		geom_errorbar(aes(ymin=empirical_ci_lower, ymax=empirical_ci_upper), width=.18, linewidth=.45, color="#111827") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		scale_x_discrete(labels=plot_df$annotation_label) +
		labs(title=plot_title) +
		xlab(xlab_title) +
		ylab(ylab) +
		figure_theme() +
		theme(
			axis.text.x=element_text(angle=35, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}

make_distance_bin_labels <- function(annotation_names) {
	default_labels = paste0("bin", seq_along(annotation_names))
	parsed_bin_ids = suppressWarnings(as.numeric(sub(".*bin([0-9]+)$", "\\1", annotation_names)))
	if (any(is.na(parsed_bin_ids))) {
		return(default_labels)
	}
	if (setequal(sort(unique(parsed_bin_ids)), c(0, 1, 2, 3, 4))) {
		label_map = c("0-1 kb", "1-5 kb", "5-25 kb", "25-50 kb", "50-100 kb")
		return(label_map[parsed_bin_ids + 1])
	}
	if (setequal(sort(unique(parsed_bin_ids)), c(0, 1, 2))) {
		label_map = c("0-5 kb", "5-50 kb", "50-100 kb")
		return(label_map[parsed_bin_ids + 1])
	}
	return(default_labels)
}

make_distance_bin_plot <- function(df, output_name, ylab, plot_title, bar_color) {
	plot_df = df[df$output_name == output_name, ]
	ordered_annotations = unique(plot_df$annotation_name)
	plot_df$annotation_name = factor(plot_df$annotation_name, levels=ordered_annotations)
	annotation_labels = make_distance_bin_labels(ordered_annotations)
	return(
		ggplot(plot_df, aes(x=annotation_name, y=mean)) +
		geom_col(fill=bar_color, color="#111827", width=.72, linewidth=.25) +
		geom_errorbar(aes(ymin=empirical_ci_lower, ymax=empirical_ci_upper), width=.18, linewidth=.45, color="#111827") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		scale_x_discrete(labels=annotation_labels) +
		labs(title=plot_title) +
		xlab("Distance to TSS bin") +
		ylab(ylab) +
		figure_theme() +
		theme(
			axis.text.x=element_text(angle=20, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}

make_variance_comparison_plot <- function(df) {
	plot_df = df[df$output_name %in% c("per_snp_eqtl_h2", "borzoi_variance"), ]
	ordered_annotations = unique(plot_df$annotation_name)
	plot_df$annotation_name = factor(plot_df$annotation_name, levels=ordered_annotations)
	annotation_labels = paste0("bin", seq_len(length(ordered_annotations)))
	plot_df$metric = factor(
		plot_df$output_name,
		levels=c("per_snp_eqtl_h2", "borzoi_variance"),
		labels=c("Per-SNP eQTL h2", "Per-SNP Borzoi variance")
	)
	return(
		ggplot(plot_df, aes(x=annotation_name, y=mean, fill=metric)) +
		geom_col(position=position_dodge(width=.78), width=.68, color="#111827", linewidth=.25) +
		geom_errorbar(aes(ymin=empirical_ci_lower, ymax=empirical_ci_upper), position=position_dodge(width=.78), width=.18, linewidth=.45, color="#111827") +
		scale_fill_manual(values=c("#3F88C5", "#7A4EAB")) +
		scale_x_discrete(labels=annotation_labels) +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		labs(title="Muscle Skeletal Variance Components", fill="Metric") +
		xlab("Borzoi magnitude bin") +
		ylab("Estimated value") +
		figure_theme() +
		theme(
			legend.position="top",
			legend.title=element_text(face="bold"),
			axis.text.x=element_text(angle=0, hjust=.5),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}

















#######################
# Command line args
#######################
ld_corr_results_output_dir = args[1]
simulation_results_dir = args[2]
cross_tissue_ld_corr_results_output_dir = args[3]
visualization_dir = args[4]

simulation_correlation_summary_df = load_simulation_correlation_summary(simulation_results_dir)
print(simulation_correlation_summary_df)
simulation_correlation_plot = make_simulation_correlation_plot(simulation_correlation_summary_df)
simulation_correlation_output_file = paste0(
	visualization_dir,
	"/simulation_correlation_estimates.pdf"
)
ggsave(simulation_correlation_output_file, simulation_correlation_plot, width=7.2, height=3.3)
print(simulation_correlation_output_file)


# Magnitude Analysis
target_tissue = "Whole_Blood"
target_sample = "GTEX-1LB8K-0005-SM-DIPED.1"
anno_method = "borzoi_magnitude_bins"
results_file = paste0(
	ld_corr_results_output_dir,
	"ld_corr_results_",
	target_tissue,
	"_",
	target_sample,
	"_",
	anno_method,
	"_bootstrap_stats.txt"
)
magnitude_results_df <- read.table(results_file, header=TRUE)



magnitude_correlation_plot = make_magnitude_correlation_plot(magnitude_results_df)
magnitude_correlation_output_file = paste0(
	visualization_dir,
	"/borzoi_magnitude_correlation.pdf"
)
ggsave(magnitude_correlation_output_file, magnitude_correlation_plot, width=7.2, height=3.3)
print(magnitude_correlation_output_file)





joint_plot <- plot_grid(
	simulation_correlation_plot,
	magnitude_correlation_plot,
	ncol=2,
	labels=c('a', 'b')
)

joint_output_file = paste0(
	visualization_dir,
	"/figure2.pdf"
)
ggsave(joint_output_file, joint_plot, width=7.4, height=2.25)

print(joint_output_file)


if (FALSE) {

xc_simulation_correlation_summary_df = load_xc_simulation_correlation_summary(simulation_results_dir)
print(xc_simulation_correlation_summary_df)
xc_simulation_correlation_plot = make_xc_simulation_correlation_plot(xc_simulation_correlation_summary_df)
xc_simulation_correlation_output_file = paste0(
	visualization_dir,
	"/xc_simulation_correlation_estimates.pdf"
)
ggsave(xc_simulation_correlation_output_file, xc_simulation_correlation_plot, width=7.2, height=3.3)
print(xc_simulation_correlation_output_file)

whole_blood_cross_tissue_df = load_whole_blood_cross_tissue_correlations(cross_tissue_ld_corr_results_output_dir)
print(whole_blood_cross_tissue_df)
whole_blood_cross_tissue_plot = make_whole_blood_cross_tissue_correlation_plot(whole_blood_cross_tissue_df)
whole_blood_cross_tissue_output_file = paste0(
	visualization_dir,
	"/whole_blood_cross_tissue_correlation.pdf"
)
joint <- plot_grid(xc_simulation_correlation_plot, whole_blood_cross_tissue_plot, ncol=2, labels=c("a", "b"))
ggsave(whole_blood_cross_tissue_output_file, joint, width=7.4, height=2.1)
print(whole_blood_cross_tissue_output_file)

}






if (FALSE) {

########################
# Load in data
########################
target_tissue = "Muscle_Skeletal"
target_sample = "GTEX-13QJ3-0726-SM-5SI68.1"
anno_method = "borzoi_magnitude_bins"

results_file = paste0(
	ld_corr_results_output_dir,
	"ld_corr_results_std_",
	target_tissue,
	"_",
	target_sample,
	"_",
	anno_method,
	"_bootstrap_stats.txt"
)

maf_anno_method = "af_bins"
maf_results_file = paste0(
	ld_corr_results_output_dir,
	"ld_corr_results_std_",
	target_tissue,
	"_",
	target_sample,
	"_",
	maf_anno_method,
	"_bootstrap_stats.txt"
)

dist_results_file = "/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/eqtl_borzoi_correlations_old/ld_corr_results/ld_corr_results_std_Muscle_Skeletal_GTEX-13QJ3-0726-SM-5SI68.1_dist_to_tss_bins_bootstrap_stats.txt"






################
# Make calibration coeficient and corelation bar plots
########################

if (file.exists(results_file) && file.exists(maf_results_file)) {
	results_df = read.table(results_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	results_df$annotation_name = factor(results_df$annotation_name, levels=unique(results_df$annotation_name))
	maf_results_df = read.table(maf_results_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	maf_results_df$annotation_name = factor(maf_results_df$annotation_name, levels=unique(maf_results_df$annotation_name))
	calibration_plot = make_bar_plot(results_df, "calibration_slope", "Calibration coefficient")
	correlation_plot = make_bar_plot(results_df, "correlation", "Correlation")
	variance_comparison_plot = make_variance_comparison_plot(results_df)
	maf_calibration_plot = make_stratified_bar_plot(maf_results_df, "calibration_slope", "Calibration coefficient", "Muscle Skeletal Calibration by MAF", "MAF bin", "#A44A3F")
	maf_correlation_plot = make_stratified_bar_plot(maf_results_df, "correlation", "Correlation", "Muscle Skeletal Correlation by MAF", "MAF bin", "#4F6D3A")

	calibration_output_file = paste0(
		visualization_dir,
		"/muscle_skeletal_calibration_coefficient.pdf"
	)
	ggsave(calibration_output_file, calibration_plot, width=7.2, height=3.3)

	correlation_output_file = paste0(
		visualization_dir,
		"/muscle_skeletal_correlation.pdf"
	)
	ggsave(correlation_output_file, correlation_plot, width=7.2, height=3.3)

	variance_output_file = paste0(
		visualization_dir,
		"/muscle_skeletal_variance_comparison.pdf"
	)
	ggsave(variance_output_file, variance_comparison_plot, width=7.2, height=3.6)

	maf_calibration_output_file = paste0(
		visualization_dir,
		"/muscle_skeletal_maf_calibration_coefficient.pdf"
	)
	ggsave(maf_calibration_output_file, maf_calibration_plot, width=7.2, height=3.3)

	maf_correlation_output_file = paste0(
		visualization_dir,
		"/muscle_skeletal_maf_correlation.pdf"
	)
	ggsave(maf_correlation_output_file, maf_correlation_plot, width=7.2, height=3.3)

	print(calibration_output_file)
	print(correlation_output_file)
	print(variance_output_file)
	print(maf_calibration_output_file)
	print(maf_correlation_output_file)
}

########################
# Make distance-to-TSS bar plots
########################
if (file.exists(dist_results_file)) {
	dist_results_df = read.table(dist_results_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dist_results_df$annotation_name = factor(dist_results_df$annotation_name, levels=unique(dist_results_df$annotation_name))
	dist_calibration_plot = make_distance_bin_plot(
		dist_results_df,
		"calibration_slope",
		"Calibration coefficient",
		"Muscle Skeletal Calibration by Distance to TSS",
		"#A44A3F"
	)
	dist_correlation_plot = make_distance_bin_plot(
		dist_results_df,
		"correlation",
		"Correlation",
		"Muscle Skeletal Correlation by Distance to TSS",
		"#4F6D3A"
	)

	dist_calibration_output_file = paste0(
		visualization_dir,
		"/muscle_skeletal_distance_to_tss_calibration_coefficient.pdf"
	)
	ggsave(dist_calibration_output_file, dist_calibration_plot, width=7.2, height=3.3)

	dist_correlation_output_file = paste0(
		visualization_dir,
		"/muscle_skeletal_distance_to_tss_correlation.pdf"
	)
	ggsave(dist_correlation_output_file, dist_correlation_plot, width=7.2, height=3.3)

	print(dist_calibration_output_file)
	print(dist_correlation_output_file)
}

########################
# Compete heatmap
########################
compete_result_files = list.files(
	ld_corr_results_output_dir,
	pattern="^ld_corr_compete_results_.*_bootstrap_stats\\.txt$",
	full.names=TRUE
)
compete_result_files = compete_result_files[grepl("_eqtl_variance_bootstrap_stats\\.txt$", compete_result_files) == FALSE]

if (length(compete_result_files) == 0) {
	stop("No compete bootstrap stats files found")
}

compete_dfs = lapply(compete_result_files, function(compete_file) {
	file_name = basename(compete_file)
	output_tissue = sub("^ld_corr_compete_results_(.*)_bootstrap_stats\\.txt$", "\\1", file_name)
	tmp_df = read.table(compete_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	track_col = "track_name"
	if ("track_name" %in% colnames(tmp_df) == FALSE) {
		if ("annotation_name" %in% colnames(tmp_df)) {
			track_col = "annotation_name"
		} else {
			stop(paste("Could not find track column in", compete_file))
		}
	}
	tmp_df = tmp_df[tmp_df$output_name == "calibration_slope", c(track_col, "mean")]
	colnames(tmp_df) = c("input_tissue", "calibration_slope")
	tmp_df$output_tissue = output_tissue
	return(tmp_df[, c("output_tissue", "input_tissue", "calibration_slope")])
})

compete_heatmap_df = do.call(rbind, compete_dfs)
compete_heatmap_plot = make_compete_calibration_heatmap(compete_heatmap_df)
compete_heatmap_output_file = paste0(
	visualization_dir,
	"/cross_tissue_calibration_heatmap.pdf"
)

n_tissues = length(sort(unique(c(compete_heatmap_df$output_tissue, compete_heatmap_df$input_tissue))))
#plot_width = max(7.5, 0.42*n_tissues + 3.5)
#plot_height = max(6.5, 0.42*n_tissues + 2.5)
joint <- plot_grid(NULL, compete_heatmap_plot, ncol=2, labels=c("A", "B"))
ggsave(compete_heatmap_output_file, joint, width=7.4, height=3.0)

}
