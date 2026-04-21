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
visualization_dir = args[2]


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
plot_width = max(7.5, 0.42*n_tissues + 3.5)
plot_height = max(6.5, 0.42*n_tissues + 2.5)
ggsave(compete_heatmap_output_file, compete_heatmap_plot, width=plot_width, height=plot_height)
print(compete_heatmap_output_file)
