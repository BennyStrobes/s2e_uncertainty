args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


pretty_tissue <- function(x) {
  x <- as.character(x)

  # Normalize weird dash characters just in case (e.g., "EBV−transformed")
  x <- gsub("\u2212|\u2013|\u2014", "-", x, perl = TRUE)

  # Fix a couple likely truncations / typos (edit as needed)
  #x[x == "Whole_Bloo"] <- "Whole_Blood"

  # Hand-tuned short labels for the long ones
  map <- c(
    "Heart_Atrial_Appendage"              = "Heart (Atrial)",
    "Heart_Left_Ventricle"                = "Heart (LV)",
    "Adipose_Subcutaneous"                = "Adipose (Subcut)",
    "Adipose_Visceral_Omentum"            = "Adipose (Visceral)",
    "Artery_Aorta"                        = "Artery (Aorta)",
    "Artery_Tibial"                       = "Artery (Tibial)",
    "Nerve_Tibial"                        = "Nerve (Tibial)",
    "Esophagus_Gastroesophageal_Junction" = "Esophagus (GEJ)",
    "Esophagus_Mucosa"                    = "Esophagus (Mucosa)",
    "Colon_Transverse"                    = "Colon (Transverse)",
    "Colon_Sigmoid"                       = "Colon (Sigmoid)",
    "Small_Intestine_Terminal_Ileum"      = "Small Intestine (Ileum)",
    "Breast_Mammary_Tissue"               = "Breast (Mammary)",
    "Kidney_Cortex"                       = "Kidney (Cortex)",
    "Muscle_Skeletal"                     = "Muscle (Skeletal)",
    "Minor_Salivary_Gland"                = "Salivary (Minor)",
    "Skin_Not_Sun_Exposed_Suprapubic"     = "Skin (No Sun)",
    "Skin_Sun_Exposed_Lower_leg"          = "Skin (Sun)",
    "Cells_Cultured_fibroblasts"          = "Cells (Fibroblasts)",
    "Cells_EBV-transformed_lymphocytes"   = "Cells (EBV Lymphs)",
    "Brain_Cerebellar_Hemisphere" = "Brain (CH)"
  )

  out <- x
  hit <- x %in% names(map)
  out[hit] <- unname(map[x[hit]])

  # Fallback: underscores -> spaces, Title Case-ish, plus a few nicer words
  miss <- !hit
  if (any(miss)) {
    y <- x[miss]
    y <- gsub("_", " ", y, fixed = TRUE)

    # light cleanup
    y <- gsub("\\bGEJ\\b", "GEJ", y)
    y <- gsub("\\bLV\\b", "LV", y)
    y <- gsub("\\bEBV\\b", "EBV", y)

    # Title Case without extra dependencies
    y <- tolower(y)
    y <- gsub("(^|\\s)([a-z])", "\\1\\U\\2", y, perl = TRUE)

    out[miss] <- y
  }

  out
}

plot_jaccard_heatmap <- function(df,
                                 tissue1_col = "tissue1",
                                 tissue2_col = "tissue2",
                                 value_col   = "jaccard_index",
                                 method = "average",
                                 show_labels = TRUE) {

  stopifnot(all(c(tissue1_col, tissue2_col, value_col) %in% names(df)))

  # Pull vectors
  t1 <- df[[tissue1_col]]
  t2 <- df[[tissue2_col]]
  v  <- df[[value_col]]

  tissues <- sort(unique(c(t1, t2)))
  n <- length(tissues)

  # Symmetric similarity matrix
  S <- matrix(NA_real_, n, n, dimnames = list(tissues, tissues))

  # Fill both directions
  idx1 <- cbind(match(t1, tissues), match(t2, tissues))
  idx2 <- cbind(match(t2, tissues), match(t1, tissues))
  S[idx1] <- v
  S[idx2] <- v

  # Check completeness (ignore diagonal)
  if (any(is.na(S))) {
    na_pairs <- which(is.na(S), arr.ind = TRUE)
    na_pairs <- na_pairs[na_pairs[,1] != na_pairs[,2], , drop = FALSE]
    if (nrow(na_pairs) > 0) {
      stop("Similarity matrix has missing entries (NA). Are all tissue pairs present?")
    }
  }

  # ---- diagonal = max off-diagonal ----
  offdiag <- S[upper.tri(S) | lower.tri(S)]
  max_offdiag <- max(offdiag, na.rm = TRUE)
  diag(S) <- max_offdiag
  # -----------------------------------

  # Cluster by distance = 1 - similarity
  D <- as.dist(1 - S)
  hc <- hclust(D, method = method)
  ord <- hc$labels[hc$order]

  # Long format for ggplot
  plot_df <- expand.grid(tissue1 = ord, tissue2 = ord, KEEP.OUT.ATTRS = FALSE)

  ri <- match(plot_df$tissue1, rownames(S))
  ci <- match(plot_df$tissue2, colnames(S))
  plot_df$value <- S[cbind(ri, ci)]

  # enforce identical ordering on axes
  plot_df$tissue1 <- factor(plot_df$tissue1, levels = ord)
  plot_df$tissue2 <- factor(plot_df$tissue2, levels = rev(ord))


lab_map <- setNames(pretty_tissue(ord), ord)



  # Plot
  p <- ggplot(plot_df, aes(x = tissue1, y = tissue2, fill = value)) +
    geom_tile() +
    coord_fixed() +
  scale_x_discrete(labels = lab_map) +
  scale_y_discrete(labels = rev(lab_map)) +
    scale_fill_viridis_c(
      option="viridis",
      limits = c(0, max_offdiag),
      name   = "Jaccard",
      oob    = scales::squish
    ) +

    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank()
    )

  p <- p + figure_theme()
  p <- p + theme(legend.position = "bottom")

p <- p + guides(
  fill = guide_colorbar(
    barheight = unit(3, "mm"),
    barwidth  = unit(40, "mm")
  )
)+
  theme(
    legend.margin = margin(t = -1, b = -1),      # trims space around legend
    legend.box.margin = margin(t = -1, b = -1),  # trims space around legend box
    legend.spacing.y = unit(0, "mm")             # removes internal vertical spacing
  )
  p <- p + theme(
  legend.title = element_text(size = 9),
  legend.text  = element_text(size = 8)
  )

  if (show_labels) {
    p <- p + theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_blank()
    )
  } else {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank()
    )
  }

  return(p)
}


load_in_borzoi_gtex_enrichments <- function(tissues, variant_annotation_enrichment_dir) {
	borzoi_tissues <- c()
	gtex_tissues <- c()
	enrichments <- c()

	for (tiss_iter in 1:length(tissues)) {
		tissue_name <- tissues[tiss_iter]
		filer <- paste0(variant_annotation_enrichment_dir, "variant_enrichment_in_fine_mapped_eqtls_",tissue_name, "_0.2_0.5_summary.txt")
		tmp_df <- read.table(filer,header=TRUE)

		tmp_df <- tmp_df[as.character(tmp_df$annotation_name) !="Bladder",]


		gtex_tissues <- c(gtex_tissues, as.character(tmp_df$annotation_name))
		enrichments <- c(enrichments, tmp_df$enrichment)
		borzoi_tissues <- c(borzoi_tissues, rep(tissue_name, length(tmp_df$annotation_name)))
	}

	df <- data.frame(
  		borzoi_tissue  = as.character(borzoi_tissues),
  		finemap_tissue = as.character(gtex_tissues),
  		enrichment     = enrichments
	)

	return(df)
}


load_in_borzoi_gtex_delta_enrichments <- function(tissues, variant_annotation_enrichment_dir) {
	tissue_name_arr <- c()
	delta_enrichment <- c()

	for (tiss_iter in 1:length(tissues)) {
		tissue_name <- tissues[tiss_iter]
		filer <- paste0(variant_annotation_enrichment_dir, "variant_enrichment_in_fine_mapped_eqtls_",tissue_name, "_0.15_0.9_summary.txt")
		tmp_df <- read.table(filer,header=TRUE)

		tmp_df <- tmp_df[as.character(tmp_df$annotation_name) %in% tissues,]

		print(tissue_name)
		print(tmp_df)


		tissue_name_arr <- c(tissue_name_arr, tissue_name)


		matched_val = tmp_df$enrichment[as.character(tmp_df$annotation_name) == tissue_name]

		non_matched_val = median(tmp_df$enrichment[as.character(tmp_df$annotation_name) != tissue_name])

		delta_val = matched_val - non_matched_val

		delta_enrichment <- c(delta_enrichment, delta_val)

	}

	df <- data.frame(
  		tissue  = as.character(tissue_name_arr),
  		delta_enrichment     = delta_enrichment
	)

	return(df)
}


make_dist_to_tss_plot_v2 <- function(dist_df) {

  eqtl_indices  <- (dist_df$method == "eQTL")   & (dist_df$significance == 1.0)
  borzoi_indices <- (dist_df$method == "borzoi") & (dist_df$significance <= 0.05)

  df <- dist_df[eqtl_indices | borzoi_indices, , drop = FALSE]
  df <- df[abs(df$dist_to_tss) < 100000, , drop = FALSE]

  # Guard: density needs >=2 points with some variation
  if (nrow(df) < 2 || length(unique(df$dist_to_tss)) < 2) {
    warning("Not enough data in the selected window to plot a density.")
    return(NULL)
  }

pp <- ggplot(df, aes(x = dist_to_tss)) +

  # eQTL: red, line only (no fill)
  geom_density(
    data = df[df$method == "eQTL", ],
    aes(color = method),
    adjust = 1,
    linewidth = 1.0
  ) +

  # Borzoi: purple, filled
  geom_density(
    data = df[df$method == "borzoi", ],
    aes(color = method, fill = method),
    adjust = 1,
    alpha = 0.3,
    linewidth = 1.3
  ) +

  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +

  scale_x_continuous(
    limits = c(-100000, 100000),
    breaks = seq(-100000, 100000, by = 20000),
    labels = function(x) x / 1000
  ) +

  scale_color_manual(values = c("eQTL" = "#4C5C68", "borzoi" = "darkorchid2")) +
  scale_fill_manual(values = c("borzoi" = "darkorchid2")) +

  labs(
    x = "Distance to TSS (kb, signed)",
    y = "Density",
    color = "",
    fill = ""
  ) +

  figure_theme() +
  theme(legend.position = "bottom")+
  guides(fill = "none")
}



make_dist_to_tss_plot <- function(dist_df, fill=TRUE) {
  eqtl_indices  <- (dist_df$method == "eQTL")   & (dist_df$significance == 1.0)
  borzoi_indices <- (dist_df$method == "borzoi") & (dist_df$significance <= 0.05)

  df <- dist_df[eqtl_indices | borzoi_indices, , drop = FALSE]
  df <- df[abs(df$dist_to_tss) < 100000, , drop = FALSE]

  # Guard: density needs >=2 points with some variation
  if (nrow(df) < 2 || length(unique(df$dist_to_tss)) < 2) {
    warning("Not enough data in the selected window to plot a density.")
    return(NULL)
  }

  if (fill == TRUE) {
  pp <- ggplot(df, aes(x = dist_to_tss, color = method, fill = method)) +
    geom_density(alpha = 0.3, adjust = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_continuous(
      limits = c(-100000, 100000),
      breaks = seq(-100000, 100000, by = 20000),
      labels = function(x) x / 1000
    ) +
    scale_fill_manual(values = c("eQTL" = "red", "borzoi" = "darkorchid2")) +
    scale_color_manual(values = c("eQTL" = "red", "borzoi" = "darkorchid2")) +
    labs(
      x = "Distance to TSS (kb, signed)",
      y = "Density",
      color = "",
      fill = ""
    ) +
    figure_theme() +
    theme(legend.position="bottom")
  } else {
  pp <- ggplot(df, aes(x = dist_to_tss, color = method)) +
    geom_density(alpha = 0.3, adjust = 1, linewidth = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_continuous(
      limits = c(-100000, 100000),
      breaks = seq(-100000, 100000, by = 20000),
      labels = function(x) x / 1000
    ) +
    scale_color_manual(values = c("eQTL" = "red", "borzoi" = "darkorchid2")) +
    labs(
      x = "Distance to TSS (kb, signed)",
      y = "Density",
      color = "",
      fill = ""
    ) +
    figure_theme() +
    theme(legend.position="bottom")
  }

  return(pp)
}


######################
# Command line args
######################
visualize_variant_anno_enrich_dir = args[1]
variant_annotation_enrichment_dir = args[2]
gtex_tissue_names_file = args[3]



# Tissue overlap heatmap between sig borzoi effects and sig borzoi effects
if (FALSE) {
tissue_overlap_df <- read.table(paste0(variant_annotation_enrichment_dir, "tisue_overlaps_0.1.txt"), header=TRUE, sep="\t")
output_file <- paste0(visualize_variant_anno_enrich_dir, "borzoi_sig_effect_jacard_heatmap.pdf")
pp <- plot_jaccard_heatmap(tissue_overlap_df)
ggsave(pp, file=output_file, width=5.2, height=5.0, units="in")
}



tissue_name="Whole_Blood"
dist_file=paste0(variant_annotation_enrichment_dir,"dist_to_tss_summary_", tissue_name, ".txt")
dist_df <- read.table(dist_file, header=TRUE, sep="\t")

pp <- make_dist_to_tss_plot_v2(dist_df)
output_file <- paste0(visualize_variant_anno_enrich_dir, "dist_around_tss_density.pdf")
ggsave(pp, file=output_file, width=5.2, height=3.8, units="in")
print(output_file)



# load in gtex tissues
if (FALSE) {
tissues <- as.character(read.table(gtex_tissue_names_file, header=TRUE)$tissue)
tissues <- tissues[tissues != "Bladder"]

filer <- paste0(variant_annotation_enrichment_dir, "variant_enrichment_in_fine_mapped_eqtls_Brain_Amygdala_0.2_0.5_summary.txt")
tmp_df <- read.table(filer,header=TRUE, sep="\t")
powered_tissues_indices = (tmp_df$aa + tmp_df$cc) > 70
powered_tissues <- as.character(tmp_df$annotation_name[powered_tissues_indices])
}

if (FALSE) {
borzoi_gtex_enrichments_df <- load_in_borzoi_gtex_enrichments(tissues, variant_annotation_enrichment_dir)
matched_tissue_indices = borzoi_gtex_enrichments_df$borzoi_tissue == borzoi_gtex_enrichments_df$finemap_tissue
non_matched_tissue_indices = borzoi_gtex_enrichments_df$borzoi_tissue != borzoi_gtex_enrichments_df$finemap_tissue
}


if (FALSE) {
borzoi_gtex_enrichments_df <- load_in_borzoi_gtex_delta_enrichments(powered_tissues, variant_annotation_enrichment_dir)
}





