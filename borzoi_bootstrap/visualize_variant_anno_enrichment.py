import numpy as np
import os
import sys
import pdb
import matplotlib.pyplot as plt





def extract_non_redundant_gwas_trait_names(non_redundant_gwas_traits_file):
	f = open(non_redundant_gwas_traits_file)
	head_count = 0
	traits = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		traits.append(data[0])
	f.close()
	return np.asarray(traits)

def parse_sldsc_results_for_enrichments(trait_sldsc_result_file, non_binary_annos):
	f = open(trait_sldsc_result_file)
	enrichments = []
	enrichment_ses = []
	anno_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'baseL2_0':
			continue
		line_anno_name = data[0].split('L2')[0]
		if line_anno_name in non_binary_annos:
			continue
		enrichments.append(float(data[4]))
		enrichment_ses.append(float(data[5]))
		anno_names.append(line_anno_name)

	f.close()


	return np.asarray(enrichments), np.asarray(enrichment_ses), np.asarray(anno_names)



def parse_sldsc_results_for_enrichments_v2(trait_sldsc_result_file, anno_subset):
	f = open(trait_sldsc_result_file)
	enrichments = []
	enrichment_ses = []
	anno_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'baseL2_0':
			continue
		line_anno_name = data[0].split('L2')[0]
		if line_anno_name not in anno_subset:
			continue
		enrichments.append(float(data[4]))
		enrichment_ses.append(float(data[5]))
		anno_names.append(anno_subset[line_anno_name])

	f.close()


	return np.asarray(enrichments), np.asarray(enrichment_ses), np.asarray(anno_names)


def meta_analyze_enrichment(enrichment_mat, enrichment_ses_mat, min_valid=1):
	"""
	Fixed-effect inverse-variance meta-analysis across samples (rows) for each feature (column).

	Parameters
	----------
	enrichment_mat : array-like, shape (N, K)
		Per-sample enrichment estimates.
	enrichment_ses_mat : array-like, shape (N, K)
		Per-sample standard errors for enrichment estimates.
	min_valid : int
		Minimum number of valid (estimate, SE) pairs required per feature to return a value.
		If fewer than min_valid are valid, returns np.nan for that feature.

	Returns
	-------
	meta_enrichment : np.ndarray, shape (K,)
		Meta-analyzed enrichment estimate for each feature.
	meta_se : np.ndarray, shape (K,)
		Standard error for the meta-analyzed estimate for each feature.
	n_valid : np.ndarray, shape (K,)
		Number of valid samples contributing to each feature.
	"""
	E = np.asarray(enrichment_mat, dtype=float)
	S = np.asarray(enrichment_ses_mat, dtype=float)

	if E.shape != S.shape:
		raise ValueError(f"Shape mismatch: enrichment_mat {E.shape} vs enrichment_ses_mat {S.shape}")

	# Valid if estimate is finite AND SE is finite and > 0
	valid = np.isfinite(E) & np.isfinite(S) & (S > 0)

	# Inverse-variance weights; invalid entries get weight 0
	W = np.zeros_like(S)
	W[valid] = 1.0 / (S[valid] ** 2)

	# Weighted sums per feature (column)
	sumW = np.sum(W, axis=0)              # (K,)
	sumWE = np.sum(W * np.where(valid, E, 0.0), axis=0)  # (K,)

	meta_enrichment = np.full(E.shape[1], np.nan)
	meta_se = np.full(E.shape[1], np.nan)

	n_valid = np.sum(valid, axis=0)

	ok = (n_valid >= min_valid) & (sumW > 0)
	meta_enrichment[ok] = sumWE[ok] / sumW[ok]
	meta_se[ok] = np.sqrt(1.0 / sumW[ok])

	return meta_enrichment, meta_se, n_valid

def make_meta_analyzed_gwas_enrichment_file(meta_gwas_enrichment_file, non_redundant_gwas_traits_file, gwas_sldsc_results_dir, non_binary_annos):
	traits = extract_non_redundant_gwas_trait_names(non_redundant_gwas_traits_file)

	# Now loop through traits and extract enrichments and enrichment_ses
	enrichment_mat = []
	enrichment_ses_mat = []
	for trait in traits:

		trait_sldsc_result_file = gwas_sldsc_results_dir + trait + '_sldsc_baselineLD_v2.2.results'

		enrichment_vec, enrichment_se_vec, annotation_names_vec = parse_sldsc_results_for_enrichments(trait_sldsc_result_file, non_binary_annos)

		# add to global
		enrichment_mat.append(enrichment_vec)
		enrichment_ses_mat.append(enrichment_se_vec)

	# reorg global
	enrichment_mat = np.asarray(enrichment_mat)
	enrichment_ses_mat = np.asarray(enrichment_ses_mat)

	meta_enrichment, meta_enrichment_se, tmp_nvalid = meta_analyze_enrichment(enrichment_mat, enrichment_ses_mat)

	t = open(meta_gwas_enrichment_file,'w')
	t.write('anno_name\tmeta_enrichment\tmeta_enrichment_se\n')

	for ii, anno_name in enumerate(annotation_names_vec):
		t.write(anno_name + '\t' + str(meta_enrichment[ii]) + '\t' + str(meta_enrichment_se[ii]) + '\n')

	return meta_enrichment, meta_enrichment - (1.96*meta_enrichment_se),  meta_enrichment + (1.96*meta_enrichment_se), annotation_names_vec



def make_meta_analyzed_gwas_enrichment_file_for_specific_annos(meta_gwas_enrichment_file, non_redundant_gwas_traits_file, gwas_sldsc_results_dir, anno_subset):
	traits = extract_non_redundant_gwas_trait_names(non_redundant_gwas_traits_file)

	# Now loop through traits and extract enrichments and enrichment_ses
	enrichment_mat = []
	enrichment_ses_mat = []
	for trait in traits:

		trait_sldsc_result_file = gwas_sldsc_results_dir + trait + '_sldsc_baselineLD_v2.2.results'

		enrichment_vec, enrichment_se_vec, annotation_names_vec = parse_sldsc_results_for_enrichments_v2(trait_sldsc_result_file, anno_subset)

		# add to global
		enrichment_mat.append(enrichment_vec)
		enrichment_ses_mat.append(enrichment_se_vec)

	# reorg global
	enrichment_mat = np.asarray(enrichment_mat)
	enrichment_ses_mat = np.asarray(enrichment_ses_mat)

	meta_enrichment, meta_enrichment_se, tmp_nvalid = meta_analyze_enrichment(enrichment_mat, enrichment_ses_mat)

	t = open(meta_gwas_enrichment_file,'w')
	t.write('anno_name\tmeta_enrichment\tmeta_enrichment_se\n')

	for ii, anno_name in enumerate(annotation_names_vec):
		t.write(anno_name + '\t' + str(meta_enrichment[ii]) + '\t' + str(meta_enrichment_se[ii]) + '\n')

	return meta_enrichment, meta_enrichment - (1.96*meta_enrichment_se),  meta_enrichment + (1.96*meta_enrichment_se), annotation_names_vec


def get_borzoi_anno_enrichments(enrichment_summary_file, non_binary_annos, filter_to_binary_anno):
	f = open(enrichment_summary_file)
	enrichment_vec = []
	enrichment_lb_vec = []
	enrichment_ub_vec = []
	enrichment_p_vec = []

	anno_name_vec = []
	head_count = 0
	f = open(enrichment_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_anno_name = data[0]

		if filter_to_binary_anno and line_anno_name in non_binary_annos:
			continue
		anno_name_vec.append(data[0])
		enrichment_vec.append(float(data[1]))
		enrichment_lb_vec.append(float(data[3]))
		enrichment_ub_vec.append(float(data[4]))
		enrichment_p_vec.append(float(data[2]))



	f.close()



	return np.asarray(enrichment_vec), np.asarray(enrichment_lb_vec), np.asarray(enrichment_ub_vec), np.asarray(anno_name_vec), np.asarray(enrichment_p_vec)

def get_borzoi_anno_enrichments_for_specific_annos(enrichment_summary_file, anno_subset):
	f = open(enrichment_summary_file)
	enrichment_vec = []
	enrichment_lb_vec = []
	enrichment_ub_vec = []
	enrichment_p_vec = []

	anno_name_vec = []
	head_count = 0
	f = open(enrichment_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_anno_name = data[0]

		if line_anno_name not in anno_subset:
			continue
		revised_anno_name = anno_subset[line_anno_name]
		anno_name_vec.append(revised_anno_name)
		enrichment_vec.append(float(data[1]))
		enrichment_lb_vec.append(float(data[3]))
		enrichment_ub_vec.append(float(data[4]))
		enrichment_p_vec.append(float(data[2]))



	f.close()



	return np.asarray(enrichment_vec), np.asarray(enrichment_lb_vec), np.asarray(enrichment_ub_vec), np.asarray(anno_name_vec), np.asarray(enrichment_p_vec)


def plot_enrichment_scatter_with_ci_log(
	meta_enrichment,
	meta_enrichment_lb,
	meta_enrichment_ub,
	borzoi_enrichment,
	borzoi_enrichment_lb,
	borzoi_enrichment_ub,
	outpath,
	*,
	title=None,
	xlabel="Meta enrichment",
	ylabel="Borzoi enrichment",
	show_identity=True,
	annotate=None,
	point_alpha=0.9,
	capsize=2.5,
	figsize=(6.5, 6.5),
	dpi=200,
):
	import numpy as np
	import matplotlib.pyplot as plt

	x   = np.asarray(meta_enrichment, dtype=float)
	xlb = np.asarray(meta_enrichment_lb, dtype=float)
	xub = np.asarray(meta_enrichment_ub, dtype=float)

	y   = np.asarray(borzoi_enrichment, dtype=float)
	ylb = np.asarray(borzoi_enrichment_lb, dtype=float)
	yub = np.asarray(borzoi_enrichment_ub, dtype=float)

	n = len(x)
	for arr, name in [
		(xlb, "meta_enrichment_lb"), (xub, "meta_enrichment_ub"),
		(y, "borzoi_enrichment"), (ylb, "borzoi_enrichment_lb"), (yub, "borzoi_enrichment_ub")
	]:
		if len(arr) != n:
			raise ValueError(f"{name} has length {len(arr)} vs {n}")

	# --- Log-scale requires everything > 0 ---
	mask = (
		np.isfinite(x) & np.isfinite(xlb) & np.isfinite(xub) &
		np.isfinite(y) & np.isfinite(ylb) & np.isfinite(yub) &
		(x > 0) & (xlb > 0) &
		(y > 0) & (ylb > 0)
	)

	x, xlb, xub = x[mask], xlb[mask], xub[mask]
	y, ylb, yub = y[mask], ylb[mask], yub[mask]

	xlb[xlb < .2] = .2

	# Error sizes (linear)
	xerr = np.vstack([x - xlb, xub - x])
	yerr = np.vstack([y - ylb, yub - y])

	if np.any(xerr < 0) or np.any(yerr < 0):
		raise ValueError("Some CI bounds appear inverted.")

	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

	ax.errorbar(
		x, y,
		xerr=xerr,
		yerr=yerr,
		fmt="o",
		capsize=capsize,
		elinewidth=1.0,
		markersize=5,
		alpha=point_alpha,
	)

	# --- LOG–LOG ---
	ax.set_xscale("log")
	ax.set_yscale("log")

	if annotate is not None:
		labels = np.asarray(annotate)[mask]
		for xi, yi, lab in zip(x, y, labels):
			ax.annotate(str(lab), (xi, yi),
						xytext=(4, 4),
						textcoords="offset points",
						fontsize=8)

	if show_identity:
		lo = max(np.min(x), np.min(y))
		hi = min(np.max(x), np.max(y))
		if lo > 0:
			ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0)

	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if title is not None:
		ax.set_title(title)

	ax.grid(True, which="both", linewidth=0.5, alpha=0.4)
	fig.tight_layout()
	fig.savefig(outpath, bbox_inches="tight")
	plt.close(fig)

def plot_enrichment_scatter_with_ci(
	meta_enrichment,
	meta_enrichment_lb,
	meta_enrichment_ub,
	borzoi_enrichment,
	borzoi_enrichment_lb,
	borzoi_enrichment_ub,
	outpath,
	*,
	title=None,
	xlabel="Meta enrichment",
	ylabel="Borzoi enrichment",
	show_identity=True,
	annotate=None,          # list of labels, same length as vectors (optional)
	point_alpha=0.9,
	capsize=2.5,
	figsize=(6.5, 6.5),
	dpi=200,
):
	"""
	Scatterplot comparing enrichments with 95% CI on both x and y axes.

	Parameters
	----------
	meta_enrichment, meta_enrichment_lb, meta_enrichment_ub : array-like
		Point estimate and lower/upper bounds for meta enrichment (x-axis).
	borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub : array-like
		Point estimate and lower/upper bounds for borzoi enrichment (y-axis).
	outpath : str
		Output filepath (e.g. "plot.png" or "plot.pdf").
	show_identity : bool
		If True, draws y=x reference line.
	annotate : array-like of str or None
		Optional labels per point.
	"""

	x = np.asarray(meta_enrichment, dtype=float)
	xlb = np.asarray(meta_enrichment_lb, dtype=float)
	xub = np.asarray(meta_enrichment_ub, dtype=float)

	y = np.asarray(borzoi_enrichment, dtype=float)
	ylb = np.asarray(borzoi_enrichment_lb, dtype=float)
	yub = np.asarray(borzoi_enrichment_ub, dtype=float)

	n = len(x)
	for arr, name in [(xlb, "meta_enrichment_lb"), (xub, "meta_enrichment_ub"),
					  (y, "borzoi_enrichment"), (ylb, "borzoi_enrichment_lb"), (yub, "borzoi_enrichment_ub")]:
		if len(arr) != n:
			raise ValueError(f"All inputs must have the same length; {name} has length {len(arr)} vs {n}")

	# Convert bounds to symmetric error lengths for matplotlib
	xerr = np.vstack([x - xlb, xub - x])
	yerr = np.vstack([y - ylb, yub - y])

	# Mask non-finite and negative error sizes
	finite = (
		np.isfinite(x) & np.isfinite(xlb) & np.isfinite(xub) &
		np.isfinite(y) & np.isfinite(ylb) & np.isfinite(yub)
	)
	if not np.all(finite):
		x, xlb, xub, y, ylb, yub, xerr, yerr = (
			x[finite], xlb[finite], xub[finite],
			y[finite], ylb[finite], yub[finite],
			xerr[:, finite], yerr[:, finite]
		)

	if np.any(xerr < 0) or np.any(yerr < 0):
		raise ValueError("Some CI bounds appear inverted (lb > estimate or ub < estimate). Check your inputs.")

	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

	ax.errorbar(
		x, y,
		xerr=xerr,
		yerr=yerr,
		fmt="o",
		capsize=capsize,
		elinewidth=1.0,
		markersize=5,
		alpha=point_alpha,
	)

	if annotate is not None:
		labels = np.asarray(annotate)
		if len(labels) != len(x):
			raise ValueError(f"annotate must have same length as data ({len(x)}), got {len(labels)}")
		for xi, yi, lab in zip(x, y, labels):
			ax.annotate(str(lab), (xi, yi), xytext=(4, 4), textcoords="offset points", fontsize=8)

	if show_identity:
		# Set limits first, then draw y=x across the visible range
		xmin = np.nanmin(x - xerr[0])
		xmax = np.nanmax(x + xerr[1])
		ymin = np.nanmin(y - yerr[0])
		ymax = np.nanmax(y + yerr[1])
		lo = np.nanmin([xmin, ymin])
		hi = np.nanmax([xmax, ymax])
		ax.set_xlim(lo, hi)
		ax.set_ylim(lo, hi)
		ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0)

	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if title is not None:
		ax.set_title(title)

	ax.grid(True, linewidth=0.5, alpha=0.4)
	fig.tight_layout()
	fig.savefig(outpath, bbox_inches="tight")
	plt.close(fig)

	return

def plot_borzoi_enrichment_ci(
	borzoi_enrichment,
	borzoi_enrichment_lb,
	borzoi_enrichment_ub,
	borzoi_anno_names,
	borzoi_enrichment_p,
	output_path,
	alpha=0.05,
	bonferroni=True,
	max_annos=None,
	figsize=(0.45, 5),
	dpi=300,
):
	"""
	Plot Borzoi enrichment with 95% CI.
	X-axis: annotation names
	Y-axis: enrichment (linear scale)
	Horizontal reference line at y = 1
	Filters to Bonferroni-significant annotations and orders by |enrichment|.
	"""

	enr   = np.asarray(borzoi_enrichment, dtype=float)
	lb    = np.asarray(borzoi_enrichment_lb, dtype=float)
	ub    = np.asarray(borzoi_enrichment_ub, dtype=float)
	pvals = np.asarray(borzoi_enrichment_p, dtype=float)
	names = np.asarray(borzoi_anno_names, dtype=object)

	if not (enr.shape == lb.shape == ub.shape == pvals.shape == names.shape):
		raise ValueError("All inputs must have the same shape")

	# Remove invalid entries
	valid = (
		np.isfinite(enr) &
		np.isfinite(lb) &
		np.isfinite(ub) &
		np.isfinite(pvals)
	)
	enr, lb, ub, pvals, names = enr[valid], lb[valid], ub[valid], pvals[valid], names[valid]

	# Bonferroni threshold
	m_tests = len(pvals)
	thr = alpha / m_tests if bonferroni else alpha

	keep = pvals < thr
	enr, lb, ub, pvals, names = enr[keep], lb[keep], ub[keep], pvals[keep], names[keep]

	if len(enr) == 0:
		raise ValueError(
			f"No annotations pass significance threshold (p < {thr:.3g})"
		)

	# Order by magnitude of enrichment
	order = np.argsort(-np.abs(enr))
	enr, lb, ub, pvals, names = enr[order], lb[order], ub[order], pvals[order], names[order]

	if max_annos is not None:
		enr, lb, ub, pvals, names = (
			enr[:max_annos], lb[:max_annos], ub[:max_annos],
			pvals[:max_annos], names[:max_annos]
		)

	n = len(enr)

	fig_w = max(6, figsize[0] * n)
	fig, ax = plt.subplots(figsize=(fig_w, figsize[1]))

	x = np.arange(n)
	yerr = np.vstack([enr - lb, ub - enr])

	ax.errorbar(
		x,
		enr,
		yerr=yerr,
		fmt="o",
		capsize=3,
		markersize=5,
		elinewidth=1,
	)

	# Reference line at enrichment = 1
	ax.axhline(1.0, linestyle="--", linewidth=1, color="black", alpha=0.7)

	ax.set_xticks(x)
	ax.set_xticklabels(names, rotation=45, ha="right")
	ax.set_ylabel("Enrichment")
	ax.set_xlabel("Annotation")

	#title = "Borzoi enrichment (Bonferroni-significant)" if bonferroni else "Borzoi enrichment (significant)"
	#title = "Enirhcment of borzoi significance within functional annotations"
	#ax.set_title(f"{title} (p < {thr:.3g}, n={n})")

	ax.grid(True, axis="y", alpha=0.4, linewidth=0.5)
	fig.tight_layout()
	fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
	plt.close(fig)

	return {
		"n_tests": m_tests,
		"threshold": thr,
		"n_kept": n,
		"kept_annotations": names.tolist(),
	}

def plot_borzoi_vs_gwas_enrichment_ci(
	borzoi_enrichment,
	borzoi_enrichment_lb,
	borzoi_enrichment_ub,
	borzoi_anno_names,
	borzoi_enrichment_p,
	gwas_enrichment,
	gwas_enrichment_lb,
	gwas_enrichment_ub,
	gwas_anno_names,
	gwas_enrichment_p,
	output_path,
	*,
	# filtering (off by default)
	filter_alpha=None,          # e.g. 0.05 to filter by p-value
	filter_on="either",         # {"either","both","borzoi","gwas"}
	# ordering
	order_by='abs_borzoi',              # None, "abs_borzoi", "abs_gwas", "abs_diff", "borzoi", "gwas"
	max_annos=None,
	figsize=(0.45, 5),
	dpi=300,
	check_name_alignment=False, # set True to assert names match exactly
	point_alpha=0.95,
	capsize=3,
	markersize=5,
	elinewidth=1,
):
	"""
	Plot Borzoi and GWAS enrichment with 95% CI.

	X-axis: annotation names
	Y-axis: enrichment (linear scale)
	Two series: Borzoi + GWAS (different colors), with slight x-offset so points don't overlap.
	Horizontal reference line at y = 1

	Assumes GWAS arrays are aligned to Borzoi by annotation order (same ordering).
	"""

	# --- coerce to arrays ---
	b_enr   = np.asarray(borzoi_enrichment, dtype=float)
	b_lb    = np.asarray(borzoi_enrichment_lb, dtype=float)
	b_ub    = np.asarray(borzoi_enrichment_ub, dtype=float)
	b_pvals = np.asarray(borzoi_enrichment_p, dtype=float)
	b_names = np.asarray(borzoi_anno_names, dtype=object)

	g_enr   = np.asarray(gwas_enrichment, dtype=float)
	g_lb    = np.asarray(gwas_enrichment_lb, dtype=float)
	g_ub    = np.asarray(gwas_enrichment_ub, dtype=float)
	g_pvals = np.asarray(gwas_enrichment_p, dtype=float)
	g_names = np.asarray(gwas_anno_names, dtype=object)

	if not (b_enr.shape == b_lb.shape == b_ub.shape == b_pvals.shape == b_names.shape):
		raise ValueError("All Borzoi inputs must have the same shape")
	if not (g_enr.shape == g_lb.shape == g_ub.shape == g_pvals.shape == g_names.shape):
		raise ValueError("All GWAS inputs must have the same shape")
	if b_enr.shape != g_enr.shape:
		raise ValueError("Borzoi and GWAS arrays must have the same shape (aligned by annotation order)")

	if check_name_alignment and not np.all(b_names == g_names):
		# show a helpful first mismatch
		mism = np.where(b_names != g_names)[0]
		i = int(mism[0])
		raise ValueError(f"Annotation name mismatch at index {i}: borzoi='{b_names[i]}' vs gwas='{g_names[i]}'")

	# --- remove invalid entries (drop any row with any non-finite across either series) ---
	valid = (
		np.isfinite(b_enr) & np.isfinite(b_lb) & np.isfinite(b_ub) & np.isfinite(b_pvals) &
		np.isfinite(g_enr) & np.isfinite(g_lb) & np.isfinite(g_ub) & np.isfinite(g_pvals)
	)
	b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[valid], b_lb[valid], b_ub[valid], b_pvals[valid], b_names[valid]
	g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[valid], g_lb[valid], g_ub[valid], g_pvals[valid], g_names[valid]

	# --- optional p-value filtering (OFF by default) ---
	if filter_alpha is not None:
		if filter_on == "either":
			keep = (b_pvals < filter_alpha) | (g_pvals < filter_alpha)
		elif filter_on == "both":
			keep = (b_pvals < filter_alpha) & (g_pvals < filter_alpha)
		elif filter_on == "borzoi":
			keep = (b_pvals < filter_alpha)
		elif filter_on == "gwas":
			keep = (g_pvals < filter_alpha)
		else:
			raise ValueError("filter_on must be one of {'either','both','borzoi','gwas'}")

		b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[keep], b_lb[keep], b_ub[keep], b_pvals[keep], b_names[keep]
		g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[keep], g_lb[keep], g_ub[keep], g_pvals[keep], g_names[keep]

	if len(b_enr) == 0:
		raise ValueError("No annotations remaining after filtering/validation")

	# --- optional ordering ---
	if order_by is not None:
		if order_by == "abs_borzoi":
			order = np.argsort(-np.abs(b_enr))
		elif order_by == "abs_gwas":
			order = np.argsort(-np.abs(g_enr))
		elif order_by == "abs_diff":
			order = np.argsort(-np.abs(b_enr - g_enr))
		elif order_by == "borzoi":
			order = np.argsort(-b_enr)
		elif order_by == "gwas":
			order = np.argsort(-g_enr)
		else:
			raise ValueError("order_by must be one of {None,'abs_borzoi','abs_gwas','abs_diff','borzoi','gwas'}")

		b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[order], b_lb[order], b_ub[order], b_pvals[order], b_names[order]
		g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[order], g_lb[order], g_ub[order], g_pvals[order], g_names[order]

	# --- cap number of annotations ---
	if max_annos is not None:
		b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[:max_annos], b_lb[:max_annos], b_ub[:max_annos], b_pvals[:max_annos], b_names[:max_annos]
		g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[:max_annos], g_lb[:max_annos], g_ub[:max_annos], g_pvals[:max_annos], g_names[:max_annos]

	n = len(b_enr)

	# --- plot ---
	fig_w = max(6, figsize[0] * n)
	fig, ax = plt.subplots(figsize=(fig_w, figsize[1]))

	x = np.arange(n)
	dx = 0.18  # horizontal dodge so two points don't overlap

	# Borzoi (blue)
	b_yerr = np.vstack([b_enr - b_lb, b_ub - b_enr])
	ax.errorbar(
		x - dx,
		b_enr,
		yerr=b_yerr,
		fmt="o",
		capsize=capsize,
		markersize=markersize,
		elinewidth=elinewidth,
		alpha=point_alpha,
		label="Borzoi",
		color="C0",
	)

	# GWAS (orange)
	g_yerr = np.vstack([g_enr - g_lb, g_ub - g_enr])
	ax.errorbar(
		x + dx,
		g_enr,
		yerr=g_yerr,
		fmt="o",
		capsize=capsize,
		markersize=markersize,
		elinewidth=elinewidth,
		alpha=point_alpha,
		label="GWAS",
		color="C1",
	)

	# Reference line at enrichment = 1
	ax.axhline(1.0, linestyle="--", linewidth=1, color="black", alpha=0.7)

	ax.set_xticks(x)
	ax.set_xticklabels(b_names, rotation=45, ha="right")
	ax.set_ylabel("Enrichment")
	ax.set_xlabel("Annotation")
	ax.legend(frameon=False)

	#ax.grid(True, axis="y", alpha=0.4, linewidth=0.5)
	fig.tight_layout()
	fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
	plt.close(fig)

	return {
		"n_tests": int(n),
		"n_kept": int(n),
		"kept_annotations": b_names.tolist(),
		"filter_alpha": filter_alpha,
		"filter_on": filter_on,
		"order_by": order_by,
	}



def plot_borzoi_vs_eqtl_enrichment_ci(
	borzoi_enrichment,
	borzoi_enrichment_lb,
	borzoi_enrichment_ub,
	borzoi_anno_names,
	borzoi_enrichment_p,
	gwas_enrichment,
	gwas_enrichment_lb,
	gwas_enrichment_ub,
	gwas_anno_names,
	gwas_enrichment_p,
	output_path,
	*,
	# filtering (off by default)
	filter_alpha=None,          # e.g. 0.05 to filter by p-value
	filter_on="either",         # {"either","both","borzoi","gwas"}
	# ordering
	order_by='abs_borzoi',              # None, "abs_borzoi", "abs_gwas", "abs_diff", "borzoi", "gwas"
	max_annos=None,
	figsize=(0.45, 5),
	dpi=300,
	check_name_alignment=False, # set True to assert names match exactly
	point_alpha=0.95,
	capsize=3,
	markersize=5,
	elinewidth=1,
):
	"""
	Plot Borzoi and GWAS enrichment with 95% CI.

	X-axis: annotation names
	Y-axis: enrichment (linear scale)
	Two series: Borzoi + GWAS (different colors), with slight x-offset so points don't overlap.
	Horizontal reference line at y = 1

	Assumes GWAS arrays are aligned to Borzoi by annotation order (same ordering).
	"""

	# --- coerce to arrays ---
	b_enr   = np.asarray(borzoi_enrichment, dtype=float)
	b_lb    = np.asarray(borzoi_enrichment_lb, dtype=float)
	b_ub    = np.asarray(borzoi_enrichment_ub, dtype=float)
	b_pvals = np.asarray(borzoi_enrichment_p, dtype=float)
	b_names = np.asarray(borzoi_anno_names, dtype=object)

	g_enr   = np.asarray(gwas_enrichment, dtype=float)
	g_lb    = np.asarray(gwas_enrichment_lb, dtype=float)
	g_ub    = np.asarray(gwas_enrichment_ub, dtype=float)
	g_pvals = np.asarray(gwas_enrichment_p, dtype=float)
	g_names = np.asarray(gwas_anno_names, dtype=object)

	if not (b_enr.shape == b_lb.shape == b_ub.shape == b_pvals.shape == b_names.shape):
		raise ValueError("All Borzoi inputs must have the same shape")
	if not (g_enr.shape == g_lb.shape == g_ub.shape == g_pvals.shape == g_names.shape):
		raise ValueError("All GWAS inputs must have the same shape")
	if b_enr.shape != g_enr.shape:
		raise ValueError("Borzoi and GWAS arrays must have the same shape (aligned by annotation order)")

	if check_name_alignment and not np.all(b_names == g_names):
		# show a helpful first mismatch
		mism = np.where(b_names != g_names)[0]
		i = int(mism[0])
		raise ValueError(f"Annotation name mismatch at index {i}: borzoi='{b_names[i]}' vs gwas='{g_names[i]}'")

	# --- remove invalid entries (drop any row with any non-finite across either series) ---
	valid = (
		np.isfinite(b_enr) & np.isfinite(b_lb) & np.isfinite(b_ub) & np.isfinite(b_pvals) &
		np.isfinite(g_enr) & np.isfinite(g_lb) & np.isfinite(g_ub) & np.isfinite(g_pvals)
	)
	b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[valid], b_lb[valid], b_ub[valid], b_pvals[valid], b_names[valid]
	g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[valid], g_lb[valid], g_ub[valid], g_pvals[valid], g_names[valid]

	# --- optional p-value filtering (OFF by default) ---
	if filter_alpha is not None:
		if filter_on == "either":
			keep = (b_pvals < filter_alpha) | (g_pvals < filter_alpha)
		elif filter_on == "both":
			keep = (b_pvals < filter_alpha) & (g_pvals < filter_alpha)
		elif filter_on == "borzoi":
			keep = (b_pvals < filter_alpha)
		elif filter_on == "gwas":
			keep = (g_pvals < filter_alpha)
		else:
			raise ValueError("filter_on must be one of {'either','both','borzoi','gwas'}")

		b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[keep], b_lb[keep], b_ub[keep], b_pvals[keep], b_names[keep]
		g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[keep], g_lb[keep], g_ub[keep], g_pvals[keep], g_names[keep]

	if len(b_enr) == 0:
		raise ValueError("No annotations remaining after filtering/validation")

	# --- optional ordering ---
	if order_by is not None:
		if order_by == "abs_borzoi":
			order = np.argsort(-np.abs(b_enr))
		elif order_by == "abs_gwas":
			order = np.argsort(-np.abs(g_enr))
		elif order_by == "abs_diff":
			order = np.argsort(-np.abs(b_enr - g_enr))
		elif order_by == "borzoi":
			order = np.argsort(-b_enr)
		elif order_by == "gwas":
			order = np.argsort(-g_enr)
		else:
			raise ValueError("order_by must be one of {None,'abs_borzoi','abs_gwas','abs_diff','borzoi','gwas'}")

		b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[order], b_lb[order], b_ub[order], b_pvals[order], b_names[order]
		g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[order], g_lb[order], g_ub[order], g_pvals[order], g_names[order]

	# --- cap number of annotations ---
	if max_annos is not None:
		b_enr, b_lb, b_ub, b_pvals, b_names = b_enr[:max_annos], b_lb[:max_annos], b_ub[:max_annos], b_pvals[:max_annos], b_names[:max_annos]
		g_enr, g_lb, g_ub, g_pvals, g_names = g_enr[:max_annos], g_lb[:max_annos], g_ub[:max_annos], g_pvals[:max_annos], g_names[:max_annos]

	n = len(b_enr)

	# --- plot ---
	fig_w = max(6, figsize[0] * n)
	fig, ax = plt.subplots(figsize=(fig_w, figsize[1]))

	x = np.arange(n)
	dx = 0.18  # horizontal dodge so two points don't overlap

	# Borzoi (blue)
	b_yerr = np.vstack([b_enr - b_lb, b_ub - b_enr])
	ax.errorbar(
		x - dx,
		b_enr,
		yerr=b_yerr,
		fmt="o",
		capsize=capsize,
		markersize=markersize,
		elinewidth=elinewidth,
		alpha=point_alpha,
		label="Borzoi",
		color="C0",
	)

	# GWAS (orange)
	g_yerr = np.vstack([g_enr - g_lb, g_ub - g_enr])
	ax.errorbar(
		x + dx,
		g_enr,
		yerr=g_yerr,
		fmt="o",
		capsize=capsize,
		markersize=markersize,
		elinewidth=elinewidth,
		alpha=point_alpha,
		label="eQTL",
		color="C1",
	)

	# Reference line at enrichment = 1
	ax.axhline(1.0, linestyle="--", linewidth=1, color="black", alpha=0.7)

	ax.set_xticks(x)
	ax.set_xticklabels(b_names, rotation=45, ha="right")
	ax.set_ylabel("Enrichment")
	ax.set_xlabel("Annotation")
	ax.legend(frameon=False)

	#ax.grid(True, axis="y", alpha=0.4, linewidth=0.5)
	fig.tight_layout()
	fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
	plt.close(fig)

	return {
		"n_tests": int(n),
		"n_kept": int(n),
		"kept_annotations": b_names.tolist(),
		"filter_alpha": filter_alpha,
		"filter_on": filter_on,
		"order_by": order_by,
	}


	
######################
# Command line args
#######################
output_dir = sys.argv[1]
variant_annotation_enrichment_dir = sys.argv[2]
gtex_tissue_names_file = sys.argv[3]
non_redundant_gwas_traits_file = sys.argv[4]
gwas_sldsc_results_dir = sys.argv[5]


'''
non_binary_annos = {}
non_binary_annos['GERP.NS'] = 1
non_binary_annos['MAF_Adj_Predicted_Allele_Age'] = 1
non_binary_annos['MAF_Adj_LLD_AFR'] = 1
non_binary_annos['Recomb_Rate_10kb'] = 1
non_binary_annos['Nucleotide_Diversity_10kb'] = 1
non_binary_annos['Backgrd_Selection_Stat'] = 1
non_binary_annos['CpG_Content_50kb'] = 1
non_binary_annos['MAF_Adj_ASMC'] = 1
non_binary_annos['GTEx_eQTL_MaxCPP'] = 1
non_binary_annos['BLUEPRINT_H3K27acQTL_MaxCPP'] = 1
non_binary_annos['BLUEPRINT_H3K4me1QTL_MaxCPP'] = 1
non_binary_annos['BLUEPRINT_DNA_methylation_MaxCPP'] = 1
non_binary_annos['Human_Enhancer_Villar_Species_Enhancer_Count'] = 1


# First create meta-analyzed gwas enrichment file
meta_gwas_enrichment_file = output_dir + 'meta_analyzed_gwas_enrichments_binary_anno_only.txt'
meta_enrichment, meta_enrichment_lb, meta_enrichment_ub, anno_names = make_meta_analyzed_gwas_enrichment_file(meta_gwas_enrichment_file, non_redundant_gwas_traits_file, gwas_sldsc_results_dir, non_binary_annos)


# Extract gtex tissue names
gtex_tissue_names = np.loadtxt(gtex_tissue_names_file,dtype=str,delimiter='\t')[1:]

tissue_name = 'Liver'
filter_to_binary_anno = True
enrichment_summary_file = variant_annotation_enrichment_dir + 'variant_enrichments_' + tissue_name + '_0.05_summary.txt'
# Extract eqtl enrichments
borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub, borzoi_anno_names, borzoi_enrichment_p = get_borzoi_anno_enrichments(enrichment_summary_file, non_binary_annos, filter_to_binary_anno)

plot_enrichment_scatter_with_ci(
	meta_enrichment, meta_enrichment_lb, meta_enrichment_ub,
	borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub,
	outpath=output_dir + "enrichment_scatter.pdf",
	title="Meta vs Borzoi enrichment (95% CI)"
)




tissue_name = 'Liver'
filter_to_binary_anno = False
enrichment_summary_file = variant_annotation_enrichment_dir + 'variant_enrichments_' + tissue_name + '_0.05_summary.txt'
# Extract eqtl enrichments
borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub, borzoi_anno_names, borzoi_enrichment_p = get_borzoi_anno_enrichments(enrichment_summary_file, non_binary_annos, filter_to_binary_anno)


res = plot_borzoi_enrichment_ci(
	borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub,
	borzoi_anno_names, borzoi_enrichment_p,
	output_path=output_dir + "borzoi_enrichment_ci.pdf",
)

'''


# Table S10 of luke o'connor polygenecity paper
# "Functional annotations"
# Though we selected 1 enhancer and 1 H3K27ac
# Further super enhancer (Vahedi) and typical enhancer (Vahedi) no longer existed
# Neither did we have low frequency UK 10K
# Also added 5' UTR

anno_subset = {}
anno_subset['Conserved_LindbladToh'] = 'Conserved'
anno_subset['TSS_Hoffman'] = 'TSS'
anno_subset['Coding_UCSC'] = 'Coding'
anno_subset['WeakEnhancer_Hoffman'] = 'Weak Enhancer'
anno_subset['SuperEnhancer_Hnisz'] = 'Super Enhancer'
anno_subset['PromoterFlanking_Hoffman'] = 'Promoter Flanking'
anno_subset['DGF_ENCODE'] = 'DGF'
anno_subset['Enhancer_Hoffman'] = 'Enhancer'
anno_subset['UTR_3_UCSC'] = 'UTR 3'
anno_subset['UTR_5_UCSC'] = 'UTR 5'
anno_subset['Promoter_UCSC'] = 'Promoter'
anno_subset['CTCF_Hoffman'] = 'CTCF'
anno_subset['TFBS_ENCODE'] = 'TFBS'
anno_subset['H3K9ac_Trynka'] = 'H3K9ac'
anno_subset['H3K4me3_Trynka'] = 'H3K4me3'
anno_subset['FetalDHS_Trynka'] = 'Fetal DHS'
anno_subset['DHS_Trynka'] = 'DHS'
anno_subset['H3K27ac_Hnisz'] = 'H3K27ac'
anno_subset['Transcr_Hoffman'] = 'Transcribed'
anno_subset['Intron_UCSC'] = 'Intron'


# GWAS
# First create meta-analyzed gwas enrichment file
meta_gwas_enrichment_file = output_dir + 'meta_analyzed_gwas_enrichments_binary_anno_only.txt'
meta_enrichment, meta_enrichment_lb, meta_enrichment_ub, anno_names = make_meta_analyzed_gwas_enrichment_file_for_specific_annos(meta_gwas_enrichment_file, non_redundant_gwas_traits_file, gwas_sldsc_results_dir, anno_subset)

# Borzoi
tissue_name = 'Whole_Blood'
enrichment_summary_file = variant_annotation_enrichment_dir + 'variant_enrichments_' + tissue_name + '_0.05_summary.txt'
# Extract eqtl enrichments
borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub, borzoi_anno_names, borzoi_enrichment_p = get_borzoi_anno_enrichments_for_specific_annos(enrichment_summary_file, anno_subset)


# eQTL
enrichment_summary_file = variant_annotation_enrichment_dir + 'fm_eqtl_variant_enrichments_' + tissue_name + '_0.9_summary.txt'
# Extract eqtl enrichments
eqtl_enrichment, eqtl_enrichment_lb, eqtl_enrichment_ub, eqtl_anno_names, eqtl_enrichment_p = get_borzoi_anno_enrichments_for_specific_annos(enrichment_summary_file, anno_subset)




plot_enrichment_scatter_with_ci_log(
	meta_enrichment, meta_enrichment_lb, meta_enrichment_ub,
	borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub,
	outpath=output_dir + "enrichment_scatter_functional_anno_only.pdf",
	title="Meta GWAS vs Borzoi enrichment (95% CI)"
)

res = plot_borzoi_enrichment_ci(
	borzoi_enrichment, borzoi_enrichment_lb, borzoi_enrichment_ub,
	borzoi_anno_names, borzoi_enrichment_p,bonferroni=False, alpha=1,
	output_path=output_dir + "borzoi_enrichment_functional_anno_only_ci.pdf",
)

plot_borzoi_vs_gwas_enrichment_ci

ordering = np.argsort(borzoi_enrichment)

res = plot_borzoi_vs_gwas_enrichment_ci(
	borzoi_enrichment[ordering], borzoi_enrichment_lb[ordering], borzoi_enrichment_ub[ordering],
	borzoi_anno_names[ordering], borzoi_enrichment_p[ordering], meta_enrichment[ordering], meta_enrichment_lb[ordering], meta_enrichment_ub[ordering], anno_names[ordering], borzoi_enrichment_p[ordering],
	output_path=output_dir + "borzoi_vs_gwas_enrichment_functional_anno_only_ci.pdf",
)


res = plot_borzoi_vs_eqtl_enrichment_ci(
	borzoi_enrichment[ordering], borzoi_enrichment_lb[ordering], borzoi_enrichment_ub[ordering],
	borzoi_anno_names[ordering], borzoi_enrichment_p[ordering], eqtl_enrichment[ordering], eqtl_enrichment_lb[ordering], eqtl_enrichment_ub[ordering], anno_names[ordering], eqtl_enrichment_p[ordering],
	output_path=output_dir + "borzoi_vs_eqtl_enrichment_functional_anno_only_ci.pdf",
)


