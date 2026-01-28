import numpy as np
import os
import sys
import pdb
import gzip
from scipy.stats import fisher_exact, norm
import pyarrow.parquet as pq
import pandas as pd




def create_mapping_from_variant_id_to_eqtl_sig(fm_eqtl_results_file, eqtl_sumstats_dir, tissue_name):
	# load in fm eqtl results
	fm_variant_gene_pairs = {}
	fm_variants = {}
	f = open(fm_eqtl_results_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		if data[0] != tissue_name:
			continue
		variant_id = data[3]
		gene_id = data[7]
		fm_variants[variant_id] = 0
		fm_variant_gene_pairs[variant_id + ':' + gene_id] = 0
	f.close()


	# Load in non-sig eqtls
	dicti = {}
	for chrom_num in range(1,23):
		print(chrom_num)
		filer = eqtl_sumstats_dir + tissue_name + '.v10.allpairs.chr' + str(chrom_num) + '.parquet'
		pf = pq.ParquetFile(filer)
		for rg in range(pf.num_row_groups):
			table = pf.read_row_group(rg)   # this is a chunk

			# Process fields and print to output
			array = np.asarray(table)

			# Loop through snp-gene pairs 
			nrows = array.shape[0]
			for row_iter in range(nrows):
				data = array[row_iter, :]
				variant_id = data[1]
				geneid = data[0]
				af = float(data[3])
				if af > .5:
					af = 1.0 - af
				if af < .05:
					continue

				sig_val = 0.0
				if variant_id in fm_variants:
					sig_val = 1.0
				variant_info = variant_id.split('_')
				alt_variant_id = variant_info[0] + '_' + variant_info[1] + '_' + variant_info[3] + '_' + variant_info[2] + '_' + variant_info[4]
				if alt_variant_id in fm_variants:
					sig_val = 1.0
					pdb.set_trace()
				if variant_id not in dicti:
					dicti[variant_id] = sig_val
				else:
					if sig_val == 1.0:
						dicti[variant_id] = sig_val

	return dicti





def create_mapping_from_borzoi_sig_to_variant_id(borzoi_variant_effect_file_stem, borzoi_effect_file_bins,sig_thresh):
	dicti = {}

	for file_bin in borzoi_effect_file_bins:

		borzoi_variant_effect_file = borzoi_variant_effect_file_stem + file_bin + '.txt'

		f = open(borzoi_variant_effect_file)
		head_count = 0
		counter = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue


			bs_vals = np.asarray(data[5].split(';')).astype(float)

			borzoi_effect = np.mean(bs_vals)
			np_p_plus = (1 + np.sum(bs_vals >= 0))/(1+len(bs_vals))
			np_p_minus = (1 + np.sum(bs_vals <=0))/(1+len(bs_vals))
			borzoi_pvalue = 2.0*np.min((np_p_plus, np_p_minus))

			if borzoi_pvalue <= sig_thresh:
				booler = 1.0
			else:
				booler = 0.0


			variant_id = data[2]

			if variant_id not in dicti:
				dicti[variant_id] = booler
			else:
				if booler == 1.0:
					dicti[variant_id] = booler
			

			counter = counter +1

		f.close()

	return dicti



def create_mapping_from_rsid_to_variant_id(genotype_dir):
	dicti = {}
	for chrom_num in range(1,23):
		chrom_bim_file = genotype_dir + str(chrom_num) + '.bim'
		f = open(chrom_bim_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			rsid = data[1]
			variant_id1 = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5] + '_b38'
			variant_id2 = 'chr' + data[0] + '_' + data[3] + '_' + data[5] + '_' + data[4] + '_b38'

			if rsid in dicti:
				print('assumption eororr')
				pdb.set_trace()

			dicti[rsid] = (variant_id1, variant_id2)


		f.close()
	return dicti

def extract_association_info(snp_anno_dir, variant_id_to_borzoi_sig, rsid_to_variantid, remove_coding_boolean):
	sig_vector = []
	annotation_mat = []
	for chrom_num in range(1,23):
		print(chrom_num)

		anno_file = snp_anno_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'

		head_count = 0
		f = gzip.open(anno_file, 'rt')
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				if chrom_num == 1:
					annotation_names = np.asarray(data[5:])
				continue
			rsid = data[2]

			if rsid not in rsid_to_variantid:
				continue
			variant_ids = rsid_to_variantid[rsid]

			if variant_ids[0] not in variant_id_to_borzoi_sig and variant_ids[1] not in variant_id_to_borzoi_sig:
				continue

			# Error check
			if variant_ids[0] in variant_id_to_borzoi_sig and variant_ids[1] in variant_id_to_borzoi_sig:
				print('assumption eororro')
				pdb.set_trace()

			if variant_ids[0] in variant_id_to_borzoi_sig:
				line_borzoi_sig = variant_id_to_borzoi_sig[variant_ids[0]]
			elif variant_ids[1] in variant_id_to_borzoi_sig:
				line_borzoi_sig = variant_id_to_borzoi_sig[variant_ids[1]]
			else:
				print('assumption erooror')
				pdb.set_trace()

			anno_vec = np.asarray(data[5:]).astype(float)

			# Skip if coding
			if remove_coding_boolean:
				if anno_vec[0] == 1.0 or anno_vec[75] == 1.0 or anno_vec[76] == 1.0:
					continue

			sig_vector.append(line_borzoi_sig)
			annotation_mat.append(anno_vec)

		f.close()


	return np.asarray(sig_vector), np.asarray(annotation_mat), np.asarray(annotation_names)


def odds_ratio_and_pvalue(aa, bb, correction=0.1):
	"""
	Compute contingency table, odds ratio, and p-value for two binary vectors.

	Parameters
	----------
	aa, bb : array-like of {0,1}
		Binary exposure and outcome vectors (same length)
	correction : float
		Haldane–Anscombe correction added to each cell (default 0.5)

	Returns
	-------
	results : dict
		a, b, c, d, odds_ratio, p_value
	"""
	aa = np.asarray(aa)
	bb = np.asarray(bb)

	if aa.shape != bb.shape:
		raise ValueError("aa and bb must have the same shape")

	a = np.sum((aa == 1) & (bb == 1))
	b = np.sum((aa == 1) & (bb == 0))
	c = np.sum((aa == 0) & (bb == 1))
	d = np.sum((aa == 0) & (bb == 0))

	# Odds ratio with correction
	or_val = ((a + correction) * (d + correction)) / ((b + correction) * (c + correction))

	# Fisher exact test (two-sided)
	table = [[a, b],
			 [c, d]]
	_, p_value = fisher_exact(table, alternative="two-sided")

	return {
		"a": a,
		"b": b,
		"c": c,
		"d": d,
		"odds_ratio": or_val,
		"p_value": p_value
	}



def enrichment_with_se(aa, bb, correction=0.0, alpha=0.05):
	"""
	Fold-enrichment (relative risk) with delta-method SE on log scale
	and CI, plus Fisher exact p-value.

	If correction>0, it is added to each cell before computing RR and SE
	(useful when any cell is 0).
	"""
	aa = np.asarray(aa)
	bb = np.asarray(bb)

	if aa.shape != bb.shape:
		raise ValueError("aa and bb must have the same shape")

	a = np.sum((aa == 1) & (bb == 1))
	b = np.sum((aa == 1) & (bb == 0))
	c = np.sum((aa == 0) & (bb == 1))
	d = np.sum((aa == 0) & (bb == 0))

	# Fisher exact p-value uses raw counts
	_, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")

	# Add pseudocounts for RR/SE if requested
	a_, b_, c_, d_ = a + correction, b + correction, c + correction, d + correction

	# Risks
	p1 = a_ / (a_ + b_)
	p0 = c_ / (c_ + d_)

	enrichment = p1 / p0
	log_enrichment = np.log(enrichment)

	# Delta-method variance of log(RR)
	# (works when a_, c_ > 0 and denominators > 0; correction helps ensure this)
	var_log = (1.0 / a_) - (1.0 / (a_ + b_)) + (1.0 / c_) - (1.0 / (c_ + d_))
	se_log = np.sqrt(var_log)

	z = norm.ppf(1 - alpha / 2)
	ci_low = np.exp(log_enrichment - z * se_log)
	ci_high = np.exp(log_enrichment + z * se_log)

	return {
		"a": int(a), "b": int(b), "c": int(c), "d": int(d),
		"rate_exposed": float(p1),
		"rate_unexposed": float(p0),
		"enrichment": float(enrichment),
		"log_enrichment": float(log_enrichment),
		"se_log_enrichment": float(se_log),
		"ci": (float(ci_low), float(ci_high)),
		"p_value": float(p_value),
	}


####################
# Command line args
####################
fm_eqtl_results_file = sys.argv[1]
eqtl_sumstats_dir= sys.argv[2]
genotype_dir = sys.argv[3]
snp_anno_dir = sys.argv[4]
sig_thresh = float(sys.argv[5])
variant_annotation_enrichment_file = sys.argv[6]
remove_coding_boolean_string = sys.argv[7]
tissue_name = sys.argv[8]
if remove_coding_boolean_string == 'True':
	remove_coding_boolean = True
elif remove_coding_boolean_string == 'False':
	remove_coding_boolean = False
else:
	print('assumption erorno')
	pdb.set_trace()






t = open(variant_annotation_enrichment_file,'w')
t.write('annotation_name\tenrichment\tpvalue\tenrichment_95_lb\tenrichment_95_ub\taa\tbb\tcc\tdd\n')


# First extract mapping from variants to borzoi significnace (in at least one gene)
variant_id_to_eqtl_sig = create_mapping_from_variant_id_to_eqtl_sig(fm_eqtl_results_file, eqtl_sumstats_dir, tissue_name)

# Now create mapping from rsid to variant id
rsid_to_variantid = create_mapping_from_rsid_to_variant_id(genotype_dir)


# Extract annotations as matrices and sig as vector
sig_vector, annotation_mat, annotation_names = extract_association_info(snp_anno_dir, variant_id_to_eqtl_sig, rsid_to_variantid, remove_coding_boolean)


for annotation_iter, annotation_name in enumerate(annotation_names):
	if len(np.unique(annotation_mat[:,annotation_iter])) != 2:
		annotation_mat[:,annotation_iter] = 1.0*(annotation_mat[:,annotation_iter] > np.mean(annotation_mat[:,annotation_iter]))


	enrichment_obj = enrichment_with_se(sig_vector, annotation_mat[:,annotation_iter])

	t.write(annotation_name + '\t' + str(enrichment_obj['enrichment']) + '\t' + str(enrichment_obj['p_value']) + '\t' + str(enrichment_obj['ci'][0]) + '\t' + str(enrichment_obj['ci'][1]) + '\t' + str(enrichment_obj['a']) + '\t' + str(enrichment_obj['b']) + '\t' + str(enrichment_obj['c']) + '\t' + str(enrichment_obj['d']) + '\n')

t.close()