import numpy as np
import os
import sys
import pdb
import gzip
from scipy.stats import fisher_exact
import statsmodels.api as sm

def create_mapping_from_borzoi_sig_to_variant_id(borzoi_variant_effect_file,sig_thresh):
	f = open(borzoi_variant_effect_file)
	dicti = {}
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

def extract_association_info(snp_anno_dir, variant_id_to_borzoi_sig, rsid_to_variantid, ct_num):
	sig_vector = []
	annotation_mat = []
	for chrom_num in range(1,23):		
		anno_file = snp_anno_dir + 'cell_type_group.' + str(ct_num) + '.' + str(chrom_num) + '.annot.gz'

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

			anno_vec = float(data[4])


			sig_vector.append(line_borzoi_sig)
			annotation_mat.append(anno_vec)

		f.close()


	return np.asarray(sig_vector), np.asarray(annotation_mat)


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
	or_val = ((a + correction) * (d + correction)) / (
		(b + correction) * (c + correction)
	)

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

def create_cell_type_mapping(filer):
	head_count = 0
	dicti = {}
	f = open(filer)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		dicti[data[0]] = data[1].split('.be')[0]
	f.close()
	return dicti



####################
# Command line args
####################
borzoi_variant_effect_file = sys.argv[1]
genotype_dir = sys.argv[2]
snp_anno_dir = sys.argv[3]
sig_thresh = float(sys.argv[4])
variant_annotation_enrichment_file = sys.argv[5]

t = open(variant_annotation_enrichment_file,'w')
t.write('annotation_name\todds_ratio\tpvalue\taa\tbb\tcc\tdd\n')


# Create mapping from cell-type number to cell type name
ct_number_to_ct_name = create_cell_type_mapping(snp_anno_dir + 'names')

# First extract mapping from variants to borzoi significnace (in at least one gene)
variant_id_to_borzoi_sig = create_mapping_from_borzoi_sig_to_variant_id(borzoi_variant_effect_file, sig_thresh)

# Now create mapping from rsid to variant id
rsid_to_variantid = create_mapping_from_rsid_to_variant_id(genotype_dir)


# Extract annotations as matrices and sig as vector
annotation_mat = []
annotation_names = []
sig_mat = []
for ct_num in range(1,11):
	annotation_name = ct_number_to_ct_name[str(ct_num)]
	sig_vector, annotation_vec = extract_association_info(snp_anno_dir, variant_id_to_borzoi_sig, rsid_to_variantid, ct_num)

	sig_mat.append(sig_vector)
	annotation_mat.append(annotation_vec)
	annotation_names.append(annotation_name)
annotation_mat = np.transpose(np.asarray(annotation_mat))
annotation_names = np.asarray(annotation_names)
sig_mat = np.transpose(np.asarray(sig_mat))

# Error checking
row_constant = np.all(sig_mat == sig_mat[:, [0]], axis=1)
if np.sum(row_constant == False) > 0:
	print('assumption eroror')

sig_vector = sig_mat[:,0]



t = open(variant_annotation_enrichment_file,'w')
t.write('annotation_name\todds_ratio\tpvalue\taa\tbb\tcc\tdd\n')
for annotation_iter, annotation_name in enumerate(annotation_names):
	if len(np.unique(annotation_mat[:,annotation_iter])) != 2:
		annotation_mat[:,annotation_iter] = 1.0*(annotation_mat[:,annotation_iter] > np.mean(annotation_mat[:,annotation_iter]))
		print('assumptioneornro')
		pdb.set_trace()

	or_obj = odds_ratio_and_pvalue(sig_vector, annotation_mat[:,annotation_iter])
	t.write(annotation_name + '\t' + str(or_obj['odds_ratio']) + '\t' + str(or_obj['p_value']) + '\t' + str(or_obj['a']) + '\t' + str(or_obj['b']) + '\t' + str(or_obj['c']) + '\t' + str(or_obj['d']) + '\n')

t.close()


'''
X = annotation_mat          # shape (N, K)
y = sig_vector.astype(int) # must be 0/1

# Add intercept explicitly
X_sm = sm.add_constant(X)

pdb.set_trace()
'''
