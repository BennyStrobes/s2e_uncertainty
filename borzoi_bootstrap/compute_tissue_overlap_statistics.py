import numpy as np
import os
import sys
import pdb
from itertools import combinations


def extract_tissue_names(gtex_tissue_names_file):
	f = open(gtex_tissue_names_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()

	return np.asarray(arr)


def extract_tissue_specific_sig_variant_gene_pairs(tissue_name, organized_borzoi_gtex_predictions_dir, sig_thresh, borzoi_effect_file_bins):
	dicti = {}
	for bin_name in borzoi_effect_file_bins:
		filer = organized_borzoi_gtex_predictions_dir + tissue_name + '_borzoi_estimates_w_uncertainty_' + bin_name + '.txt'
		f = open(filer)
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

			counter = counter +1

			if borzoi_pvalue > sig_thresh:
				continue

			variant_id = data[2]
			gene_id = data[3]
			var_gene_pair = variant_id + ':' + gene_id
			dicti[var_gene_pair] = 1
			
		f.close()


	return dicti



def jaccard_index_from_dicts(dict1, dict2):
	"""
	Compute Jaccard index between the keys of two dictionaries.

	J(A, B) = |A ∩ B| / |A ∪ B|
	"""
	k1 = set(dict1.keys())
	k2 = set(dict2.keys())

	union = k1 | k2
	if not union:
		return 0.0

	return len(k1 & k2) / len(union)


####################
# Command line args
sig_thresh = float(sys.argv[1])
tissue_overlap_output_file = sys.argv[2]
gtex_tissue_names_file = sys.argv[3]
organized_borzoi_gtex_predictions_dir = sys.argv[4] # Output file

borzoi_effect_file_bins = np.asarray(['0','1'])



# Extract ordered tissue names
tissue_names = extract_tissue_names(gtex_tissue_names_file)



# Extract significant variant gene pairs in each tissue
tissue_to_dictionary = {}
for tissue_name in tissue_names:
	print(tissue_name)
	tissue_specific_sig_dicti = extract_tissue_specific_sig_variant_gene_pairs(tissue_name, organized_borzoi_gtex_predictions_dir, sig_thresh, borzoi_effect_file_bins)
	tissue_to_dictionary[tissue_name] = tissue_specific_sig_dicti

t = open(tissue_overlap_output_file,'w')
t.write('tissue1\ttissue2\tjaccard_index\n')

# Now loop through all pairs of tissues
for t1, t2 in combinations(tissue_names, 2):
	dicti_t1 = tissue_to_dictionary[t1]
	dicti_t2 = tissue_to_dictionary[t2]

	overlap = jaccard_index_from_dicts(dicti_t1, dicti_t2)
	

	t.write(t1 + '\t' + t2 + '\t' + str(overlap) + '\n')

t.close()

