import numpy as np
import os
import sys
import pdb



def extract_ordered_target_names(full_gtex_target_file):
	f = open(full_gtex_target_file)
	sample_tissues = []
	tissues = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[-1]
		sample_name = data[1]

		sample_tissues.append(sample_name + ':' + tissue_name)
		tissues.append(tissue_name)
	f.close()
	sample_tissues = np.asarray(sample_tissues)
	tissues = np.asarray(tissues)

	# Quick error checking
	if len(sample_tissues) != len(np.unique(sample_tissues)):
		print('assumptin eorororo')
		pdb.set_trace()

	return sample_tissues, tissues

def extract_dictionary_list_of_fine_mapped_vg_pairs(processed_fm_eqtl_output_file):
	f = open(processed_fm_eqtl_output_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[2]
		gene_id = data[5]
		dicti[variant_id + ':' + gene_id] = 1
	f.close()

	return dicti






full_gtex_target_file = sys.argv[1]
processed_fm_eqtl_output_file = sys.argv[2]
organized_borzoi_gtex_predictions = sys.argv[3]


# Extract target_names
target_names, tissue_names = extract_ordered_target_names(full_gtex_target_file)
# Extract mapping from tissue name to target indices
ordered_tissues = np.sort(np.unique(tissue_names))


# Extract dictionary list of fine-mapped variant gene pairs
fm_variant_gene_pairs = extract_dictionary_list_of_fine_mapped_vg_pairs(processed_fm_eqtl_output_file)



vgt_to_bs_borzoi_info = {}
for tissue_name in ordered_tissues:
	print(tissue_name)
	filer = organized_borzoi_gtex_predictions + tissue_name + '_borzoi_estimates_w_uncertainty.txt'
	f = open(filer)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[2]
		gene_id = data[3]
		vgt = variant_id + ':' + gene_id + ':' + tissue_name
		vg = variant_id + ':' + gene_id
		if vg not in fm_variant_gene_pairs:
			continue
		meany = data[4]
		bs = data[5]
		if vgt in vgt_to_bs_borzoi_info:
			print('assumptioneororor')
			pdb.set_trace()
		vgt_to_bs_borzoi_info[vgt] = (meany, bs)
	f.close()

f = open(processed_fm_eqtl_output_file)
t = open(organized_borzoi_gtex_predictions + 'fm_organized_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects.txt','w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write('\t'.join(data) + '\t' + 'mean_borzoi_log_sed' + '\t' + 'sdev_borzoi_log_sed' + '\t' + 'bs_borzoi_log_sed' + '\n')
		continue
	variant = data[2]
	gene = data[5]
	tissue = 'Whole_Blood'
	vgt = variant + ':' + gene + ':' + tissue
	if vgt not in vgt_to_bs_borzoi_info:
		continue
	meany_string, bs_string = vgt_to_bs_borzoi_info[vgt]
	bs_vals = np.asarray(bs_string.split(';')).astype(float)
	sdevs = np.std(bs_vals)
	se_boot = np.std(bs_vals, ddof=1)


	t.write('\t'.join(data) + '\t' + meany_string + '\t' + str(se_boot) + '\t' + bs_string + '\n')

f.close()
t.close()


print(organized_borzoi_gtex_predictions + 'fm_organized_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects.txt')

