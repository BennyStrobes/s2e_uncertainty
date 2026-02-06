import numpy as np
import os
import pdb
import sys
import time

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



def create_mapping_from_tissue_name_to_target_indices(tissue_names, ordered_tissues):
	mapping = []
	for tissue_name in ordered_tissues:
		indices = np.where(tissue_names==tissue_name)[0]
		mapping.append(np.asarray(indices))

	return mapping

def create_mapping_from_ensamble_id_to_gene_tss_info(dist_to_tss_summary_file):
	dicti = {}
	f = open(dist_to_tss_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		chrom_num = data[0]
		pos = float(data[1])
		geneid = data[3].split('.')[0]
		strand = data[5]

		if geneid.startswith('ENSG') == False:
			print('assumptione roronro')
			pdb.set_trace()

		if strand != '-' and strand != '+':
			print('assumption oerororor')
			pdb.set_trace()


		dicti[geneid] = (chrom_num, pos, strand)
	f.close()

	return dicti

# Command line args
borzoi_gtex_predictions = sys.argv[1]
full_gtex_target_file = sys.argv[2]
organized_borzoi_gtex_predictions = sys.argv[3]


# Command line args
orig_bootstraps = np.arange(1,101)
bootstraps = []
for bs in orig_bootstraps:
	bootstraps.append(bs)
bootstraps = np.asarray(bootstraps)

parallel_batch_names = np.asarray(['0'])

# Extract target_names
target_names, tissue_names = extract_ordered_target_names(full_gtex_target_file)

# Extract mapping from tissue name to target indices
ordered_tissues = np.sort(np.unique(tissue_names))
tissue_name_to_target_indices = create_mapping_from_tissue_name_to_target_indices(tissue_names, ordered_tissues)

# Open output handles
t = []
for ii, tissue_name in enumerate(ordered_tissues):
	t.append(open(organized_borzoi_gtex_predictions + tissue_name + '_gene_span_borzoi_estimates_w_uncertainty.txt','w'))
# Print headers
for ii, tissue_name in enumerate(ordered_tissues):
	t[ii].write('chrom_num\tvariant_position\tvariant_name\tgene_name\tborzoi_mean_effect\tborzoi_bootstrapped_mean_effect\n')

# Now loop through parallelized batches
for parallel_batch_name in parallel_batch_names:

	# Open input file handles for this batch
	files = []
	for jj, bootstrap_iter in enumerate(bootstraps):
		files.append(open(borzoi_gtex_predictions + 'bs' + str(bootstrap_iter) + '_PIP_0.9_borzoi_pred_eqtl_effects_gene_span_borzoi_sed_results.txt'))

	head_counter = 0
	# Lines is across bootstraps
	for lines in zip(*files):
		if head_counter == 0:
			head_counter = head_counter + 1
			continue
		#timer1 = time.time()
		# Organize results across bootstraps into matrix
		borzoi_log_seds_arr = []
		for jj, line in enumerate(lines):
			line = line.rstrip()
			data = line.split('\t')
			if jj == 0:
				row_identifier = '\t'.join(data[:4])
			else:
				if row_identifier != '\t'.join(data[:4]):
					print('assumption erroorro')
					pdb.set_trace()
			
			borzoi_log_seds = np.asarray(data[4].split(';')).astype(float)
			borzoi_log_seds_arr.append(borzoi_log_seds)
		borzoi_log_seds_arr = np.asarray(borzoi_log_seds_arr)


		# Loop through tissues
		for ii, tissue_name in enumerate(ordered_tissues):
			tissue_log_seds = np.mean(borzoi_log_seds_arr[:, tissue_name_to_target_indices[ii]],axis=1)
			t[ii].write(row_identifier + '\t' + str(np.mean(tissue_log_seds)) + '\t' + ';'.join(tissue_log_seds.astype(str)) + '\n')
		#timer2 = time.time()
		#print(timer2-timer1)

	# Close input file handles for this batch
	for f in files:
		f.close()


# Close output file handles
for ii, tissue_name in enumerate(ordered_tissues):
	t[ii].close()
