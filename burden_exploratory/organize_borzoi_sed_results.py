import numpy as np
import os
import sys
import pdb
import h5py
import re



def initialize_vg_pairs_to_pred_dictionary(vg_pairs_to_test_file):
	f = open(vg_pairs_to_test_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		variant_id = data[1]
		ens_id = data[8].split('.')[0]

		vg_pair = variant_id + ':' + ens_id

		if vg_pair in dicti:
			print('assumption eornornro')
			pdb.set_trace()
		dicti[vg_pair] = []
	f.close()

	return dicti





######################
# Commmand line args
######################
borzoi_pred_input_stem = sys.argv[1]
output_file = sys.argv[2]
vg_pairs_to_test_file = sys.argv[3]

# Initially mapping from vg pairs to array
vg_pairs_to_preds = initialize_vg_pairs_to_pred_dictionary(vg_pairs_to_test_file)

# Loop through model reps
reps = np.arange(0,4)
for rep in reps:

	sed_h5 = borzoi_pred_input_stem + str(rep) + '/sed.h5'
	ff = h5py.File(sed_h5, "r")


	genes = np.asarray(ff['gene']).astype(str)
	snp_indices = np.asarray(ff['si']).astype(int)
	snp_ids = (np.asarray(ff['snp'])[snp_indices]).astype(str)
	log_seds = np.asarray(ff['logSED']).astype(float)

	for ii, gene_id in enumerate(genes):
		snp_id = snp_ids[ii]
		snp_gene_pair = snp_id + ':' + gene_id.split('.')[0]

		if snp_gene_pair in vg_pairs_to_preds:
			vg_pairs_to_preds[snp_gene_pair].append(log_seds[ii])
	ff.close()


f = open(vg_pairs_to_test_file)
t = open(output_file,'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'avg_log_sed\n')
		continue

	variant_id = data[1]
	ens_id = data[8].split('.')[0]

	vg_pair = variant_id + ':' + ens_id

	effect = np.mean(np.asarray(vg_pairs_to_preds[vg_pair]))

	t.write(line + '\t' + str(effect) + '\n')
	if np.abs(effect) > .05:
		print(line + '\t' + str(effect))


f.close()
t.close()












