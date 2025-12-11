import numpy as np
import os
import sys
import pdb
import h5py
import re
import scipy.stats


def to_underscore(s: str) -> str:
	s = s.replace(")", "")
	# replace " (" OR any run of spaces and hyphens with "_"
	return re.sub(r"(?:\s*\(|[ -]+)", "_", s)


fm_eqtls = sys.argv[1]
sed_h5 = sys.argv[2]
output = sys.argv[3]
gtex_sample_attribute_file = sys.argv[4]
targets_file = sys.argv[5]


#####################
# First create mapping from gtex sample name to tissue name
#####################
gtex_sample_to_tissue_name = {}
gtex_tissues = {}
f = open(gtex_sample_attribute_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gtex_sample_id = data[0]
	tissue_type = to_underscore(data[6])
	if tissue_type == 'Cells_EBV_transformed_lymphocytes':
		tissue_type = 'Cells_EBV-transformed_lymphocytes'
	elif tissue_type == 'Brain_Spinal_cord_cervical_c_1':
		tissue_type = 'Brain_Spinal_cord_cervical_c-1'
	gtex_sample_to_tissue_name[gtex_sample_id] = tissue_type
	gtex_tissues[tissue_type] = 1
f.close()



#####################
# Second create mapping from gtex tissue name to targets
#####################
gtex_tissue_name_to_target_indices = {}
f = open(targets_file)
head_count = 0
target_counter = 0
target_tissue_names_arr = []
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	short_gtex_sample_name = data[1].split('.')[0]
	if short_gtex_sample_name not in gtex_sample_to_tissue_name:
		print('assumption eroror')
		pdb.set_trace()
	target_tissue_names_arr.append(gtex_sample_to_tissue_name[short_gtex_sample_name])
f.close()
target_tissue_names_arr = np.asarray(target_tissue_names_arr)
unique_target_tissues = np.unique(target_tissue_names_arr)

for target_tissue in unique_target_tissues:
	indices = np.where(target_tissue_names_arr == target_tissue)[0]
	#indices = np.where(target_tissue_names_arr == np.random.choice(unique_target_tissues))[0]
	gtex_tissue_name_to_target_indices[target_tissue] = indices


#####################
# Now fill in variant effect sizes
#####################


# Variant gene pair to effect size mapping
mapping = {}
ff = h5py.File(sed_h5, "r")


genes = np.asarray(ff['gene']).astype(str)
snp_indices = np.asarray(ff['si']).astype(int)
snp_ids = (np.asarray(ff['snp'])[snp_indices]).astype(str)
log_seds = np.asarray(ff['logSED']).astype(float)

for ii, gene_id in enumerate(genes):
	snp_id = snp_ids[ii]
	snp_gene_pair = snp_id + ':' + gene_id
	mapping[snp_gene_pair] = log_seds[ii]

f = open(fm_eqtls)
head_count = 0
t = open(output,'w')
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write('\t'.join(np.asarray(data)) + '\t' + 'borzoi_log_sed\n')
		continue
	tissue_type = data[0]
	snp_id = data[3]
	gene_id = data[7]
	snp_gene_pair = snp_id + ':' + gene_id

	if snp_gene_pair in mapping and tissue_type in gtex_tissue_name_to_target_indices:
		all_scores = mapping[snp_gene_pair]
		target_indices = gtex_tissue_name_to_target_indices[tissue_type]

		t.write(line + '\t' + str(np.mean(all_scores[target_indices])) + '\n')

t.close()
f.close()
ff.close()



aa = np.loadtxt(output,dtype=str,delimiter='\t')

susie = aa[1:,-3].astype(float)
borzoi = aa[1:,-1].astype(float)


print(scipy.stats.pearsonr(susie,borzoi))
print(scipy.stats.pearsonr(susie[np.abs(borzoi) > .1],borzoi[np.abs(borzoi) > .1]))
print(scipy.stats.pearsonr(susie[np.abs(borzoi) > .2],borzoi[np.abs(borzoi) > .2]))
print(scipy.stats.pearsonr(susie[np.abs(borzoi) > .3],borzoi[np.abs(borzoi) > .3]))
print(scipy.stats.pearsonr(susie[np.abs(borzoi) > .5],borzoi[np.abs(borzoi) > .5]))

