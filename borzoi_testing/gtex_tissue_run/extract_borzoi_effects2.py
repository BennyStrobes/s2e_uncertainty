import numpy as np
import os
import sys
import pdb
import h5py




fm_eqtls = sys.argv[1]
sed_h5 = sys.argv[2]
output = sys.argv[3]

# Variant gene pair to effect size mapping
mapping = {}
f = h5py.File(sed_h5, "r")

genes = np.asarray(f['gene']).astype(str)
snp_indices = np.asarray(f['si']).astype(int)
snp_ids = (np.asarray(f['snp'])[snp_indices]).astype(str)
log_seds = np.mean(np.asarray(f['logSED'][:,9:12]).astype(float),axis=1)

for ii, gene_id in enumerate(genes):
	snp_id = snp_ids[ii]
	snp_gene_pair = snp_id + ':' + gene_id
	mapping[snp_gene_pair] = log_seds[ii]
f.close()

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
	snp_id = data[2]
	gene_id = data[6]
	snp_gene_pair = snp_id + ':' + gene_id
	if snp_gene_pair in mapping:
		score = mapping[snp_gene_pair]
		t.write(line + '\t' + str(score) + '\n')
t.close()
f.close()