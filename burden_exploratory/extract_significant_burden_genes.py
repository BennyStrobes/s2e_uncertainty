import numpy as np
import os
import sys
import pdb
import gzip











###################
# Command line args
####################
trait_names_file = sys.argv[1]
burden_test_data_dir = sys.argv[2]
sig_burden_genes_dir = sys.argv[3]
gene_annotation_summary_file = sys.argv[4]


# First create dictionary mapping from gene id to (ensamble_id, chrom, tss (hg19))
dicti = {}
f = open(gene_annotation_summary_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	gene_id = data[0]
	chrom_num = data[1]
	ensamble_id = data[2]
	tss = data[4]


	dicti[gene_id] = (ensamble_id, chrom_num, tss)
f.close()



# Extract burden test trait names
burden_trait_names = np.loadtxt(trait_names_file, dtype=str,delimiter='\t')[1:, 1]


burden_test_type = 'PTV_0.01'
pval_thresh = .05/20000

# Open output file
output_file = sig_burden_genes_dir + 'sig_burden_genes_' + burden_test_type + '.txt'
t = open(output_file,'w')
t.write('trait_name\tgene_name\tensamble_id\tchrom_num\tgene_tss\teffect_size\teffect_size_se\tpvalue\n')

# Loop through traits
used_gt_pairs = {}
used_genes = {}
for trait_name in burden_trait_names:

	trait_file = burden_test_data_dir + 'bolt_337K_unrelStringentBrit.' + trait_name + '.bgen.stats.gz'

	f = gzip.open(trait_file, 'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0].startswith(burden_test_type + '_') == False:
			continue
		pvalue = float(data[-1])
		if pvalue > pval_thresh:
			continue
		effect_size = float(data[10])
		effect_size_se = float(data[11])
		gene_name_info = data[0].split('_')
		if len(gene_name_info) != 3:
			print('assumpitnoenrornr')
			pdb.set_trace()
		gene_name = gene_name_info[2]

		if gene_name not in dicti:
			print('skip gene' + gene_name)
			continue

		gene_info = dicti[gene_name]

		gt_pair = gene_name + ';' + trait_name
		if gt_pair in used_gt_pairs:
			print('errorr')
			pdb.set_trace()
		if gene_name in used_genes:
			continue
		used_gt_pairs[gt_pair] = 1
		used_genes[gene_name] = 1
		t.write(trait_name + '\t' + gene_name + '\t' + gene_info[0] + '\t' + gene_info[1] + '\t' + gene_info[2] + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(pvalue) + '\n')

	f.close()
t.close()

print(output_file)