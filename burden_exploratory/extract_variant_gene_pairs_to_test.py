import numpy as np
import os
import sys
import pdb
from bgen.reader import BgenFile







def extract_burden_genes_on_this_chrom(chrom_string, burden_genes_summary_file):
	arr = []
	used = {}
	f = open(burden_genes_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[1]
		line_chrom_string = data[3]
		if line_chrom_string != chrom_string:
			continue

		if gene_name in used:
			continue
		used[gene_name] = 1

		arr.append((gene_name, data[4], data[2]))
	f.close()
	return arr





######################
# command line args
######################
burden_genes_summary_file = sys.argv[1]
genotype_dir = sys.argv[2]
variant_gene_pairs_dir = sys.argv[3]

distance_window = 25000


t = open(variant_gene_pairs_dir + 'burden_variant_gene_pairs_dist_' + str(distance_window) + '.txt','w')
t.write('chr\tvariant_id\tvariant_position\tA1\tA2\talt_AF\tgene_name\tgene_tss\tensamble_id\n')

for chrom_num in range(1,23):
	print(chrom_num)

	# Extract burden genes and tss for this chromosome
	chrom_burden_genes = extract_burden_genes_on_this_chrom('chr' + str(chrom_num), burden_genes_summary_file)

	# Skip chrom with no burden genes
	if len(chrom_burden_genes) == 0:
		continue

	chrom_arr = ['NULL']*250000000


	# Loop through burden genes on this chromosome
	gene_to_tss = {}
	gene_to_ens = {}
	for burden_gene_info in chrom_burden_genes:
		# Extract relevent info
		burden_gene_name = burden_gene_info[0]
		burden_gene_tss = int(burden_gene_info[1])
		burden_gene_ens = burden_gene_info[2]
		gene_to_tss[burden_gene_name] = burden_gene_tss
		gene_to_ens[burden_gene_name] = burden_gene_ens

		start = burden_gene_tss - distance_window
		end = burden_gene_tss + distance_window

		for chrom_pos in np.arange(start, end+1):
			if chrom_arr[chrom_pos] == 'NULL':
				chrom_arr[chrom_pos] = burden_gene_name
			else:
				chrom_arr[chrom_pos] = chrom_arr[chrom_pos] + ';' + burden_gene_name

	# Open chrom bgen file
	chrom_bgen_file = genotype_dir + 'UKB_MAF0.001_v3.' + str(chrom_num) + '.bgen'

	bfile = BgenFile(chrom_bgen_file, delay_parsing=True)

	for var in bfile:
		# var.position is 1-based genomic position
		variant_position = var.pos

		# Quick error check
		if var.varid.startswith(str(chrom_num) + ':'):
			tmp_a1 = var.varid.split('_')[1]
			tmp_a2 = var.varid.split('_')[2]
			if var.alleles[0] != tmp_a1:
				print('skippy2')
				continue
			if var.alleles[1] != tmp_a2:
				print('skippy')
				continue

		if chrom_arr[variant_position] == 'NULL':
			continue

		af = np.sum(var.alt_dosage)/(2.0*len(var.alt_dosage))
		# Filter out common variants
		if af > 0.01 and af < 0.99:
			continue

		if len(var.alleles[0]) != 1 or len(var.alleles[1]) != 1:
			continue

		var_genes = np.asarray(chrom_arr[variant_position].split(';'))
	
		new_var_id = 'chr' + str(chrom_num) + '_' + str(variant_position) + '_' + var.alleles[0] + '_' + var.alleles[1] + '_b37'

		for var_gene in var_genes:
			t.write('chr' + str(chrom_num) + '\t' + new_var_id + '\t' + str(variant_position) + '\t' + var.alleles[0] + '\t' + var.alleles[1] + '\t' + str(af) + '\t' + var_gene + '\t' + str(gene_to_tss[var_gene]) + '\t' + gene_to_ens[var_gene] + '\n')

	bfile.close()
	t.flush()
t.close()
