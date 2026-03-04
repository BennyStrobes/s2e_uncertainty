import numpy as np
import os
import sys
import pdb







####################
# Command line args
###################

variant_gene_pairs_dir = sys.argv[1]

output_file = variant_gene_pairs_dir + 'burden_variant_gene_pairs_dist_25000.txt'
t = open(output_file,'w')

for chrom_num in range(1,23):

	input_file = variant_gene_pairs_dir + 'burden_variant_gene_pairs_dist_25000_' + str(chrom_num) + '.txt'

	head_count = 0
	f = open(input_file)
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			if chrom_num == 1:
				t.write(line + '\n')
			continue
		t.write(line + '\n')
	f.close()

t.close()