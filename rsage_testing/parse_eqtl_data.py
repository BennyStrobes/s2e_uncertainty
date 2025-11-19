import numpy as np
import os
import sys
import pdb










#########################
# Command line args
#########################
raw_fine_mapping_eqtl_results_file = sys.argv[1]
tissue_name = sys.argv[2]
pip_threshold = float(sys.argv[3])
processed_fm_eqtl_output_file = sys.argv[4]

# Get dictionary of valid, autosomal chromosomes
valid_chroms= {}
for chrom_num in range(1,23):
	valid_chroms['chr' + str(chrom_num)] = 1

f = open(raw_fine_mapping_eqtl_results_file)
t = open(processed_fm_eqtl_output_file,'w')
weird_cases = 0

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if len(data) != 20:
		weird_cases = weird_cases + 1
		continue
	if head_count == 0:
		head_count =head_count + 1
		t.write(data[0] + '\t' + 'variant_position' + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\t' + data[7] + '\t' + data[11] + '\t' + data[12] + '\t' + data[16] + '\t' + data[18] + '\t' + data[19] + '\n')
		continue

	# Filter out irrelevent lines
	if data[8] != 'GTEx':
		continue
	if data[9] != 'SUSIE':
		continue
	if data[10] != tissue_name:
		continue
	if float(data[16]) < pip_threshold:
		continue

	# Limit to snps only
	if len(data[5]) > 1 or len(data[6]) > 1:
		continue

	if data[4].split('_')[2] != data[5] or data[4].split('_')[3] != data[6]:
		if data[4].split('_')[2] != data[6] or data[4].split('_')[3] != data[5]:
			continue	
		continue

	if data[0] not in valid_chroms:
		continue


	hg38_variant_pos = data[4].split('_')[1]
	t.write(data[0] + '\t' + hg38_variant_pos + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\t' + data[7] + '\t' + data[11] + '\t' + data[12] + '\t' + data[16] + '\t' + data[18] + '\t' + data[19] + '\n')

f.close()
t.close()
print(processed_fm_eqtl_output_file)