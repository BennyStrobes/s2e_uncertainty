import numpy as np
import os
import sys
import pdb
import pyarrow.parquet as pq
import pickle










#########################
# Command line args
#########################
raw_fine_mapping_eqtl_results_dir = sys.argv[1]
pip_threshold = float(sys.argv[2])
processed_fm_eqtl_output_file = sys.argv[3]

# Get dictionary of valid, autosomal chromosomes
valid_chroms= {}
for chrom_num in range(1,23):
	valid_chroms['chr' + str(chrom_num)] = 1


t = open(processed_fm_eqtl_output_file,'w')
t.write('tissue\tchromosome\tvariant_position\tvariant_hg38\tallele1\tallele2\tminor_allele\tgene\tpip\tbeta_posterior\tsd_posterior\tgene_type\n')

for file_name in os.listdir(raw_fine_mapping_eqtl_results_dir):
	if file_name.endswith('parquet') == False:
		continue

	full_file_name = raw_fine_mapping_eqtl_results_dir + file_name
	
	tissue_name =file_name.split('.v10')[0]
	print(tissue_name)

	pf = pq.ParquetFile(full_file_name)
	for rg in range(pf.num_row_groups):
		table = pf.read_row_group(rg)   # this is a chunk

		# Process fields and print to output
		array = np.asarray(table)
		# Loop through snp-gene pairs 
		nrows = array.shape[0]
		for row_iter in range(nrows):
			data = array[row_iter, :]
			gene_id = data[0]
			gene_type = data[2]
			variant_id = data[3]
			pip = float(data[4])
			if pip < pip_threshold:
				continue
			af = float(data[5])
			if np.isnan(af):
				continue

			afc = float(data[8])
			afc_se = float(data[9])

			variant_info = variant_id.split('_')

			variant_chrom = variant_info[0]

			if variant_chrom not in valid_chroms:
				continue
			
			variant_pos = variant_info[1]
			A1 = variant_info[2]
			A2 = variant_info[3]

			if len(A1) != 1:
				continue
			if len(A2) != 1:
				continue

			if af > .5:
				minor_allele = A1
			else:
				minor_allele = A2
			

			t.write(tissue_name + '\t' + variant_chrom + '\t' + variant_pos + '\t' + variant_id + '\t' + A1 + '\t' + A2 + '\t' + minor_allele + '\t' + gene_id + '\t' + str(pip) + '\t' + str(afc) + '\t' + str(afc_se) + '\t' + gene_type + '\n')

t.close()

'''
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
		t.write('tissue\t' + data[0] + '\t' + 'variant_position' + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\t' + data[7] + '\t' + data[11] + '\t' + data[12] + '\t' + data[16] + '\t' + data[18] + '\t' + data[19] + '\n')
		continue

	# Filter out irrelevent lines
	if data[8] != 'GTEx':
		continue
	if data[9] != 'SUSIE':
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

	if data[0] != data[4].split('_')[0]:
		continue
	t.write(data[10] + '\t' + data[0] + '\t' + hg38_variant_pos + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\t' + data[7] + '\t' + data[11] + '\t' + data[12] + '\t' + data[16] + '\t' + data[18] + '\t' + data[19] + '\n')

f.close()
t.close()
print(processed_fm_eqtl_output_file)

'''