import numpy as np
import os
import sys
import pdb
import pyarrow.parquet as pq
import pandas as pd





def extract_ordered_list_of_all_variant_ids_on_this_chromosome(gtex_sumstats_dir, chrom_string):
	variant_ids = {}
	counter = 0
	for file_name in os.listdir(gtex_sumstats_dir):
		if file_name.endswith('chr' + str(chrom_num) + '.parquet') == False:
			continue

		counter = counter + 1
		print(file_name)
		pf = pq.ParquetFile(gtex_sumstats_dir + file_name)
		for rg in range(pf.num_row_groups):
			table = pf.read_row_group(rg)   # this is a chunk

			# Process fields and print to output
			array = np.asarray(table)

			# Loop through snp-gene pairs 
			nrows = array.shape[0]
			for row_iter in range(nrows):
				data = array[row_iter, :]
				variant_id = data[1]
				af = float(data[3])
				if af > .5:
					af = 1.0 - af
				if af < .05:
					continue
				variant_ids[variant_id] = 1

	arr_tupler = []
	for variant_id in [*variant_ids]:
		info = variant_id.split('_')
		arr_tupler.append((info[0], int(info[1]), variant_id, info[2], info[3]))

	sorted_list = sorted(arr_tupler, key=lambda x: x[1])

	return sorted_list





# Command line args
gtex_sumstats_dir = sys.argv[1]
output_vcf_file = sys.argv[2]


t = open(output_vcf_file,'w')

for chrom_num in range(1,23):
	print(chrom_num)

	# Extract ordered list of variant ids in this chromosome
	variant_tuples = extract_ordered_list_of_all_variant_ids_on_this_chromosome(gtex_sumstats_dir, str(chrom_num))

	for variant_tuple in variant_tuples:
		# Limit to SNVs
		if len(variant_tuple[3]) != 1 or len(variant_tuple[4]) != 1:
			continue
		t.write(variant_tuple[0] + '\t' + str(variant_tuple[1]) + '\t' + variant_tuple[2] + '\t' + variant_tuple[3] + '\t' + variant_tuple[4] + '\n')


t.close()

