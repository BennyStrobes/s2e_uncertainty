import numpy as np
import os
import sys
import pdb
import pyarrow.parquet as pq
import pandas as pd
import pyarrow.compute as pc
import h5py
import gzip


def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf):
	valid_chroms = {}
	for chrom_num in range(1,23):
		valid_chroms['chr' + str(chrom_num)] = 1


	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in valid_chroms:
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = 1

	f.close()

	return dicti








#######################
# Command line args
#######################
gtex_v10_pc_genes_gtf = sys.argv[1]
eqtl_sumstats_dir = sys.argv[2]
tissue_name = sys.argv[3]
eqtl_ss_output_dir = sys.argv[4]

# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)


# Open output file handle
output_file = eqtl_ss_output_dir + tissue_name + '_eqtl_sumstats.txt.gz'
t = gzip.open(output_file,'wt')
t.write('gene\tvariant\tchr\tsnp_pos\tA0\tA1\teqtl_effect_size\teqtl_effect_size_se\tN\tmaf\n')

for chrom_num in range(1,23):
	print(chrom_num)
	filer = eqtl_sumstats_dir + tissue_name + '.v10.allpairs.chr' + str(chrom_num) + '.parquet'
	pf = pq.ParquetFile(filer)

	for rg in range(pf.num_row_groups):
		table = pf.read_row_group(
			rg,
			columns=['gene_id', 'variant_id', 'tss_distance', 'af', 'slope', 'slope_se', 'ma_count']
		)

		if table.num_columns == 0:
			continue


		gene_col = table['gene_id']
		variant_col = table['variant_id']
		dist_col = table['tss_distance']
		af_col = table['af']
		slope_col = table['slope']
		slope_se_col = table['slope_se']


		gene_ids = table['gene_id'].to_pylist()
		variant_ids = table['variant_id'].to_pylist()
		tss_dists = table['tss_distance'].to_pylist()
		afs = table['af'].to_pylist()
		ma_counts = table['ma_count'].to_pylist()
		slopes = table['slope'].to_pylist()
		slope_ses = table['slope_se'].to_pylist()


		for ii, gene_id in enumerate(gene_ids):
			if gene_id.split('.')[0] not in pc_genes:
				continue

			maf = afs[ii]
			if maf > 0.5:
				maf = 1.0 - maf
			if maf == 0.0:
				continue
			pred_N = (ma_counts[ii]/maf)/2.0
			round_pred_N = np.round(pred_N)

			line_tss_dist = tss_dists[ii]
			if np.abs(line_tss_dist) > 100000:
				continue
			variant_id = variant_ids[ii]
			var_info = variant_id.split('_')
			if len(var_info[2]) != 1 or len(var_info[3]) != 1:
				continue

			af = afs[ii]
			if slopes[ii] is None or slope_ses[ii] is None:
				continue


			t.write(gene_id.split('.')[0] + '\t' + variant_id + '\t' + str(chrom_num) + '\t' + var_info[1] + '\t' + var_info[2] + '\t' + var_info[3] + '\t' + str(slopes[ii]) + '\t' + str(slope_ses[ii]) + '\t' + str(round_pred_N) + '\t' + str(maf) + '\n')
	t.flush()
t.close()