import numpy as np
import os
import sys
import pdb
import pyarrow.parquet as pq
import pandas as pd

def create_mapping_from_ensamble_id_to_gene_tss_info(dist_to_tss_summary_file):
	dicti = {}
	f = open(dist_to_tss_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		chrom_num = data[0]
		pos = float(data[1])
		geneid = data[3].split('.')[0]
		strand = data[5]

		if geneid.startswith('ENSG') == False:
			print('assumptione roronro')
			pdb.set_trace()

		if strand != '-' and strand != '+':
			print('assumption oerororor')
			pdb.set_trace()


		dicti[geneid] = (chrom_num, pos, strand)
	f.close()

	return dicti



def create_mapping_from_variant_gene_pair_to_borzoi_sig(borzoi_variant_effect_file_stem, borzoi_effect_file_bins):
	dicti = {}

	for file_bin in borzoi_effect_file_bins:

		borzoi_variant_effect_file = borzoi_variant_effect_file_stem + file_bin + '.txt'

		f = open(borzoi_variant_effect_file)
		head_count = 0
		counter = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue


			bs_vals = np.asarray(data[5].split(';')).astype(float)

			borzoi_effect = np.mean(bs_vals)
			np_p_plus = (1 + np.sum(bs_vals >= 0))/(1+len(bs_vals))
			np_p_minus = (1 + np.sum(bs_vals <=0))/(1+len(bs_vals))
			borzoi_pvalue = 2.0*np.min((np_p_plus, np_p_minus))


			variant_id = data[2]
			gene_id = data[3]
			variant_gene_pair = variant_id + ':' + gene_id

			if variant_gene_pair in dicti:
				print('assumptin eornornor')
				pdb.set_trace()

			dicti[variant_gene_pair] = borzoi_pvalue
			
			counter = counter +1



		f.close()

	return dicti




def create_mapping_from_variant_gene_to_eqtl_sig(fm_eqtl_results_file, eqtl_sumstats_dir, tissue_name):
	# load in fm eqtl results
	fm_variant_gene_pairs = {}
	fm_variants = {}
	f = open(fm_eqtl_results_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] != tissue_name:
			continue
		variant_id = data[3]
		gene_id = data[7]
		fm_variants[variant_id] = 0
		fm_variant_gene_pairs[variant_id + ':' + gene_id] = 0
	f.close()


	# Load in non-sig eqtls
	dicti = {}
	for chrom_num in range(1,23):
		print(chrom_num)
		filer = eqtl_sumstats_dir + tissue_name + '.v10.allpairs.chr' + str(chrom_num) + '.parquet'
		pf = pq.ParquetFile(filer)
		for rg in range(pf.num_row_groups):
			table = pf.read_row_group(rg)   # this is a chunk

			# Process fields and print to output
			array = np.asarray(table)

			# Loop through snp-gene pairs 
			nrows = array.shape[0]
			for row_iter in range(nrows):
				data = array[row_iter, :]
				variant_id = data[1]
				geneid = data[0]
				af = float(data[3])

				sig_val = 0.0
				vg_pair = variant_id + ':' + geneid
				if vg_pair in fm_variant_gene_pairs:
					sig_val = 1.0

				
				#vg_pair1 = variant_id + ':' + geneid
				#vg_pair2 = alt_variant_id + ':' + geneid

				dicti[vg_pair] = sig_val
		print(len(dicti))

	return dicti



##########################
# Command line args
##########################
borzoi_results_file_stem= sys.argv[1]
processed_fm_eqtl_output_file = sys.argv[2]
eqtl_sumstats_dir = sys.argv[3]
gene_tss_file = sys.argv[4]
dist_to_tss_summary_file = sys.argv[5]
tissue_name = sys.argv[6]

# First create mapping from ensamble id to gene (chrom, tss, strand)
gene_tss_info = create_mapping_from_ensamble_id_to_gene_tss_info(gene_tss_file)



borzoi_effect_file_bins = np.asarray(['0','1'])



# First extract mapping from variant-gene pairs to borzoi significnace
variant_gene_to_borzoi_sig = create_mapping_from_variant_gene_pair_to_borzoi_sig(borzoi_results_file_stem, borzoi_effect_file_bins)

# Now create mapping from variant-gene pairs to eqtl sig
variant_gene_to_eqtl_sig = create_mapping_from_variant_gene_to_eqtl_sig(processed_fm_eqtl_output_file, eqtl_sumstats_dir, tissue_name)




t = open(dist_to_tss_summary_file,'w')
t.write('method\tvariant_id\tgene_id\tdist_to_tss\tsignificance\n')

for variant_gene_pair in [*variant_gene_to_borzoi_sig]:

	variant_id, gene_id = variant_gene_pair.split(':')

	if variant_gene_pair not in variant_gene_to_eqtl_sig:

		variant_info = variant_id.split('_')
		alt_variant_id = variant_info[0] + '_' + variant_info[1] + '_' + variant_info[3] + '_' + variant_info[2] + '_' + variant_info[4]
		alt_variant_gene_pair = alt_variant_id + ':' + gene_id
		if alt_variant_gene_pair in variant_gene_to_eqtl_sig:
			print('assumption eorroro')
			pdb.set_trace()

		continue

	gene_id = gene_id.split('.')[0]
	if gene_id not in gene_tss_info:
		continue

	gene_chrom, gene_tss, gene_strand = gene_tss_info[gene_id]

	var_pos = float(variant_id.split('_')[1])
	var_chrom = variant_id.split('_')[0]

	if var_chrom != gene_chrom:
		print('chrom mismatch error')
		continue

	if gene_strand == '+':
		distance = var_pos - gene_tss
	elif gene_strand == '-':
		distance = gene_tss - var_pos
	else:
		print('assumptioneorrornon on gene strand')
		pdb.set_trace()

	borzoi_sig = variant_gene_to_borzoi_sig[variant_gene_pair]
	eqtl_sig = variant_gene_to_eqtl_sig[variant_gene_pair]

	used_borzio = 0
	used_eqtl = 0
	if borzoi_sig < 0.05:
		used_borzio =1
		t.write('borzoi\t' + variant_id + '\t' + gene_id + '\t' + str(distance) + '\t' + str(borzoi_sig) + '\n')
	if eqtl_sig == 1:
		used_eqtl = 1
		t.write('eQTL\t' + variant_id + '\t' + gene_id + '\t' + str(distance) + '\t' + str(eqtl_sig) + '\n')
	
	if np.random.rand() < 0.01:
		if used_borzio == 0:
			t.write('borzoi\t' + variant_id + '\t' + gene_id + '\t' + str(distance) + '\t' + str(borzoi_sig) + '\n')
		if used_eqtl == 0:
			t.write('eQTL\t' + variant_id + '\t' + gene_id + '\t' + str(distance) + '\t' + str(eqtl_sig) + '\n')


t.close()
print(dist_to_tss_summary_file)




