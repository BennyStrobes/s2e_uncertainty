import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink
import time






def load_in_cis_variants(gene_cis_variant_file):
	f = open(gene_cis_variant_file)
	ordered_variant_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ordered_variant_names.append(data[1])
	f.close()

	return np.asarray(ordered_variant_names)


def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	mapping = {}

	n_snps = len(ordered_snps)
	for snp_iter in range(n_snps):
		snp_name = ordered_snps[snp_iter]

		if snp_name in mapping:
			print('asssumption erororo')
			pdb.set_trace()
		mapping[snp_name] = snp_iter

	return mapping


def standardize_geno(geno_mat):
	std_geno_mat = np.copy(geno_mat)
	n_snps = std_geno_mat.shape[0]

	for snp_iter in range(n_snps):
		std_geno_mat[snp_iter,:] = (std_geno_mat[snp_iter,:] - np.mean(std_geno_mat[snp_iter,:]))/np.std(std_geno_mat[snp_iter,:], ddof=1)

	return std_geno_mat


######################
# Command line args
######################
onek_genomes_plink_filestem = sys.argv[1]
gene_summary_file = sys.argv[2]
ld_dir = sys.argv[3]

# Open global gene summary file
ld_gene_summary_file = ld_dir + 'gene_summary_w_ld_full2.txt'
t_glob = open(ld_gene_summary_file,'w')
t_glob.write('gene_name\tchrom_num\tgene_tss\tn_cis_variants\tcis_variant_bim_file\tld_summary_file\traw_ld_file\traw_std_genotype_file\n')

# Loop through chromosomes
for chrom_num in range(1,23):
	print(chrom_num)
	# string of chromosome name
	chrom_string = 'chr' + str(chrom_num)

	# Load in chromosome plink data
	(bim, fam, G) = read_plink(onek_genomes_plink_filestem + str(chrom_num))

	# Create mapping from variant id to index
	rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))

	# Loop through genes / filter to only genes on this chromosome
	f = open(gene_summary_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_chrom_num = data[1]
		if line_chrom_num != chrom_string:
			continue
		gene_name = data[0]
		counter = counter + 1
		gene_cis_variant_file = data[4]
		n_cis_variants = int(data[3])
		if n_cis_variants < 11:
			continue
		cis_variant_names = load_in_cis_variants(gene_cis_variant_file)

		# Extract variant indices corresponding to cis_variant names
		cis_variant_indices = []
		for variant_name in cis_variant_names:
			cis_variant_indices.append(rsid_to_genotype_index[variant_name])
		cis_variant_indices = np.asarray(cis_variant_indices)

		# Extract genotype matrix
		geno_mat = G[cis_variant_indices,:].compute()
		n_samp = geno_mat.shape[1]
		
		# Compute LD
		LD = np.corrcoef(geno_mat)

		# Get squared LD
		squared_LD = np.square(LD)
		squared_LD_adj = squared_LD - ((1.0-squared_LD)/(n_samp-2.0))

		# Compute Ld-scores
		ld_scores = np.sum(squared_LD_adj,axis=1)
		mean_ld = np.mean(LD,axis=0)

		# Standardize genotype
		std_geno_mat = standardize_geno(geno_mat)

		# Print to output file
		# Output file has a row for each variant and column for ld score and average mean
		gene_ld_summary_file = ld_dir + gene_name + '_ld_summary.txt'
		t = open(gene_ld_summary_file,'w')
		t.write('variant_id\tld_score\tmean_ld\n')
		for var_iter, var_name in enumerate(cis_variant_names):
			t.write(var_name + '\t' + str(ld_scores[var_iter]) + '\t' + str(mean_ld[var_iter]) + '\n')
		t.close()

		# Raw LD
		gene_raw_ld_file = ld_dir + gene_name + '_raw_ld.npy'
		np.save(gene_raw_ld_file, LD)

		# Raw genotype
		gene_raw_geno_file = ld_dir + gene_name + '_raw_geno.npy'
		np.save(gene_raw_geno_file, np.transpose(std_geno_mat))		

		# update global file
		t_glob.write('\t'.join(data) + '\t' + gene_ld_summary_file + '\t' + gene_raw_ld_file + '\t' + gene_raw_geno_file + '\n')

	f.close()


t_glob.close()