import numpy as np
import os
import sys
import pdb
import gzip
from pandas_plink import read_plink




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





def create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file):
	f = gzip.open(est_borzoi_effect_size_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		var_id = data[1]
		effect = float(data[2])
		if gene_id not in mapping:
			mapping[gene_id] = []
		mapping[gene_id].append((var_id, effect))
	f.close()
	return mapping

def create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file, eqtl_sample_size):
	desired_se = np.sqrt(1.0/eqtl_sample_size)
	f = gzip.open(est_eqtl_effect_size_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		var_id = data[1]
		effect = float(data[2])
		se = float(data[3])
		if gene_id not in mapping:
			mapping[gene_id] = []
		zed = effect/se
		std_effect = zed*desired_se
		mapping[gene_id].append((var_id, std_effect))
	f.close()
	return mapping

def load_in_effects(tupler):
	var_arr = []
	effect_arr = []

	for tuply in tupler:
		var_arr.append(tuply[0])
		effect_arr.append(tuply[1])



	return np.asarray(effect_arr), np.asarray(var_arr)


def extract_ldsc_vectors(gene_ld_summary_file, onek_genomes_plink_filestem, gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects):
	marginal_eqtl_vec = []
	marginal_borzoi_vec = []
	ld_score_vec = []

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
		f = open(gene_ld_summary_file)
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

			cis_variant_names = load_in_cis_variants(gene_cis_variant_file)

			# Extract variant indices corresponding to cis_variant names
			cis_variant_indices = []
			for variant_name in cis_variant_names:
				cis_variant_indices.append(rsid_to_genotype_index[variant_name])
			cis_variant_indices = np.asarray(cis_variant_indices)

			# Extract genotype matrix
			geno_mat = G[cis_variant_indices,:].compute()
			n_samp = geno_mat.shape[1]

			# Create LD matrix
			LD = np.corrcoef(geno_mat)

			# Get squared LD
			squared_LD = np.square(LD)
			squared_LD_adj = squared_LD - ((1.0-squared_LD)/(n_samp-2.0))

			# Compute Ld-scores
			ld_scores = np.sum(squared_LD_adj,axis=1)
			
			# Load in est eqtl effect sizes
			est_eqtl_effect_sizes, tmp_variant_names1 = load_in_effects(gene_id_to_est_eqtl_effects[gene_name])

			# Load in est borzoi causal effect sizes
			est_borzoi_causal_effect_sizes, tmp_variant_names2 = load_in_effects(gene_id_to_est_borzoi_effects[gene_name])
			
			# Quick error checking
			if np.array_equal(tmp_variant_names1,cis_variant_names) == False or np.array_equal(tmp_variant_names2,cis_variant_names) == False:
				print('variant mismatch error')
				pdb.set_trace()

			# Get borzoi predicted marginal effect sizes
			est_borzoi_marginal_effects = np.dot(LD, est_borzoi_causal_effect_sizes)

			for var_iter, variant_name in enumerate(cis_variant_names):
				marginal_eqtl_vec.append(est_eqtl_effect_sizes[var_iter])
				marginal_borzoi_vec.append(est_borzoi_marginal_effects[var_iter])
				ld_score_vec.append(ld_scores[var_iter])

		f.close()
	return np.asarray(marginal_eqtl_vec), np.asarray(marginal_borzoi_vec), np.asarray(ld_score_vec)


#####################
# Command line args
#####################
est_borzoi_effect_size_file = sys.argv[1]
est_eqtl_effect_size_file = sys.argv[2]
gene_ld_summary_file = sys.argv[3]
onek_genomes_plink_filestem = sys.argv[4]
eqtl_sample_size = float(sys.argv[5])
results_file = sys.argv[6]


# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file)

# Create mapping from gene id to vector of est (standardized) eqtl effect sizes
gene_id_to_est_eqtl_effects = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file, eqtl_sample_size)


# Extract big vector (concatented across genes, variants) of standardized eqtl effects size, marginal pred borzoi effect sizes, and ldscores
eqtl_marginal_effect_sizes, borzoi_marginal_effect_sizes, ldscores = extract_ldsc_vectors(gene_ld_summary_file, onek_genomes_plink_filestem, gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects)


# Run regression
deltas = eqtl_marginal_effect_sizes - borzoi_marginal_effect_sizes
y = np.square(deltas) - (1.0/eqtl_sample_size)
x = ldscores

tau2 = np.sum(x * y) / np.sum(x * x)
t = open(results_file,'w')
t.write('borzoi_sampling_variance\n')
t.write(str(tau2) + '\n')
t.close()








