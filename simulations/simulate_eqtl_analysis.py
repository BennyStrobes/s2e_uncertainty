import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink
import gzip





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


def create_mapping_from_gene_id_to_causal_effects(causal_variant_gene_effect_size_file):
	f = gzip.open(causal_variant_gene_effect_size_file,'rt')
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

def marginal_snp_effects_and_se(gene_expression, geno_mat):
	"""
	Compute marginal SNP effect sizes and standard errors from univariate
	linear regressions with an intercept.

	Parameters
	----------
	gene_expression : ndarray, shape (n_samples,)
		Expression vector.
	geno_mat : ndarray, shape (n_snps, n_samples)
		Genotype matrix.

	Returns
	-------
	betas : ndarray, shape (n_snps,)
		Marginal OLS effect size for each SNP.
	ses : ndarray, shape (n_snps,)
		Standard error of each marginal effect size.
	"""
	y = np.asarray(gene_expression, dtype=float)
	G = np.asarray(geno_mat, dtype=float)

	if y.ndim != 1:
		raise ValueError("gene_expression must be 1D")
	if G.ndim != 2:
		raise ValueError("geno_mat must be 2D")
	if G.shape[1] != y.shape[0]:
		raise ValueError("geno_mat must have shape (n_snps, n_samples)")

	n = y.shape[0]
	if n <= 2:
		raise ValueError("Need at least 3 samples")

	# Center y and each SNP genotype vector to account for intercept
	y_c = y - np.mean(y)
	G_c = G - np.mean(G, axis=1, keepdims=True)

	# Sxx_j = sum_i (g_ji - mean_g_j)^2
	sxx = np.sum(G_c**2, axis=1)

	# Identify monomorphic / zero-variance SNPs
	valid = sxx > 0

	betas = np.full(G.shape[0], np.nan)
	ses = np.full(G.shape[0], np.nan)

	if not np.any(valid):
		return betas, ses

	# beta_j = sum_i Gc_ji * yc_i / sum_i Gc_ji^2
	numer = G_c[valid] @ y_c
	betas_valid = numer / sxx[valid]

	# Residual sum of squares for each SNP:
	# RSS_j = sum(y_c^2) - beta_j^2 * Sxx_j
	# because with intercept and centered variables, fitted values are beta_j * G_cj
	y_ss = np.sum(y_c**2)
	rss = y_ss - (betas_valid**2) * sxx[valid]

	# Numerical guard against tiny negative values from floating point
	rss = np.maximum(rss, 0.0)

	sigma2 = rss / (n - 2)
	ses_valid = np.sqrt(sigma2 / sxx[valid])

	betas[valid] = betas_valid
	ses[valid] = ses_valid

	return betas, ses

######################
# Command line args
######################
causal_variant_gene_effect_size_file = sys.argv[1]
est_eqtl_effect_size_file = sys.argv[2]
gene_ld_summary_file = sys.argv[3]
onek_genomes_plink_filestem = sys.argv[4]
eqtl_sample_size = int(sys.argv[5])
simulation_iter = int(sys.argv[6])



# set random seed
np.random.choice(simulation_iter)

# Create mapping from gene id to vector of causal effects
gene_id_to_causal_effects = create_mapping_from_gene_id_to_causal_effects(causal_variant_gene_effect_size_file)


# Get sample indices
sample_indices = np.random.choice(np.arange(489), size=eqtl_sample_size, replace=False)


# Open output file handle
t = gzip.open(est_eqtl_effect_size_file, 'wt')
t.write('gene\tvariant\teqtl_effect_size\teqtl_effect_size_se\n')

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
		geno_mat = geno_mat[:, sample_indices]
		n_samp = geno_mat.shape[1]

		# Standardize genotype matrix
		std_geno_mat = standardize_geno(geno_mat)

		# extract causal effects for gene
		causal_effect_tupler = gene_id_to_causal_effects[gene_name]
		tmp_variant_names = []
		tmp_causal_effects = []
		for tuply in causal_effect_tupler:
			tmp_variant_names.append(tuply[0])
			tmp_causal_effects.append(tuply[1])
		tmp_variant_names = np.asarray(tmp_variant_names)
		tmp_causal_effects = np.asarray(tmp_causal_effects)

		# QUick error check
		if np.array_equal(cis_variant_names, tmp_variant_names) == False:
			print('assumption erororo')
			pdb.set_trace()

		# Get genetic component of gene expression
		genetic_gene = np.dot(np.transpose(std_geno_mat), tmp_causal_effects)

		# Simulate gene expression
		gene_var = np.var(genetic_gene,ddof=1)
		if gene_var > 1:
			print(gene_var)
			resid_var = 1e-5
		else:
			resid_var = 1.0 - gene_var
		gene_expression = np.random.normal(loc=genetic_gene, scale=np.sqrt(resid_var))

		# Run eqtl mapping
		eqtl_beta, eqtl_beta_se = marginal_snp_effects_and_se(gene_expression, geno_mat)

		# Print to output file
		for var_iter, variant_id in enumerate(cis_variant_names):
			t.write(gene_name + '\t' + variant_id + '\t' + str(eqtl_beta[var_iter]) + '\t' + str(eqtl_beta_se[var_iter]) + '\n')

	f.close()
t.close()

print(est_eqtl_effect_size_file)
