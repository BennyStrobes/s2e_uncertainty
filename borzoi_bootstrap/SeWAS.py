import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink1_bin
import pickle
import argparse
import gzip











def print_SeWAS_bear():
	# Credit: https://github.com/ivolo/animals/blob/master/data/animals.txt
	print(r"""
   _,-""`""-~`)
(`~           \
 |     a   a   \
 ;        o     ; ___  _,,,,_     _.-~'.
  \      `^`    /`_.-"~      `~-;`      \
   \_      _  .'                 `,     |
	 |`-                           \'__/
	/                      ,_       \  `'-.
   /    .-""~~--.            `"-,   ;_    /
  |              \               \  | `""`
   \__.--'`"-.   /_               |'
			  `"`  `~~~---..,     |
 SeWAS                       \ _.-'`-.
							  \       \
							   '.     /
								 `"~"`
	""")
	return



def extract_dictionary_list_of_genes_from_file(protein_coding_genes_file, chrom_string):
	f = open(protein_coding_genes_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != chrom_string:
			continue
		dicti[data[3]] = 1

	f.close()
	return dicti

def extract_borzoi_effect_sizes(borzoi_effect_size_file, chrom_num, valid_genes):
	chrom_string = 'chr' + str(chrom_num)
	gene_obj = {}

	f = open(borzoi_effect_size_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count+1
			continue
		if chrom_string != data[0]:
			continue

		variant_name = data[2]
		gene_name = data[3]
		if gene_name not in valid_genes:
			continue

		mean_borzoi_effect = float(data[4])
		bs_borzoi_effects = np.asarray(data[5].split(';')).astype(float)

		bs_minor_count = np.min([np.sum(bs_borzoi_effects > 0.0), np.sum(bs_borzoi_effects < 0.0)])

		if gene_name not in gene_obj:
			gene_obj[gene_name] = []
		gene_obj[gene_name].append((variant_name, mean_borzoi_effect, np.std(bs_borzoi_effects), bs_minor_count, bs_borzoi_effects))

		counter = counter + 1

	f.close()

	return gene_obj


def load_in_gwas_sumstats(sumstat_file):
	if sumstat_file.endswith('gz'):
		f = gzip.open(sumstat_file,'rt')
	else:
		f = open(sumstat_file)


	variant_to_z = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 16:
			print('assumption eororo')
			pdb.set_trace()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			header = np.asarray(data)
			continue
		# Extract variant info from line
		chrom_num = data[1]
		var_pos = data[2]
		alt_allele = data[4]
		ref_allele = data[5]

		# Filter to snvs
		if len(ref_allele) != 1 or len(alt_allele) != 1:
			continue

		# Create variant ID (in gtex format)
		variant_id = 'chr' + chrom_num + '_' + var_pos + '_' + ref_allele + '_' + alt_allele + '_b38'

		# Get sumstat info
		chi_sq_stat = float(data[-2])
		ols_beta = float(data[-6])

		if chi_sq_stat <= 0.0:
			continue

		# Convert from chi-sq to z-score
		z_score = np.sqrt(chi_sq_stat)
		if ols_beta < 0.0:
			z_score = z_score*-1.0
		
		# Add to global dictionary
		variant_to_z[variant_id] = z_score

	f.close()

	return variant_to_z

def load_in_genotype_data_for_this_chromosome(plink_genotype_stem, chrom_num):
	G_obj = read_plink1_bin(plink_genotype_stem + str(chrom_num) + ".bed", plink_genotype_stem + str(chrom_num) +".bim", plink_genotype_stem + str(chrom_num) +".fam", verbose=False)
	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	G_obj_gtex_ids_dicti = {}
	gtex_variant_to_hm3_bool = {}
	for ii, posser in enumerate(G_obj_pos):
		tmp_gtex_id = 'chr' + G_obj_chrom[ii] + '_' + str(posser) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii] + '_b38'
		G_obj_gtex_ids_dicti[tmp_gtex_id] = ii

	return G_obj_geno, G_obj_gtex_ids_dicti

def extract_borzoi_effects_for_gene(list_of_tuples):
	variants = []
	betas = []
	beta_ses = []
	beta_mab = []
	beta_bs = []
	for tupler in list_of_tuples:
		variants.append(tupler[0])
		betas.append(tupler[1])
		beta_ses.append(tupler[2])
		beta_mab.append(tupler[3])
		beta_bs.append(tupler[4])

	return np.asarray(variants), np.asarray(betas), np.asarray(beta_bs)


def extract_gwas_zeds_for_variant_array(gene_borzoi_variants, variant_to_gwas_z):
	zeds = []
	for variant_id in gene_borzoi_variants:
		var_info = variant_id.split('_')
		alt_variant_id = var_info[0] + '_' + var_info[1] + '_' + var_info[3] + '_' + var_info[2] + '_b38'

		if variant_id in variant_to_gwas_z:
			zed_score = variant_to_gwas_z[variant_id]
		elif alt_variant_id in variant_to_gwas_z:
			zed_score = variant_to_gwas_z[alt_variant_id]
		else: # Missing
			zed_score = np.nan
		zeds.append(zed_score)
	zeds = np.asarray(zeds)

	return zeds


def extract_boolean_on_whether_we_have_genotype_data_for_each_variant(gene_borzoi_variants, genotype_variant_to_position):
	boolers = []

	for variant_id in gene_borzoi_variants:
		var_info = variant_id.split('_')
		alt_variant_id = var_info[0] + '_' + var_info[1] + '_' + var_info[3] + '_' + var_info[2] + '_b38'

		if variant_id in genotype_variant_to_position or alt_variant_id in genotype_variant_to_position:
			boolers.append(True)
		else:
			boolers.append(False)

	return np.asarray(boolers)


def extract_genotype_matrix_for_subset_of_snps(G_obj_geno, G_obj_gtex_ids_dicti, gene_eqtl_variants):
	geno = []
	mafs = []
	for variant_id in gene_eqtl_variants:
		var_info = variant_id.split('_')
		alt_variant_id = var_info[0] + '_' + var_info[1] + '_' + var_info[3] + '_' + var_info[2] + '_b38'

		if variant_id in G_obj_gtex_ids_dicti:
			tmp_geno = G_obj_geno[:, G_obj_gtex_ids_dicti[variant_id]]
			geno.append(tmp_geno)
			mafs.append(np.sum(tmp_geno)/(2.0*len(tmp_geno)))
		elif alt_variant_id in G_obj_gtex_ids_dicti:
			tmp_geno = G_obj_geno[:, G_obj_gtex_ids_dicti[alt_variant_id]]
			geno.append(-tmp_geno)
			mafs.append(np.sum(tmp_geno)/(2.0*len(tmp_geno)))
		else:
			print('assumption erroror')
			pdb.set_trace()

	geno = np.asarray(geno)
	mafs = np.asarray(mafs)

	return geno, mafs

def compute_twas_zed(gene_gwas_zeds, std_gene_borzoi_effects, R_mat):
	"""
	Compute a TWAS-style gene Z score:
		Z_TWAS = (w^T z) / sqrt(w^T R w)

	Parameters
	----------
	gene_gwas_zeds : array-like, shape (K,)
		SNP-level GWAS Z-scores for the K variants in the gene window.
	std_gene_borzoi_effects : array-like, shape (K,)
		Weights w for the same K variants (often standardized).
	R_mat : array-like, shape (K, K)
		LD correlation/covariance matrix among the K variants.

	Returns
	-------
	z_twas : float
		TWAS-style gene Z score.
	"""
	z = np.asarray(gene_gwas_zeds, dtype=float).reshape(-1)
	w = np.asarray(std_gene_borzoi_effects, dtype=float).reshape(-1)
	R = np.asarray(R_mat, dtype=float)

	if z.shape != w.shape:
		raise ValueError(f"gene_gwas_zeds and std_gene_borzoi_effects must have same shape. "
						 f"Got {z.shape} vs {w.shape}")
	if R.ndim != 2 or R.shape[0] != R.shape[1]:
		raise ValueError(f"R_mat must be square (KxK). Got {R.shape}")
	if R.shape[0] != z.shape[0]:
		raise ValueError(f"R_mat dimension must match length K. Got R {R.shape} vs K={z.shape[0]}")

	# Drop any SNPs with missing z or w (and subset R accordingly)
	mask = np.isfinite(z) & np.isfinite(w)
	if not np.all(mask):
		idx = np.where(mask)[0]
		if idx.size == 0:
			return np.nan
		z = z[idx]
		w = w[idx]
		R = R[np.ix_(idx, idx)]

	# Numerator: w^T z
	num = float(np.dot(w, z))

	# Denominator: sqrt(w^T R w)
	den_sq = float(np.dot(w, R @ w))

	# Numerical safety
	if not np.isfinite(den_sq) or den_sq <= 0.0:
		return np.nan

	return num / np.sqrt(den_sq)


def compute_twas_zed_bs(gene_gwas_zeds, std_gene_borzoi_effects_bs, R_mat):
	"""
	Compute TWAS-style gene Z scores for B bootstrap weight vectors.

	For each bootstrap b:
		Z_b = (w_b^T z) / sqrt(w_b^T R w_b)

	Parameters
	----------
	gene_gwas_zeds : array-like, shape (K,)
		SNP-level GWAS Z-scores for the K variants.
	std_gene_borzoi_effects_bs : array-like, shape (K, B)
		Bootstrap weight matrix. Each column is a weight vector w_b.
	R_mat : array-like, shape (K, K)
		LD correlation/covariance matrix among the K variants.

	Returns
	-------
	z_twas_bs : ndarray, shape (B,)
		TWAS Z score for each bootstrap.
	"""
	z = np.asarray(gene_gwas_zeds, dtype=float).reshape(-1)
	W = np.asarray(std_gene_borzoi_effects_bs, dtype=float)
	R = np.asarray(R_mat, dtype=float)

	if W.ndim != 2:
		raise ValueError("std_gene_borzoi_effects_bs must be a K x B matrix")

	K, B = W.shape

	if z.shape[0] != K:
		raise ValueError(f"Length of gene_gwas_zeds ({z.shape[0]}) "
						 f"must match K ({K})")

	if R.shape != (K, K):
		raise ValueError(f"R_mat must be shape ({K}, {K}), got {R.shape}")

	# Drop SNPs with missing z or any missing weights across bootstraps
	mask = np.isfinite(z) & np.all(np.isfinite(W), axis=1)
	if not np.all(mask):
		idx = np.where(mask)[0]
		if idx.size == 0:
			return np.full(B, np.nan)
		z = z[idx]
		W = W[idx, :]
		R = R[np.ix_(idx, idx)]

	# Numerator: w_b^T z  -> shape (B,)
	num = W.T @ z

	# Denominator: sqrt(diag(W^T R W)) -> shape (B,)
	RW = R @ W                  # (K, B)
	den_sq = np.sum(W * RW, axis=0)

	# Numerical safety
	z_twas_bs = np.full(B, np.nan)
	valid = np.isfinite(den_sq) & (den_sq > 0)
	z_twas_bs[valid] = num[valid] / np.sqrt(den_sq[valid])

	return z_twas_bs

def rank_based_pvalue(bs_corry):
	np_p_plus = (1 + np.sum(bs_corry >= 0))/(1+len(bs_corry))
	np_p_minus = (1 + np.sum(bs_corry <=0))/(1+len(bs_corry))
	borzoi_pvalues = 2.0*np.min((np_p_plus, np_p_minus))
	return borzoi_pvalues

def compute_uncertainty_weighted_twas_zed(gene_gwas_zeds, std_gene_borzoi_effects_bs, R_mat):
	RW = R_mat @ std_gene_borzoi_effects_bs                  # (K, B)
	sampled_variances = np.sum(std_gene_borzoi_effects_bs * RW, axis=0)
	total_var = np.mean(sampled_variances)

	meany = np.mean(std_gene_borzoi_effects_bs,axis=1)

	z_twas_tot_var = np.dot(gene_gwas_zeds, meany)/np.sqrt(total_var)
	

	return z_twas_tot_var

################################################################################
# main
################################################################################
def main():
	######################
	# Command line args
	######################
	# Necessary
	parser = argparse.ArgumentParser()
	parser.add_argument('--sumstat-file', default='', type=str,
						help='File containing GWAS summary statistics')
	parser.add_argument('--plink-genotype-stem', default='', type=str,
						help='Path to plink genotype stem')
	parser.add_argument('--borzoi-result-file', default='', type=str,
						help='Path to File name containing bootstrapped borzoi results')
	parser.add_argument('--gene-list-file', default='', type=str,
						help='Path to File name containing list of genes to run on')
	parser.add_argument('--output-file', default='', type=str,
						help='Path to output file stem')

	# Defaults
	parser.add_argument('--min-snps-per-gene', default=10, type=int,
						help='Minimum number of snps per gene. else we throw out the gene.')
	args = parser.parse_args()

	print_SeWAS_bear()

	# Open output file handle
	t = open(args.output_file,'w')
	t.write('gene_id\tchrom_num\tTWAS_z\tbootstrapped_cov_pvalue\tUncertainty_TWAS_Z\n')

	# Load in GWAS sumstat z-scores into dictionary mapping from variant to z-score
	variant_to_gwas_z = load_in_gwas_sumstats(args.sumstat_file)

	# Loop through chromosomes
	for chrom_num in range(1,5):

		print(chrom_num)

		# Extract dictionary list of genes to use on this chromosome
		valid_genes = extract_dictionary_list_of_genes_from_file(args.gene_list_file, 'chr' + str(chrom_num))

		# Extract borzoi effect sizes and uncertainties into gene-based dictionary data structure
		borzoi_effect_size_gene_obj = extract_borzoi_effect_sizes(args.borzoi_result_file, chrom_num, valid_genes)

		# Load in genotype data for this chromosome
		G_obj_geno, genotype_variant_to_position = load_in_genotype_data_for_this_chromosome(args.plink_genotype_stem, chrom_num)

		# Loop through genes
		for geneid in [*borzoi_effect_size_gene_obj]:

			# Extract borzoi effects for this gene
			gene_borzoi_variants, gene_borzoi_effects, gene_borzoi_effect_bs = extract_borzoi_effects_for_gene(borzoi_effect_size_gene_obj[geneid])

			# Extract GWAS z-scores for these variants corresponding to this gene
			# return nans for missing variants
			gene_gwas_zeds = extract_gwas_zeds_for_variant_array(gene_borzoi_variants, variant_to_gwas_z)

			# Extract boolean vector represnting whether we have genotype data for each variant
			gene_genotype_boolean = extract_boolean_on_whether_we_have_genotype_data_for_each_variant(gene_borzoi_variants, genotype_variant_to_position)
	
			# Get variant indices
			valid_variant_indices = (np.isnan(gene_gwas_zeds) == False) & (gene_genotype_boolean)

			# SKip genes with too few snps
			if np.sum(valid_variant_indices) < args.min_snps_per_gene:
				continue

			# Filter stuff to valid variants
			gene_borzoi_variants = gene_borzoi_variants[valid_variant_indices]
			gene_borzoi_effects = gene_borzoi_effects[valid_variant_indices]
			gene_borzoi_effect_bs = gene_borzoi_effect_bs[valid_variant_indices, :]
			gene_gwas_zeds = gene_gwas_zeds[valid_variant_indices]

			# Extract genotype matrix for gene-snps
			gene_geno, gene_geno_maf = extract_genotype_matrix_for_subset_of_snps(G_obj_geno, genotype_variant_to_position, gene_borzoi_variants)

			# Get borzoi effects in standardized genotype space
			std_gene_borzoi_effects = gene_borzoi_effects*np.sqrt(2.0*gene_geno_maf*(1.0-gene_geno_maf))
			std_gene_borzoi_effects_bs = np.transpose(np.transpose(gene_borzoi_effect_bs)*np.sqrt(2.0*gene_geno_maf*(1.0-gene_geno_maf)))


			# Compute LD Matrix
			R_mat = np.corrcoef(gene_geno)

			# Compute TWAS z-scores
			twas_zed = compute_twas_zed(gene_gwas_zeds, std_gene_borzoi_effects, R_mat)
			
			# Compute TWAS z-scores for each bootstrapped sample
			twas_zed_bs = compute_twas_zed_bs(gene_gwas_zeds, std_gene_borzoi_effects_bs, R_mat)
			
			# Compute rank based pvalue of covariance not including zero
			bs_cov_pvalue = rank_based_pvalue(twas_zed_bs)


			# Compute uncertainty_weighted_twas
			uncertainty_weighted_twas_zed = compute_uncertainty_weighted_twas_zed(gene_gwas_zeds, std_gene_borzoi_effects_bs, R_mat)



			# Print to output
			t.write(geneid + '\t' + str(chrom_num) + '\t' + str(twas_zed) + '\t' + str(bs_cov_pvalue) + '\t' + str(uncertainty_weighted_twas_zed) + '\n')
			t.flush()
	t.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
	main()