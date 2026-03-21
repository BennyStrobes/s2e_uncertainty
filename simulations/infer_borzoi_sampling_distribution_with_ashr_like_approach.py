import numpy as np
import os
import sys
import pdb
import gzip
from pandas_plink import read_plink
from scipy.special import logsumexp




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


def extract_ldsc_vectors(gene_ld_summary_file, onek_genomes_plink_filestem, gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, use_adjusted_LD=False):
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
			if use_adjusted_LD:
				ld_scores = np.sum(squared_LD_adj,axis=1)
			else:
				ld_scores = np.sum(squared_LD,axis=1)
			
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


def fit_ldscore_scaled_error_mixture(
	eqtl_marginal_effect_sizes,
	borzoi_marginal_effect_sizes,
	ldscores,
	eqtl_sample_size,
	omega2_grid=None,
	max_iter=500,
	tol=1e-8,
	min_pi=1e-12,
	verbose=True,
):
	"""
	Fit a working variance-mixture model to residuals

		r_k = eqtl_marginal_effect_sizes[k] - borzoi_marginal_effect_sizes[k]

	using the likelihood

		r_k ~ sum_h pi_h * N(0, 1/N + ldscore_k * omega2_h)

	where omega2_h is a fixed variance grid and pi_h are learned by EM.

	Parameters
	----------
	eqtl_marginal_effect_sizes : np.ndarray, shape (M,)
	borzoi_marginal_effect_sizes : np.ndarray, shape (M,)
	ldscores : np.ndarray, shape (M,)
	eqtl_sample_size : float or int
		GWAS / eQTL sample size giving noise variance 1/N.
	omega2_grid : np.ndarray, optional
		Fixed variance grid for eta_k = borzoi_hat_k - beta_k.
		If None, a default log-spaced grid is created.
	max_iter : int
	tol : float
	min_pi : float
		Small floor for mixture weights to avoid exact zeros.
	verbose : bool

	Returns
	-------
	results : dict with keys
		residuals
		omega2_grid
		pi
		posterior_resp
		log_likelihood_trace
		mean_error_variance
		fitted_variance_per_snp
	"""
	eqtl_marginal_effect_sizes = np.asarray(eqtl_marginal_effect_sizes)
	borzoi_marginal_effect_sizes = np.asarray(borzoi_marginal_effect_sizes)
	ldscores = np.asarray(ldscores)

	if eqtl_marginal_effect_sizes.ndim != 1:
		raise ValueError("eqtl_marginal_effect_sizes must be 1D")
	if borzoi_marginal_effect_sizes.ndim != 1:
		raise ValueError("borzoi_marginal_effect_sizes must be 1D")
	if ldscores.ndim != 1:
		raise ValueError("ldscores must be 1D")
	if len(eqtl_marginal_effect_sizes) != len(borzoi_marginal_effect_sizes):
		raise ValueError("effect size arrays must have same length")
	if len(eqtl_marginal_effect_sizes) != len(ldscores):
		raise ValueError("effect size arrays and ldscores must have same length")
	if np.any(ldscores < 0):
		raise ValueError("ldscores must be nonnegative")

	# Working residual
	residuals = eqtl_marginal_effect_sizes - borzoi_marginal_effect_sizes
	M = residuals.shape[0]
	noise_var = 1.0 / float(eqtl_sample_size)

	# Default grid
	if omega2_grid is None:
		# crude data-driven upper scale
		approx_total_var = max(np.mean(residuals**2) - noise_var, 1e-12)
		upper = max(approx_total_var, 1e-8)
		omega2_grid = np.concatenate([
			np.array([0.0]),
			np.exp(np.linspace(np.log(upper * 1e-4), np.log(upper * 1e2), 30))
		])
	omega2_grid = np.asarray(omega2_grid, dtype=float)

	if omega2_grid.ndim != 1:
		raise ValueError("omega2_grid must be 1D")
	if np.any(omega2_grid < 0):
		raise ValueError("omega2_grid must be nonnegative")

	H = len(omega2_grid)

	# Initialize mixture weights uniformly
	pi = np.ones(H) / H
	loglik_trace = []

	# Precompute pieces with broadcasting:
	# var_mat[k, h] = 1/N + ldscore[k] * omega2_grid[h]
	var_mat = noise_var + np.outer(ldscores, omega2_grid)

	# numerical safety
	var_mat = np.maximum(var_mat, 1e-300)

	for it in range(max_iter):
		old_pi = pi.copy()

		# E-step:
		# log p(r_k | h) + log pi_h
		# Normal(0, var): -0.5*(log(2pi*var) + r^2/var)
		log_comp = (
			np.log(pi[None, :]) -
			0.5 * (np.log(2.0 * np.pi * var_mat) + (residuals[:, None] ** 2) / var_mat)
		)

		log_denom = logsumexp(log_comp, axis=1)
		resp = np.exp(log_comp - log_denom[:, None])

		# M-step:
		pi = resp.mean(axis=0)

		# avoid exact zeros
		pi = np.maximum(pi, min_pi)
		pi = pi / np.sum(pi)

		loglik = np.sum(log_denom)
		loglik_trace.append(loglik)

		if verbose and (it % 10 == 0 or it < 5):
			mean_error_var = np.sum(pi * omega2_grid)
			print(
				f"iter={it:4d}  loglik={loglik:.6f}  "
				f"E[error_var]={mean_error_var:.6e}"
			)

		# convergence
		if np.max(np.abs(pi - old_pi)) < tol:
			if verbose:
				print(f"Converged at iteration {it}")
			break

	mean_error_variance = np.sum(pi * omega2_grid)
	fitted_variance_per_snp = noise_var + ldscores * mean_error_variance

	return {
		"residuals": residuals,
		"omega2_grid": omega2_grid,
		"pi": pi,
		"posterior_resp": resp,
		"log_likelihood_trace": np.asarray(loglik_trace),
		"mean_error_variance": mean_error_variance,
		"fitted_variance_per_snp": fitted_variance_per_snp,
	}



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
eqtl_marginal_effect_sizes, borzoi_marginal_effect_sizes, ldscores = extract_ldsc_vectors(gene_ld_summary_file, onek_genomes_plink_filestem, gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, use_adjusted_LD=False)
#ldscores[ldscores<=0.0] = 0.0


# Run regression
deltas = eqtl_marginal_effect_sizes - borzoi_marginal_effect_sizes


results = fit_ldscore_scaled_error_mixture(
	eqtl_marginal_effect_sizes=eqtl_marginal_effect_sizes,
	borzoi_marginal_effect_sizes=borzoi_marginal_effect_sizes,
	ldscores=ldscores,
	eqtl_sample_size=eqtl_sample_size,
	omega2_grid=None,   # or pass your own fixed grid
	max_iter=500,
	tol=1e-8,
	verbose=True,
)


t = open(results_file,'w')
t.write('omega2\tpi\n')
for omega2, pi in zip(results["omega2_grid"], results["pi"]):
	t.write(str(omega2) + '\t' + str(pi) + '\n')
t.close()



'''
y = np.square(deltas) - (1.0/eqtl_sample_size)
x = ldscores

tau2 = np.sum(x * y) / np.sum(x * x)
t = open(results_file,'w')
t.write('borzoi_sampling_variance\n')
t.write(str(tau2) + '\n')
t.close()
'''


'''
pis = np.asarray([.4, .6])
sig_sqs = np.asarray([.1, .9])
scaler = sampling_variance/np.sum(pis*sig_sqs)
sig_sqs = sig_sqs*scaler
'''


