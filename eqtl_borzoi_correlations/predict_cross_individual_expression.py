import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time
from scipy import stats
from scipy.optimize import minimize_scalar


try:
	from numba import njit
except ImportError:
	def njit(*args, **kwargs):
		def decorator(func):
			return func
		return decorator







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
		chrom_num = data[2]
		snp_pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		effect = float(data[6])

		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()

		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, effect)
	f.close()
	return mapping



def create_mapping_from_gene_id_to_variant_gene_annotations(sim_variant_gene_annotation_file):
	f = gzip.open(sim_variant_gene_annotation_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[6:])
			continue
		gene_id = data[0]
		var_id = data[1]
		chrom_num = data[2]
		snp_pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		anno = np.asarray(data[6:]).astype(float)
		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, anno)
	f.close()
	return mapping, anno_names



def create_mapping_from_gene_id_to_expression_vector(expr_file):
	dicti = {}
	head_count = 0
	f = open(expr_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[3].split('.')[0]
		expr = np.asarray(data[4:]).astype(float)
		if gene_id in dicti:
			print('assumption erororr')
			pdb.set_trace()
		dicti[gene_id] = expr
	f.close()
	return dicti

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

def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
	if len(snp_array) != len(a0_arr):
		print('assumption eorroro')
		pdb.set_trace()
	if len(snp_array) != len(a1_arr):
		print('assumption eorroro')
		pdb.set_trace()

	dicti = {}

	for ii, snp_id in enumerate(snp_array):
		if snp_id in dicti:
			print('assumpationoenroer')
			pdb.set_trace()
		dicti[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
	return dicti

def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num

def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_borzoi_effects):
	unique_vars = np.unique([*var_to_est_borzoi_effects])
	final_vars = []
	for var in unique_vars:
		if var not in rsid_to_genotype_index:
			continue
		geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])

		passing = True
		if var in var_to_est_borzoi_effects:
			borzoi_alleles = var_to_est_borzoi_effects[var][4:6]
			if set(geno_alleles) != set(borzoi_alleles):
				passing = False
		if passing == False:
			continue
		final_vars.append(var)
	return np.asarray(final_vars)

def load_in_snp_gene_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.nan)
			alleles.append(('nan', 'nan'))
			print('assumption erororr')
			pdb.set_trace()
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles)

def load_in_snp_gene_anno_data(ordered_cis_variants, var_to_variant_gene_anno, n_anno):
	annos = []
	alleles = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_variant_gene_anno:
			annos.append(np.full(n_anno, np.nan))
			alleles.append(('nan', 'nan'))
			print('assumptioneorrnorn')
			pdb.set_trace()
		else:
			var_info = var_to_variant_gene_anno[variant_id]
			annos.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.vstack(annos), np.asarray(alleles)

def sample_causal_effects_from_gaussian_distribution(borzoi_vec, anno, distribution_obj, n_samp=100):
	if len(borzoi_vec) != anno.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	avg_mu = distribution_obj['avg_mu']
	avg_var = distribution_obj['avg_var']
	bin_indices = distribution_obj['bin_indices']

	if anno.shape[1] != len(avg_mu):
		print('assumption eroror')
		pdb.set_trace()
	if anno.shape[1] != len(avg_var):
		print('assumption eroror')
		pdb.set_trace()
	if np.array_equal(bin_indices, np.arange(len(avg_mu))) == False:
		print('assumption eroror')
		pdb.set_trace()

	snp_bins = np.argmax(anno, axis=1)
	if np.all(np.sum(anno, axis=1) == 1.0) == False:
		print('assumption eroror')
		pdb.set_trace()

	per_snp_mu = avg_mu[snp_bins]*borzoi_vec
	per_snp_sd = np.sqrt(avg_var[snp_bins])

	sampled_betas = np.random.normal(loc=per_snp_mu[None, :], scale=per_snp_sd[None, :], size=(n_samp, len(borzoi_vec)))
	return sampled_betas

def create_continuous_piecewise_linear_squared_basis(delta, knots):
	delta_sq = np.square(np.asarray(delta))
	basis = [np.ones(len(delta_sq)), delta_sq]
	for knot in knots:
		basis.append(np.maximum(delta_sq - knot, 0.0))
	return np.transpose(np.vstack(basis))

def create_continuous_piecewise_linear_basis(abs_delta, knots):
	abs_delta = np.asarray(abs_delta)
	basis = [np.ones(len(abs_delta)), abs_delta]
	for knot in knots:
		basis.append(np.maximum(abs_delta - knot, 0.0))
	return np.transpose(np.vstack(basis))

def compute_ldscore_grid_prior_variance(borzoi_vec, distribution_obj):
	grid_basis = create_continuous_piecewise_linear_basis(np.abs(borzoi_vec), distribution_obj['grid_knots'])
	return np.dot(grid_basis, distribution_obj['grid_coefs'])

def compute_ldscore_grid_squared_prior_variance(borzoi_vec, distribution_obj):
	grid_basis = create_continuous_piecewise_linear_squared_basis(borzoi_vec, distribution_obj['grid_knots'])
	return np.dot(grid_basis, distribution_obj['grid_coefs'])

def sample_causal_effects_from_ldscore_grid_distribution(borzoi_vec, distribution_obj, n_samp=100):
	per_snp_mu = distribution_obj['a_prior']*borzoi_vec
	per_snp_var = compute_ldscore_grid_prior_variance(borzoi_vec, distribution_obj)
	if np.any(per_snp_var < 0.0):
		print('warning: ldscore_grid predicted negative prior variances; flooring at 0.0 for sampling')
	per_snp_sd = np.sqrt(np.maximum(per_snp_var, 0.0))
	sampled_betas = np.random.normal(loc=per_snp_mu[None, :], scale=per_snp_sd[None, :], size=(n_samp, len(borzoi_vec)))
	return sampled_betas

def sample_causal_effects_from_ldscore_grid_squared_distribution(borzoi_vec, distribution_obj, n_samp=100):
	per_snp_mu = distribution_obj['a_prior']*borzoi_vec
	per_snp_var = compute_ldscore_grid_squared_prior_variance(borzoi_vec, distribution_obj)
	if np.any(per_snp_var < 0.0):
		print('warning: ldscore_grid_squared predicted negative prior variances; flooring at 0.0 for sampling')
	per_snp_sd = np.sqrt(np.maximum(per_snp_var, 0.0))
	sampled_betas = np.random.normal(loc=per_snp_mu[None, :], scale=per_snp_sd[None, :], size=(n_samp, len(borzoi_vec)))
	return sampled_betas

def estimate_cis_snp_heritability_with_lrt(genotype_mat, expr_vec):
	X = np.asarray(genotype_mat, dtype=float)
	y = np.asarray(expr_vec, dtype=float)
	y = y - np.mean(y)

	X = X - np.mean(X, axis=0)
	snp_sdevs = np.std(X, axis=0)
	valid_snps = np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
	X = X[:, valid_snps]
	snp_sdevs = snp_sdevs[valid_snps]
	if X.shape[1] == 0:
		return np.nan, np.nan

	X = X/snp_sdevs[None, :]
	n_samples = X.shape[0]
	grm = np.dot(X, np.transpose(X))/X.shape[1]
	eigenvalues, eigenvectors = np.linalg.eigh(grm)
	transformed_y = np.dot(np.transpose(eigenvectors), y)

	def log_likelihood(h2):
		variance_scale = h2*eigenvalues + (1.0 - h2)
		if np.any(variance_scale <= 0.0):
			return -np.inf
		residual_var = np.mean(np.square(transformed_y)/variance_scale)
		if residual_var <= 0.0:
			return -np.inf
		return -0.5*(n_samples*np.log(2.0*np.pi) + n_samples*np.log(residual_var) + np.sum(np.log(variance_scale)) + n_samples)

	null_log_likelihood = log_likelihood(0.0)
	opt = minimize_scalar(lambda h2: -log_likelihood(h2), bounds=(0.0, 0.999999), method='bounded')
	h2 = opt.x
	alt_log_likelihood = -opt.fun
	lrt_stat = np.maximum(2.0*(alt_log_likelihood - null_log_likelihood), 0.0)
	lrt_pvalue = 0.5*stats.chi2.sf(lrt_stat, df=1)
	return h2, lrt_pvalue

def compute_mean_causal_effects_from_distribution(borzoi_vec, anno, distribution_obj):
	if len(borzoi_vec) != anno.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	if distribution_obj.get('model_effect_scale') in ('allelic_grid', 'allelic_grid_squared'):
		return distribution_obj['a_prior']*borzoi_vec

	bin_indices = distribution_obj['bin_indices']
	snp_bins = np.argmax(anno, axis=1)
	if np.all(np.sum(anno, axis=1) == 1.0) == False:
		print('assumption eroror')
		pdb.set_trace()

	if 'avg_mu' in distribution_obj:
		if anno.shape[1] != len(distribution_obj['avg_mu']):
			print('assumption eroror')
			pdb.set_trace()
		if np.array_equal(bin_indices, np.arange(len(distribution_obj['avg_mu']))) == False:
			print('assumption eroror')
			pdb.set_trace()
		return distribution_obj['avg_mu'][snp_bins]*borzoi_vec

	if 'effect_size_grid_pis' in distribution_obj:
		effect_size_grid_pis = distribution_obj['effect_size_grid_pis']
		if anno.shape[1] != effect_size_grid_pis.shape[0]:
			print('assumption eroror')
			pdb.set_trace()
		if np.array_equal(bin_indices, np.arange(effect_size_grid_pis.shape[0])) == False:
			print('assumption eroror')
			pdb.set_trace()
		if 'effect_size_grid_lower' in distribution_obj and 'effect_size_grid_upper' in distribution_obj:
			per_grid_mean = (distribution_obj['effect_size_grid_lower'] + distribution_obj['effect_size_grid_upper'])/2.0
		else:
			per_grid_mean = distribution_obj['fixed_effect_size_grid']
		per_bin_mean = np.dot(effect_size_grid_pis, per_grid_mean)
		return per_bin_mean[snp_bins]

	print('assumption eroror')
	pdb.set_trace()

def sample_causal_effects_from_ashr_distribution(borzoi_vec, anno, distribution_obj, n_samp=100):
	if len(borzoi_vec) != anno.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	avg_mu = distribution_obj['avg_mu']
	prior_var_pis = distribution_obj['prior_var_pis']
	prior_var_fixed_grid = distribution_obj['prior_var_fixed_grid']
	bin_indices = distribution_obj['bin_indices']

	if anno.shape[1] != len(avg_mu):
		print('assumption eroror')
		pdb.set_trace()
	if anno.shape[1] != prior_var_pis.shape[0]:
		print('assumption eroror')
		pdb.set_trace()
	if np.array_equal(bin_indices, np.arange(len(avg_mu))) == False:
		print('assumption eroror')
		pdb.set_trace()

	snp_bins = np.argmax(anno, axis=1)
	if np.all(np.sum(anno, axis=1) == 1.0) == False:
		print('assumption eroror')
		pdb.set_trace()

	n_snps = len(borzoi_vec)
	sampled_betas = np.zeros((n_samp, n_snps))
	for snp_iter, snp_bin in enumerate(snp_bins):
		per_snp_mu = avg_mu[snp_bin]*borzoi_vec[snp_iter]
		latent_states = np.random.choice(len(prior_var_fixed_grid), size=n_samp, p=prior_var_pis[snp_bin, :])
		per_snp_sd = np.sqrt(prior_var_fixed_grid[latent_states])
		sampled_betas[:, snp_iter] = np.random.normal(loc=per_snp_mu, scale=per_snp_sd)
	return sampled_betas

def sample_causal_effects_from_effect_size_grid_distribution(borzoi_vec, anno, distribution_obj, n_samp=100):
	if len(borzoi_vec) != anno.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	effect_size_grid_pis = distribution_obj['effect_size_grid_pis']
	fixed_effect_size_grid = distribution_obj['fixed_effect_size_grid']
	bin_indices = distribution_obj['bin_indices']

	if anno.shape[1] != effect_size_grid_pis.shape[0]:
		print('assumption eroror')
		pdb.set_trace()
	if np.array_equal(bin_indices, np.arange(effect_size_grid_pis.shape[0])) == False:
		print('assumption eroror')
		pdb.set_trace()

	snp_bins = np.argmax(anno, axis=1)
	if np.all(np.sum(anno, axis=1) == 1.0) == False:
		print('assumption eroror')
		pdb.set_trace()

	n_snps = len(borzoi_vec)
	sampled_betas = np.zeros((n_samp, n_snps))
	for snp_iter, snp_bin in enumerate(snp_bins):
		latent_states = np.random.choice(len(fixed_effect_size_grid), size=n_samp, p=effect_size_grid_pis[snp_bin, :])
		sampled_betas[:, snp_iter] = fixed_effect_size_grid[latent_states]
	return sampled_betas


def sample_causal_effects_from_uniform_prior_grid_distribution(borzoi_vec, anno, distribution_obj, n_samp=100):
	if len(borzoi_vec) != anno.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	effect_size_grid_pis = distribution_obj['effect_size_grid_pis']
	effect_size_grid_lower = distribution_obj['effect_size_grid_lower']
	effect_size_grid_upper = distribution_obj['effect_size_grid_upper']
	bin_indices = distribution_obj['bin_indices']

	if anno.shape[1] != effect_size_grid_pis.shape[0]:
		print('assumption eroror')
		pdb.set_trace()
	if np.array_equal(bin_indices, np.arange(effect_size_grid_pis.shape[0])) == False:
		print('assumption eroror')
		pdb.set_trace()

	snp_bins = np.argmax(anno, axis=1)
	if np.all(np.sum(anno, axis=1) == 1.0) == False:
		print('assumption eroror')
		pdb.set_trace()

	n_snps = len(borzoi_vec)
	sampled_betas = np.zeros((n_samp, n_snps))
	for snp_iter, snp_bin in enumerate(snp_bins):
		latent_states = np.random.choice(len(effect_size_grid_lower), size=n_samp, p=effect_size_grid_pis[snp_bin, :])
		sampled_betas[:, snp_iter] = np.random.uniform(low=effect_size_grid_lower[latent_states], high=effect_size_grid_upper[latent_states])
	return sampled_betas



def run_expression_prediction(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, distribution_obj, distribution, output_file, standardize_geno):
	# Initialize output file
	t = open(output_file,'w')
	t.write('gene_id\traw_correlation\trecalibrated_correlation\tmax_abs_borzoi\tmax_abs_standardized_borzoi\tself_corr_positive_fraction\tcorr_distribution_positive_prob\tSNR\tavg_signal\tdirectional_FSR\tdirectional_FSR_v2\traw_correlation_pvalue\traw_correlation_prob_negative\tcis_snp_h2\tcis_snp_h2_pvalue\n')


	# Loop through chromsomes
	for chrom_num in range(1,23):
		print(chrom_num)

		##################################
		# Load in per-chrom-genotype data
		##################################
		# string of chromosome name
		chrom_string = 'chr' + str(chrom_num)
		# Load in chromosome plink data
		(bim, fam, G) = read_plink(genotype_stem + str(chrom_num))
		# Create mapping from variant id to index
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		# Create mapping from rsid to a0, a1
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))


		##################################
		# Loop through genes on this chromosome
		# (Analysis done seperately for each gene)
		##################################
		for gene_id in [*gene_id_to_est_borzoi_effects]:

			# Limit to genes on this chromosome
			gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
			if str(gene_chrom_num) != str(chrom_num):
				continue

			# Gene needs both borzoi effects AND expression and variant gene anno
			if gene_id not in gene_id_to_expression_vector:
				continue
			if distribution not in ('ldscore_grid', 'ldscore_grid_squared') and gene_id not in gene_id_to_variant_gene_anno:
				continue

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])
			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 10:
				continue

			# Load in data for gene
			# Borzoi
			borzoi_effects_unstandardized, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])

			# Anno
			if distribution in ('ldscore_grid', 'ldscore_grid_squared'):
				variant_anno = np.zeros((len(ordered_cis_variants), 0))
				borzoi_anno_variant_alleles = borzoi_variant_alleles
			else:
				variant_anno, borzoi_anno_variant_alleles = load_in_snp_gene_anno_data(ordered_cis_variants, gene_id_to_variant_gene_anno[gene_id], len(anno_names))

			# Load in LD
			cis_genotype_indices = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				
				# Also flip signs of borzoi effects to match LD
				if np.isnan(borzoi_effects_unstandardized[var_index]) == False:
					if borzoi_variant_alleles[var_index,:][0] == geno_alleles[0]:
						borzoi_effects_unstandardized[var_index] = -1.0*borzoi_effects_unstandardized[var_index]
						pdb.set_trace()
				if distribution not in ('ldscore_grid', 'ldscore_grid_squared') and borzoi_variant_alleles[var_index,:][0] != borzoi_anno_variant_alleles[var_index,:][0]:
					print('annotation alllele assumption erororo')
					pdb.set_trace()
			
			# Extract genotype
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = (G[cis_genotype_indices,:].compute())[:, genotype_sample_indices]
			row_means = np.nanmean(geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(geno_mat))
			geno_mat[nan_rows, nan_cols] = row_means[nan_rows]

			geno_mat = 2.0 - geno_mat

			snp_means = np.mean(geno_mat, axis=1)
			snp_sdevs = np.std(geno_mat, axis=1)
			zero_var_snps = snp_sdevs == 0.0
			valid_snps = np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
			snp_sdevs[zero_var_snps] = 1.0
			if standardize_geno == 'True':
				geno_mat = (geno_mat - snp_means[:, None])/snp_sdevs[:, None]
				borzoi_effects_standardized = borzoi_effects_unstandardized*snp_sdevs
			else:
				geno_mat = (geno_mat - snp_means[:, None])
				borzoi_effects_standardized = borzoi_effects_unstandardized

			expr_vec = gene_id_to_expression_vector[gene_id]
			anno = variant_anno[valid_snps, :]
			tmp_standardized_borzoi_vec = (borzoi_effects_unstandardized*snp_sdevs)[valid_snps]
			borzoi_vec = borzoi_effects_standardized[valid_snps]
			borzoi_vec_unstandardized = borzoi_effects_unstandardized[valid_snps]
			genotype_mat = np.transpose(geno_mat[valid_snps, :])


			cis_snp_h2, cis_snp_h2_pvalue = estimate_cis_snp_heritability_with_lrt(genotype_mat, expr_vec)


			default_pred_expr = np.dot(genotype_mat, borzoi_vec)
			default_corry = np.corrcoef(expr_vec, default_pred_expr)[0,1]
			model_mean_borzois = compute_mean_causal_effects_from_distribution(borzoi_vec, anno, distribution_obj)
			recalibrated_pred_expr = np.dot(genotype_mat, model_mean_borzois)
			recalibrated_corry = np.corrcoef(expr_vec, recalibrated_pred_expr)[0,1]

			if distribution == 'gaussian':
				sampled_borzois = sample_causal_effects_from_gaussian_distribution(borzoi_vec, anno, distribution_obj, n_samp=300)
			elif distribution == 'ashr':
				sampled_borzois = sample_causal_effects_from_ashr_distribution(borzoi_vec, anno, distribution_obj, n_samp=300)
			elif distribution == 'effect_size_grid':
				sampled_borzois = sample_causal_effects_from_effect_size_grid_distribution(borzoi_vec, anno, distribution_obj, n_samp=300)
			elif distribution.startswith('uniform_prior_grid'):
				sampled_borzois = sample_causal_effects_from_uniform_prior_grid_distribution(borzoi_vec, anno, distribution_obj, n_samp=300)
			elif distribution == 'ldscore_grid':
				sampled_borzois = sample_causal_effects_from_ldscore_grid_distribution(borzoi_vec, distribution_obj, n_samp=1000)
			elif distribution == 'ldscore_grid_squared':
				sampled_borzois = sample_causal_effects_from_ldscore_grid_squared_distribution(borzoi_vec, distribution_obj, n_samp=1000)
			else:
				print('distribution assumption eroror')
				pdb.set_trace()

			sampled_pred_expr = np.dot(genotype_mat, np.transpose(sampled_borzois))
			self_corr_mat = np.corrcoef(np.transpose(sampled_pred_expr))
			self_corrs = self_corr_mat[np.triu_indices(self_corr_mat.shape[0], k=1)]
			self_corr_positive_fraction = np.sum(self_corrs > 0)/len(self_corrs)


			signal = np.var(recalibrated_pred_expr,ddof=1)      # Var_i(E_s[X_i beta^(s)])
			noise = np.mean(np.var(sampled_pred_expr, axis=1, ddof=1))  # E_i(Var_s[X_i beta^(s)])
			R_g = signal / (signal + noise + 1e-12)


			avg_signal = signal/len(model_mean_borzois)



			# Reported score: X mu
			reported_score = recalibrated_pred_expr.copy()

			# Samples of true genetic expression: X beta^(s)
			# sampled_pred_expr has shape: n_individuals x n_samples

			# Center over individuals so this corresponds to covariance/correlation direction
			reported_score_centered = reported_score - np.mean(reported_score)
			sampled_pred_expr_centered = sampled_pred_expr - np.mean(sampled_pred_expr, axis=0)

			# Alignment between reported score and each sampled true score
			# proportional to (X mu)^T (X beta_true^(s))
			alignment_samples = np.dot(reported_score_centered, sampled_pred_expr_centered)

			# Monte Carlo FSR: probability reported score points wrong way
			directional_fsr_mc = np.mean(alignment_samples < 0)

			default_corry2, default_corry_p_value = stats.pearsonr(expr_vec, default_pred_expr)
			z = np.arctanh(default_corry)
			se = 1.0/np.sqrt(len(expr_vec) - 3)
			prob_negative = stats.norm.cdf(0.0, loc=z, scale=se)


			corr_distr = []
			for col_iter in range(sampled_pred_expr.shape[1]):
				corr_distr.append(np.corrcoef(expr_vec, sampled_pred_expr[:, col_iter])[0,1])
			corr_distr = np.asarray(corr_distr)
			corr_dist_positive_prob = np.sum(corr_distr > 0)/len(corr_distr)

			finite_sample_corrs = corr_distr[np.isfinite(corr_distr)]
			if len(finite_sample_corrs) == 0 or len(expr_vec) <= 3:
				directional_fsr_mc_v2 = np.nan
			else:
				finite_sample_z = np.arctanh(np.clip(finite_sample_corrs, -0.999999, 0.999999))
				finite_sample_se = 1.0/np.sqrt(len(expr_vec) - 3)
				directional_fsr_mc_v2 = np.mean(stats.norm.cdf(0.0, loc=finite_sample_z, scale=finite_sample_se))


			t.write(gene_id + '\t' + str(default_corry) + '\t' + str(recalibrated_corry) + '\t' + str(np.max(np.abs(borzoi_vec_unstandardized))) + '\t' + str(np.max(np.abs(tmp_standardized_borzoi_vec))) + '\t' + str(self_corr_positive_fraction) + '\t' + str(corr_dist_positive_prob) + '\t' + str(R_g) + '\t' + str(avg_signal) + '\t' + str(directional_fsr_mc) + '\t' + str(directional_fsr_mc_v2) + '\t' + str(default_corry_p_value) + '\t' + str(prob_negative) + '\t' + str(cis_snp_h2) + '\t' + str(cis_snp_h2_pvalue) + '\n')
			t.flush()

	t.close()


	return 



def load_in_gaussian_distribution_data(borzoi_based_prior_output_stem):
	posterior_summary_file = borzoi_based_prior_output_stem + '_posterior_summary.txt'

	mu_mapping = {}
	var_mapping = {}

	f = open(posterior_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		param_name = data[0]
		param_mean = float(data[1])

		if param_name.startswith('mu_prior_'):
			bin_index = int(param_name.split('mu_prior_')[1])
			mu_mapping[bin_index] = param_mean
		elif param_name.startswith('sig_sq_prior_'):
			bin_index = int(param_name.split('sig_sq_prior_')[1])
			var_mapping[bin_index] = param_mean
	f.close()

	ordered_mu_bins = np.sort(np.asarray([*mu_mapping]))
	ordered_var_bins = np.sort(np.asarray([*var_mapping]))
	if np.array_equal(ordered_mu_bins, ordered_var_bins) == False:
		print('assumption erororor')
		pdb.set_trace()

	avg_mu = np.asarray([mu_mapping[bin_index] for bin_index in ordered_mu_bins])
	avg_var = np.asarray([var_mapping[bin_index] for bin_index in ordered_var_bins])

	distribution_obj = {}
	distribution_obj['avg_mu'] = avg_mu
	distribution_obj['avg_var'] = avg_var
	distribution_obj['bin_indices'] = ordered_mu_bins
	return distribution_obj

def load_in_ashr_distribution_data(borzoi_based_prior_output_stem):
	posterior_summary_file = borzoi_based_prior_output_stem + '_posterior_summary.txt'
	prior_variance_grid_file = borzoi_based_prior_output_stem + '_prior_variance_grid.txt'

	mu_mapping = {}
	pi_mapping = {}

	f = open(posterior_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		param_name = data[0]
		param_mean = float(data[1])

		if param_name.startswith('mu_prior_'):
			bin_index = int(param_name.split('mu_prior_')[1])
			mu_mapping[bin_index] = param_mean
		elif param_name.startswith('prior_var_pi_'):
			param_suffix = param_name.split('prior_var_pi_')[1]
			bin_index = int(param_suffix.split('_')[0])
			grid_index = int(param_suffix.split('_')[1])
			if bin_index not in pi_mapping:
				pi_mapping[bin_index] = {}
			pi_mapping[bin_index][grid_index] = param_mean
	f.close()

	grid_vals = []
	f = open(prior_variance_grid_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		grid_vals.append(float(data[1]))
	f.close()
	prior_var_fixed_grid = np.asarray(grid_vals)

	ordered_mu_bins = np.sort(np.asarray([*mu_mapping]))
	ordered_pi_bins = np.sort(np.asarray([*pi_mapping]))
	if np.array_equal(ordered_mu_bins, ordered_pi_bins) == False:
		print('assumption erororor')
		pdb.set_trace()

	avg_mu = np.asarray([mu_mapping[bin_index] for bin_index in ordered_mu_bins])
	prior_var_pis = np.zeros((len(ordered_mu_bins), len(prior_var_fixed_grid)))
	for row_iter, bin_index in enumerate(ordered_mu_bins):
		ordered_grid_indices = np.sort(np.asarray([*pi_mapping[bin_index]]))
		if np.array_equal(ordered_grid_indices, np.arange(len(prior_var_fixed_grid))) == False:
			print('assumption erororor')
			pdb.set_trace()
		prior_var_pis[row_iter, :] = np.asarray([pi_mapping[bin_index][grid_index] for grid_index in ordered_grid_indices])

	distribution_obj = {}
	distribution_obj['avg_mu'] = avg_mu
	distribution_obj['prior_var_pis'] = prior_var_pis
	distribution_obj['prior_var_fixed_grid'] = prior_var_fixed_grid
	distribution_obj['bin_indices'] = ordered_mu_bins
	return distribution_obj

def load_in_effect_size_grid_distribution_data(borzoi_based_prior_output_stem):
	posterior_summary_file = borzoi_based_prior_output_stem + '_posterior_summary.txt'
	fixed_effect_size_grid_file = borzoi_based_prior_output_stem + '_fixed_effect_size_grid.txt'

	pi_mapping = {}

	f = open(posterior_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		param_name = data[0]
		param_mean = float(data[1])

		if param_name.startswith('effect_size_grid_pi_'):
			param_suffix = param_name.split('effect_size_grid_pi_')[1]
			bin_index = int(param_suffix.split('_')[0])
			grid_index = int(param_suffix.split('_')[1])
			if bin_index not in pi_mapping:
				pi_mapping[bin_index] = {}
			pi_mapping[bin_index][grid_index] = param_mean
	f.close()

	grid_vals = []
	f = open(fixed_effect_size_grid_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		grid_vals.append(float(data[1]))
	f.close()
	fixed_effect_size_grid = np.asarray(grid_vals)

	ordered_pi_bins = np.sort(np.asarray([*pi_mapping]))
	effect_size_grid_pis = np.zeros((len(ordered_pi_bins), len(fixed_effect_size_grid)))
	for row_iter, bin_index in enumerate(ordered_pi_bins):
		ordered_grid_indices = np.sort(np.asarray([*pi_mapping[bin_index]]))
		if np.array_equal(ordered_grid_indices, np.arange(len(fixed_effect_size_grid))) == False:
			print('assumption erororor')
			pdb.set_trace()
		effect_size_grid_pis[row_iter, :] = np.asarray([pi_mapping[bin_index][grid_index] for grid_index in ordered_grid_indices])

	distribution_obj = {}
	distribution_obj['effect_size_grid_pis'] = effect_size_grid_pis
	distribution_obj['fixed_effect_size_grid'] = fixed_effect_size_grid
	distribution_obj['bin_indices'] = ordered_pi_bins
	return distribution_obj


def load_in_uniform_prior_grid_distribution_data(borzoi_based_prior_output_stem):
	posterior_summary_file = borzoi_based_prior_output_stem + '_posterior_summary.txt'
	fixed_effect_size_grid_file = borzoi_based_prior_output_stem + '_fixed_effect_size_grid.txt'

	pi_mapping = {}

	f = open(posterior_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		param_name = data[0]
		param_mean = float(data[1])

		if param_name.startswith('effect_size_grid_pi_'):
			param_suffix = param_name.split('effect_size_grid_pi_')[1]
			bin_index = int(param_suffix.split('_')[0])
			grid_index = int(param_suffix.split('_')[1])
			if bin_index not in pi_mapping:
				pi_mapping[bin_index] = {}
			pi_mapping[bin_index][grid_index] = param_mean
	f.close()

	grid_vals = []
	grid_lower = []
	grid_upper = []
	f = open(fixed_effect_size_grid_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		grid_vals.append(float(data[1]))
		grid_lower.append(float(data[2]))
		grid_upper.append(float(data[3]))
	f.close()
	fixed_effect_size_grid = np.asarray(grid_vals)
	effect_size_grid_lower = np.asarray(grid_lower)
	effect_size_grid_upper = np.asarray(grid_upper)

	ordered_pi_bins = np.sort(np.asarray([*pi_mapping]))
	effect_size_grid_pis = np.zeros((len(ordered_pi_bins), len(fixed_effect_size_grid)))
	for row_iter, bin_index in enumerate(ordered_pi_bins):
		ordered_grid_indices = np.sort(np.asarray([*pi_mapping[bin_index]]))
		if np.array_equal(ordered_grid_indices, np.arange(len(fixed_effect_size_grid))) == False:
			print('assumption erororor')
			pdb.set_trace()
		effect_size_grid_pis[row_iter, :] = np.asarray([pi_mapping[bin_index][grid_index] for grid_index in ordered_grid_indices])

	distribution_obj = {}
	distribution_obj['effect_size_grid_pis'] = effect_size_grid_pis
	distribution_obj['fixed_effect_size_grid'] = fixed_effect_size_grid
	distribution_obj['effect_size_grid_lower'] = effect_size_grid_lower
	distribution_obj['effect_size_grid_upper'] = effect_size_grid_upper
	distribution_obj['bin_indices'] = ordered_pi_bins
	return distribution_obj


def load_in_ldscore_grid_distribution_data_from_files(borzoi_based_prior_output_stem, distribution_name, expected_model_effect_scale):
	summary_file = borzoi_based_prior_output_stem + '_' + distribution_name + '_summary.txt'
	param_file = borzoi_based_prior_output_stem + '_' + distribution_name + '_params.txt'

	grid_knots = None
	model_effect_scale = None
	f = open(summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'grid_knots':
			grid_knots = np.asarray(data[1].split(',')).astype(float)
		elif data[0] == 'model_effect_scale':
			model_effect_scale = data[1]
	f.close()

	if grid_knots is None:
		print(distribution_name + ' prior loading error: could not find grid_knots in ' + summary_file)
		pdb.set_trace()
	if model_effect_scale != expected_model_effect_scale:
		print(distribution_name + ' prior loading error: expected model_effect_scale=' + expected_model_effect_scale)
		pdb.set_trace()

	a_prior = None
	grid_coef_mapping = {}
	f = open(param_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		param_name = data[0]
		param_index = int(data[1])
		param_val = float(data[2])
		if param_name == 'a_prior':
			a_prior = param_val
		elif param_name == 'grid_coef':
			grid_coef_mapping[param_index] = param_val
	f.close()

	if a_prior is None:
		print(distribution_name + ' prior loading error: could not find a_prior in ' + param_file)
		pdb.set_trace()

	ordered_coef_indices = np.sort(np.asarray([*grid_coef_mapping]))
	if np.array_equal(ordered_coef_indices, np.arange(len(grid_knots) + 2)) == False:
		print(distribution_name + ' prior loading error: grid coefficient indices do not match expected basis size')
		pdb.set_trace()
	grid_coefs = np.asarray([grid_coef_mapping[index] for index in ordered_coef_indices])

	distribution_obj = {}
	distribution_obj['model_effect_scale'] = model_effect_scale
	distribution_obj['a_prior'] = a_prior
	distribution_obj['grid_knots'] = grid_knots
	distribution_obj['grid_coefs'] = grid_coefs
	return distribution_obj


def load_in_ldscore_grid_distribution_data(borzoi_based_prior_output_stem):
	return load_in_ldscore_grid_distribution_data_from_files(borzoi_based_prior_output_stem, 'ldscore_grid', 'allelic_grid')

def load_in_ldscore_grid_squared_distribution_data(borzoi_based_prior_output_stem):
	return load_in_ldscore_grid_distribution_data_from_files(borzoi_based_prior_output_stem, 'ldscore_grid_squared', 'allelic_grid_squared')


########################
# Command line args
########################
borzoi_effect_file = sys.argv[1]
borzoi_annotation_file = sys.argv[2]
genotype_stem = sys.argv[3]
genotype_sample_mapping_file = sys.argv[4]
expr_file = sys.argv[5]
distribution = sys.argv[6]
borzoi_based_prior_output_stem = sys.argv[7]
output_file = sys.argv[8]
standardize_geno = sys.argv[9]


###########################
# Load in data
###########################

if distribution == 'gaussian':
	distribution_obj = load_in_gaussian_distribution_data(borzoi_based_prior_output_stem)
elif distribution == 'ashr':
	distribution_obj = load_in_ashr_distribution_data(borzoi_based_prior_output_stem)
elif distribution == 'effect_size_grid':
	distribution_obj = load_in_effect_size_grid_distribution_data(borzoi_based_prior_output_stem)
elif distribution.startswith('uniform_prior_grid'):
	distribution_obj = load_in_uniform_prior_grid_distribution_data(borzoi_based_prior_output_stem)
elif distribution == 'ldscore_grid':
	if standardize_geno != 'False':
		print('ldscore_grid prediction is currently implemented only for standardize_geno=False because the prior is on the allelic effect scale')
		sys.exit(1)
	distribution_obj = load_in_ldscore_grid_distribution_data(borzoi_based_prior_output_stem)
elif distribution == 'ldscore_grid_squared':
	if standardize_geno != 'False':
		print('ldscore_grid_squared prediction is currently implemented only for standardize_geno=False because the prior is on the allelic effect scale')
		sys.exit(1)
	distribution_obj = load_in_ldscore_grid_squared_distribution_data(borzoi_based_prior_output_stem)
else:
	print('distribution assumption eroror')
	pdb.set_trace()

# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)

# Create mapping from gene id to vector of variant-gene annotations
if distribution in ('ldscore_grid', 'ldscore_grid_squared'):
	gene_id_to_variant_gene_anno = None
	anno_names = np.asarray([])
else:
	gene_id_to_variant_gene_anno, anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(borzoi_annotation_file)

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

# Create mapping from gene id to expression vector
gene_id_to_expression_vector = create_mapping_from_gene_id_to_expression_vector(expr_file)


# run expression prediction
run_expression_prediction(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, distribution_obj, distribution, output_file, standardize_geno)
