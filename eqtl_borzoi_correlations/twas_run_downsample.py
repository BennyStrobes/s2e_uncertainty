import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time
from scipy.optimize import minimize_scalar
from scipy import stats
from scipy.special import logsumexp

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

def compute_ldscore_grid_squared_prior_variance(borzoi_vec, distribution_obj):
	grid_basis = create_continuous_piecewise_linear_squared_basis(borzoi_vec, distribution_obj['grid_knots'])
	return np.dot(grid_basis, distribution_obj['grid_coefs'])

def compute_ldscore_grid_squared_prior_distribution(borzoi_vec, distribution_obj):
	prior_mean_betas = distribution_obj['a_prior']*borzoi_vec
	prior_var_betas = compute_ldscore_grid_squared_prior_variance(borzoi_vec, distribution_obj)
	if np.any(prior_var_betas < 0.0):
		print('warning: ldscore_grid_squared predicted negative prior variances; flooring at 0.0 for prior covariance')
	prior_var_betas = np.maximum(prior_var_betas, 0.0)
	return prior_mean_betas, np.diag(prior_var_betas)

def sample_causal_effects_from_ldscore_grid_squared_distribution(borzoi_vec, distribution_obj, n_samp=100):
	per_snp_mu = distribution_obj['a_prior']*borzoi_vec
	per_snp_var = compute_ldscore_grid_squared_prior_variance(borzoi_vec, distribution_obj)
	if np.any(per_snp_var < 0.0):
		print('warning: ldscore_grid_squared predicted negative prior variances; flooring at 0.0 for sampling')
	per_snp_sd = np.sqrt(np.maximum(per_snp_var, 0.0))
	sampled_betas = np.random.normal(loc=per_snp_mu[None, :], scale=per_snp_sd[None, :], size=(n_samp, len(borzoi_vec)))
	return sampled_betas

def compute_ldscore_grid_squared_posterior_distribution(training_genotype_mat, training_expr_vec, borzoi_vec, distribution_obj, full_covariance_snp_limit=2000):
	prior_mean_betas = distribution_obj['a_prior']*borzoi_vec
	prior_var_betas = compute_ldscore_grid_squared_prior_variance(borzoi_vec, distribution_obj)
	if np.any(prior_var_betas < 0.0):
		print('warning: ldscore_grid_squared predicted negative prior variances; flooring at 0.0 for posterior update')
	prior_var_betas = np.maximum(prior_var_betas, 1e-12)

	training_expr_vec_centered = training_expr_vec - np.mean(training_expr_vec)
	training_prior_pred_expr = np.dot(training_genotype_mat, prior_mean_betas)
	training_residual_var = np.var(training_expr_vec_centered - training_prior_pred_expr, ddof=1)
	training_residual_var = np.maximum(training_residual_var, 1e-12)

	training_genotype_prior_scaled = training_genotype_mat*prior_var_betas[None, :]
	posterior_sample_covariance = np.dot(training_genotype_prior_scaled, np.transpose(training_genotype_mat)) + training_residual_var*np.eye(training_genotype_mat.shape[0])
	posterior_alpha = np.linalg.solve(posterior_sample_covariance, training_expr_vec_centered - training_prior_pred_expr)
	posterior_mean_betas = prior_mean_betas + prior_var_betas*np.dot(np.transpose(training_genotype_mat), posterior_alpha)

	posterior_sample_covariance_inv_training_genotype = np.linalg.solve(posterior_sample_covariance, training_genotype_mat)
	posterior_var_betas = prior_var_betas - np.square(prior_var_betas)*np.sum(training_genotype_mat*posterior_sample_covariance_inv_training_genotype, axis=0)
	posterior_var_betas = np.maximum(posterior_var_betas, 0.0)

	if len(borzoi_vec) > full_covariance_snp_limit:
		posterior_covariance_mat = np.diag(posterior_var_betas)
	else:
		posterior_covariance_mat = np.diag(prior_var_betas) - prior_var_betas[:, None]*np.dot(np.transpose(training_genotype_mat), posterior_sample_covariance_inv_training_genotype)*prior_var_betas[None, :]
		posterior_covariance_mat[np.diag_indices_from(posterior_covariance_mat)] = posterior_var_betas

	return posterior_mean_betas, posterior_covariance_mat

def sample_effects_from_gaussian_distribution(effect_mean, effect_covariance_mat, n_samp=300):
	effect_covariance_mat = (effect_covariance_mat + np.transpose(effect_covariance_mat))/2.0
	effect_sdev = np.sqrt(np.maximum(np.diag(effect_covariance_mat), 0.0))
	if np.all(np.isfinite(effect_mean)) == False or np.all(np.isfinite(effect_sdev)) == False or np.all(np.isfinite(effect_covariance_mat)) == False:
		print('warning: non-finite effect covariance; using diagonal covariance only')
		return np.random.normal(loc=effect_mean[None, :], scale=effect_sdev[None, :], size=(n_samp, len(effect_mean)))
	if len(effect_mean) > 2000:
		return np.random.normal(loc=effect_mean[None, :], scale=effect_sdev[None, :], size=(n_samp, len(effect_mean)))
	try:
		return np.random.multivariate_normal(mean=effect_mean, cov=effect_covariance_mat, size=n_samp)
	except np.linalg.LinAlgError:
		print('warning: effect covariance SVD failed; using diagonal covariance only')
		return np.random.normal(loc=effect_mean[None, :], scale=effect_sdev[None, :], size=(n_samp, len(effect_mean)))

def compute_prediction_correlation_and_directional_fsr(validation_genotype_mat, validation_expr_vec, effect_mean, effect_covariance_mat, n_samp=300):
	predicted_expr = np.dot(validation_genotype_mat, effect_mean)
	prediction_correlation = np.corrcoef(validation_expr_vec, predicted_expr)[0,1]
	_, correlation_pvalue = stats.pearsonr(validation_expr_vec, predicted_expr)
	fisher_z = np.arctanh(prediction_correlation)
	fisher_z_se = 1.0/np.sqrt(len(validation_expr_vec) - 3)
	correlation_prob_negative = stats.norm.cdf(0.0, loc=fisher_z, scale=fisher_z_se)

	sampled_effects = sample_effects_from_gaussian_distribution(effect_mean, effect_covariance_mat, n_samp=n_samp)
	sampled_pred_expr = np.dot(validation_genotype_mat, np.transpose(sampled_effects))

	predicted_expr_centered = predicted_expr - np.mean(predicted_expr)
	sampled_pred_expr_centered = sampled_pred_expr - np.mean(sampled_pred_expr, axis=0)
	alignment_samples = np.dot(predicted_expr_centered, sampled_pred_expr_centered)
	directional_fsr_mc = np.mean(alignment_samples < 0)

	return prediction_correlation, directional_fsr_mc, correlation_pvalue, correlation_prob_negative

def compute_directional_fsr_from_genotype(genotype_mat, effect_mean, effect_covariance_mat, n_samp=300):
	predicted_expr = np.dot(genotype_mat, effect_mean)
	sampled_effects = sample_effects_from_gaussian_distribution(effect_mean, effect_covariance_mat, n_samp=n_samp)
	sampled_pred_expr = np.dot(genotype_mat, np.transpose(sampled_effects))

	predicted_expr_centered = predicted_expr - np.mean(predicted_expr)
	sampled_pred_expr_centered = sampled_pred_expr - np.mean(sampled_pred_expr, axis=0)
	alignment_samples = np.dot(predicted_expr_centered, sampled_pred_expr_centered)
	directional_fsr_mc = np.mean(alignment_samples < 0)
	return directional_fsr_mc

def ld_impute_missing_gwas_z_scores(gwas_z_scores, ld_mat, ridge=0.1, min_observed_snps=1, min_imputation_r2=0.0):
	gwas_z_scores = np.asarray(gwas_z_scores, dtype=float)
	ld_mat = np.asarray(ld_mat, dtype=float)
	observed_snps = np.isfinite(gwas_z_scores)
	missing_snps = observed_snps == False
	imputed_snps = np.zeros(len(gwas_z_scores), dtype=bool)
	imputation_r2 = np.full(len(gwas_z_scores), np.nan)
	if np.sum(missing_snps) == 0 or np.sum(observed_snps) < min_observed_snps:
		return gwas_z_scores, imputed_snps, imputation_r2

	ld_mat = (ld_mat + np.transpose(ld_mat))/2.0
	ld_mat[np.diag_indices_from(ld_mat)] = 1.0
	R_oo = ld_mat[np.ix_(observed_snps, observed_snps)]
	R_mo = ld_mat[np.ix_(missing_snps, observed_snps)]
	R_om = np.transpose(R_mo)
	R_oo = R_oo + ridge*np.eye(R_oo.shape[0])
	try:
		imputation_weights = np.linalg.solve(R_oo, R_om)
		imputation_alpha = np.linalg.solve(R_oo, gwas_z_scores[observed_snps])
	except np.linalg.LinAlgError:
		imputation_weights = np.linalg.lstsq(R_oo, R_om, rcond=None)[0]
		imputation_alpha = np.linalg.lstsq(R_oo, gwas_z_scores[observed_snps], rcond=None)[0]
	missing_z_scores = np.dot(R_mo, imputation_alpha)
	missing_imputation_r2 = np.sum(R_mo*np.transpose(imputation_weights), axis=1)
	missing_imputation_r2 = np.clip(missing_imputation_r2, 0.0, 1.0)
	pass_r2 = missing_imputation_r2 >= min_imputation_r2

	imputed_gwas_z_scores = np.copy(gwas_z_scores)
	missing_indices = np.where(missing_snps)[0]
	imputed_gwas_z_scores[missing_indices[pass_r2]] = missing_z_scores[pass_r2]
	imputed_snps[missing_indices[pass_r2]] = True
	imputation_r2[missing_indices] = missing_imputation_r2
	return imputed_gwas_z_scores, imputed_snps, imputation_r2

def run_twas_from_per_allele_sumstats(gwas_effect_sizes, gwas_effect_ses, eqtl_effect_sizes, ld_mat, snp_sdevs, impute_missing_gwas_z_scores=True, gwas_z_imputation_ridge=0.1):
	"""Run TWAS from aligned per-allele GWAS effects, eQTL weights, LD, and SNP SDs.

	Assumes gwas_effect_sizes and eqtl_effect_sizes are on the same allele dosage scale.
	Inputs are converted to standardized genotype scale before computing the TWAS statistic.
	ld_mat is the SNP-by-SNP genotype correlation matrix, and snp_sdevs are allele dosage
	standard deviations from the LD/reference panel. If impute_missing_gwas_z_scores=True,
	missing GWAS Z-scores are LD-imputed and the returned effect is the TWAS Z-score with SE=1.
	"""
	gwas_effect_sizes = np.asarray(gwas_effect_sizes, dtype=float)
	gwas_effect_ses = np.asarray(gwas_effect_ses, dtype=float)
	eqtl_effect_sizes = np.asarray(eqtl_effect_sizes, dtype=float)
	ld_mat = np.asarray(ld_mat, dtype=float)
	snp_sdevs = np.asarray(snp_sdevs, dtype=float)

	n_snps = len(eqtl_effect_sizes)
	if len(gwas_effect_sizes) != n_snps or len(gwas_effect_ses) != n_snps or len(snp_sdevs) != n_snps:
		print('twas input dimension assumption error')
		pdb.set_trace()
	if ld_mat.shape != (n_snps, n_snps):
		print('twas LD dimension assumption error')
		pdb.set_trace()

	valid_snps = np.isfinite(eqtl_effect_sizes) & np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
	if np.sum(valid_snps) == 0:
		return np.nan, np.nan, np.nan

	gwas_effect_sizes = gwas_effect_sizes[valid_snps]
	gwas_effect_ses = gwas_effect_ses[valid_snps]
	eqtl_effect_sizes = eqtl_effect_sizes[valid_snps]
	snp_sdevs = snp_sdevs[valid_snps]
	ld_mat = ld_mat[np.ix_(valid_snps, valid_snps)]
	if np.all(np.isfinite(ld_mat)) == False:
		return np.nan, np.nan, np.nan
	ld_mat = (ld_mat + np.transpose(ld_mat))/2.0
	ld_mat[np.diag_indices_from(ld_mat)] = 1.0

	eqtl_effect_sizes_standardized = eqtl_effect_sizes*snp_sdevs

	if impute_missing_gwas_z_scores:
		observed_gwas_snps = np.isfinite(gwas_effect_sizes) & np.isfinite(gwas_effect_ses) & (gwas_effect_ses > 0.0)
		if np.sum(observed_gwas_snps) == 0:
			return np.nan, np.nan, np.nan
		gwas_z_scores = np.full(len(gwas_effect_sizes), np.nan)
		gwas_z_scores[observed_gwas_snps] = gwas_effect_sizes[observed_gwas_snps]/gwas_effect_ses[observed_gwas_snps]
		gwas_z_scores, _, _ = ld_impute_missing_gwas_z_scores(gwas_z_scores, ld_mat, ridge=gwas_z_imputation_ridge)
		valid_z_snps = np.isfinite(gwas_z_scores)
		if np.sum(valid_z_snps) == 0:
			return np.nan, np.nan, np.nan
		gwas_z_scores = gwas_z_scores[valid_z_snps]
		eqtl_effect_sizes_standardized = eqtl_effect_sizes_standardized[valid_z_snps]
		ld_mat = ld_mat[np.ix_(valid_z_snps, valid_z_snps)]
		predicted_expression_variance = np.dot(eqtl_effect_sizes_standardized, np.dot(ld_mat, eqtl_effect_sizes_standardized))
		if np.isfinite(predicted_expression_variance) == False or predicted_expression_variance <= 0.0:
			return np.nan, np.nan, np.nan
		twas_z = np.dot(eqtl_effect_sizes_standardized, gwas_z_scores)/np.sqrt(predicted_expression_variance)
		twas_pvalue = 2.0*stats.norm.sf(np.abs(twas_z))
		return twas_z, 1.0, twas_pvalue

	valid_gwas_snps = np.isfinite(gwas_effect_sizes) & np.isfinite(gwas_effect_ses) & (gwas_effect_ses > 0.0)
	if np.sum(valid_gwas_snps) == 0:
		return np.nan, np.nan, np.nan
	gwas_effect_sizes = gwas_effect_sizes[valid_gwas_snps]
	gwas_effect_ses = gwas_effect_ses[valid_gwas_snps]
	snp_sdevs = snp_sdevs[valid_gwas_snps]
	eqtl_effect_sizes_standardized = eqtl_effect_sizes_standardized[valid_gwas_snps]
	ld_mat = ld_mat[np.ix_(valid_gwas_snps, valid_gwas_snps)]

	gwas_effect_sizes_standardized = gwas_effect_sizes*snp_sdevs
	gwas_effect_ses_standardized = gwas_effect_ses*snp_sdevs

	predicted_expression_variance = np.dot(eqtl_effect_sizes_standardized, np.dot(ld_mat, eqtl_effect_sizes_standardized))
	if np.isfinite(predicted_expression_variance) == False or predicted_expression_variance <= 0.0:
		return np.nan, np.nan, np.nan

	twas_numerator = np.dot(eqtl_effect_sizes_standardized, gwas_effect_sizes_standardized)
	gwas_effect_covariance_mat = ld_mat*(gwas_effect_ses_standardized[:, None]*gwas_effect_ses_standardized[None, :])
	twas_numerator_variance = np.dot(eqtl_effect_sizes_standardized, np.dot(gwas_effect_covariance_mat, eqtl_effect_sizes_standardized))
	if np.isfinite(twas_numerator_variance) == False or twas_numerator_variance <= 0.0:
		return np.nan, np.nan, np.nan

	twas_effect = twas_numerator/predicted_expression_variance
	twas_effect_se = np.sqrt(twas_numerator_variance)/predicted_expression_variance
	twas_z = twas_effect/twas_effect_se
	twas_pvalue = 2.0*stats.norm.sf(np.abs(twas_z))
	return twas_effect, twas_effect_se, twas_pvalue

def compute_twas_positive_probability_from_effect_distribution(gwas_effect_sizes, gwas_effect_ses, effect_mean, effect_covariance_mat, ld_mat, snp_sdevs, n_samp=200, gwas_z_imputation_ridge=0.1):
	"""Average P(TWAS effect > 0) across sampled eQTL effect vectors."""
	gwas_effect_sizes = np.asarray(gwas_effect_sizes, dtype=float)
	gwas_effect_ses = np.asarray(gwas_effect_ses, dtype=float)
	effect_mean = np.asarray(effect_mean, dtype=float)
	effect_covariance_mat = np.asarray(effect_covariance_mat, dtype=float)
	ld_mat = np.asarray(ld_mat, dtype=float)
	snp_sdevs = np.asarray(snp_sdevs, dtype=float)

	n_snps = len(effect_mean)
	if len(gwas_effect_sizes) != n_snps or len(gwas_effect_ses) != n_snps or len(snp_sdevs) != n_snps:
		print('twas positive probability input dimension assumption error')
		pdb.set_trace()
	if ld_mat.shape != (n_snps, n_snps) or effect_covariance_mat.shape != (n_snps, n_snps):
		print('twas positive probability matrix dimension assumption error')
		pdb.set_trace()

	valid_snps = np.isfinite(effect_mean) & np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
	if np.sum(valid_snps) == 0:
		return np.nan

	gwas_effect_sizes = gwas_effect_sizes[valid_snps]
	gwas_effect_ses = gwas_effect_ses[valid_snps]
	effect_mean = effect_mean[valid_snps]
	effect_covariance_mat = effect_covariance_mat[np.ix_(valid_snps, valid_snps)]
	ld_mat = ld_mat[np.ix_(valid_snps, valid_snps)]
	snp_sdevs = snp_sdevs[valid_snps]
	if np.all(np.isfinite(ld_mat)) == False:
		return np.nan
	ld_mat = (ld_mat + np.transpose(ld_mat))/2.0
	ld_mat[np.diag_indices_from(ld_mat)] = 1.0

	observed_gwas_snps = np.isfinite(gwas_effect_sizes) & np.isfinite(gwas_effect_ses) & (gwas_effect_ses > 0.0)
	if np.sum(observed_gwas_snps) == 0:
		return np.nan
	gwas_z_scores = np.full(len(gwas_effect_sizes), np.nan)
	gwas_z_scores[observed_gwas_snps] = gwas_effect_sizes[observed_gwas_snps]/gwas_effect_ses[observed_gwas_snps]
	gwas_z_scores, _, _ = ld_impute_missing_gwas_z_scores(gwas_z_scores, ld_mat, ridge=gwas_z_imputation_ridge)

	valid_z_snps = np.isfinite(gwas_z_scores)
	if np.sum(valid_z_snps) == 0:
		return np.nan
	gwas_z_scores = gwas_z_scores[valid_z_snps]
	effect_mean = effect_mean[valid_z_snps]
	effect_covariance_mat = effect_covariance_mat[np.ix_(valid_z_snps, valid_z_snps)]
	ld_mat = ld_mat[np.ix_(valid_z_snps, valid_z_snps)]
	snp_sdevs = snp_sdevs[valid_z_snps]

	sampled_effects = sample_effects_from_gaussian_distribution(effect_mean, effect_covariance_mat, n_samp=n_samp)
	sampled_effects_standardized = sampled_effects*snp_sdevs[None, :]
	predicted_expression_variances = np.sum(sampled_effects_standardized*np.dot(sampled_effects_standardized, ld_mat), axis=1)
	valid_samples = np.isfinite(predicted_expression_variances) & (predicted_expression_variances > 0.0)
	if np.sum(valid_samples) == 0:
		return np.nan
	twas_z_scores = np.dot(sampled_effects_standardized[valid_samples, :], gwas_z_scores)/np.sqrt(predicted_expression_variances[valid_samples])
	twas_z_scores = twas_z_scores[np.isfinite(twas_z_scores)]
	if len(twas_z_scores) == 0:
		return np.nan
	return np.mean(stats.norm.cdf(twas_z_scores))

def compute_coloc_log_abf(effect_sizes, effect_ses, prior_mean=0.0, prior_var=0.15*0.15):
	effect_sizes = np.asarray(effect_sizes, dtype=float)
	effect_vars = np.square(np.asarray(effect_ses, dtype=float))
	prior_mean = np.broadcast_to(np.asarray(prior_mean, dtype=float), effect_sizes.shape)
	prior_var = np.broadcast_to(np.asarray(prior_var, dtype=float), effect_sizes.shape)
	prior_var = np.maximum(prior_var, 0.0)

	log_abf = np.full(len(effect_sizes), -np.inf)
	valid = np.isfinite(effect_sizes) & np.isfinite(effect_vars) & np.isfinite(prior_mean) & np.isfinite(prior_var) & (effect_vars > 0.0)
	if np.any(valid) == False:
		return log_abf

	v = effect_vars[valid]
	w = prior_var[valid]
	beta = effect_sizes[valid]
	mu = prior_mean[valid]
	log_abf[valid] = 0.5*np.log(v/(v + w)) + 0.5*(np.square(beta)/v - np.square(beta - mu)/(v + w))
	return log_abf

def compute_coloc_pph_from_log_abfs(eqtl_log_abf, gwas_log_abf, p1=1e-4, p2=1e-4, p12=1e-5):
	eqtl_log_abf = np.asarray(eqtl_log_abf, dtype=float)
	gwas_log_abf = np.asarray(gwas_log_abf, dtype=float)
	valid = np.isfinite(eqtl_log_abf) & np.isfinite(gwas_log_abf)
	eqtl_log_abf = eqtl_log_abf[valid]
	gwas_log_abf = gwas_log_abf[valid]
	if len(eqtl_log_abf) == 0:
		return np.full(5, np.nan)

	log_sum_eqtl = logsumexp(eqtl_log_abf)
	log_sum_gwas = logsumexp(gwas_log_abf)
	log_sum_shared = logsumexp(eqtl_log_abf + gwas_log_abf)
	log_h3_product = log_sum_eqtl + log_sum_gwas
	if log_sum_shared >= log_h3_product:
		log_distinct = -np.inf
	else:
		log_distinct = log_h3_product + np.log1p(-np.exp(log_sum_shared - log_h3_product))

	log_hypotheses = np.asarray([
		0.0,
		np.log(p1) + log_sum_eqtl,
		np.log(p2) + log_sum_gwas,
		np.log(p1) + np.log(p2) + log_distinct,
		np.log(p12) + log_sum_shared,
	])
	return np.exp(log_hypotheses - logsumexp(log_hypotheses))

def run_default_coloc(eqtl_effect_sizes, eqtl_ses, gwas_effect_sizes, gwas_ses, eqtl_prior_var=0.15*0.15, gwas_prior_var=0.15*0.15, p1=1e-4, p2=1e-4, p12=1e-5):
	eqtl_log_abf = compute_coloc_log_abf(eqtl_effect_sizes, eqtl_ses, prior_var=eqtl_prior_var)
	gwas_log_abf = compute_coloc_log_abf(gwas_effect_sizes, gwas_ses, prior_var=gwas_prior_var)
	return compute_coloc_pph_from_log_abfs(eqtl_log_abf, gwas_log_abf, p1=p1, p2=p2, p12=p12)

def run_borzoi_infused_coloc(eqtl_effect_sizes, eqtl_ses, gwas_effect_sizes, gwas_ses, borzoi_mean_betas, borzoi_covariance_mat, gwas_prior_var=0.15*0.15, p1=1e-4, p2=1e-4, p12=1e-5):
	borzoi_prior_vars = np.diag(borzoi_covariance_mat)
	eqtl_log_abf = compute_coloc_log_abf(eqtl_effect_sizes, eqtl_ses, prior_mean=borzoi_mean_betas, prior_var=borzoi_prior_vars)
	gwas_log_abf = compute_coloc_log_abf(gwas_effect_sizes, gwas_ses, prior_var=gwas_prior_var)
	return compute_coloc_pph_from_log_abfs(eqtl_log_abf, gwas_log_abf, p1=p1, p2=p2, p12=p12)

def compute_posterior_mean_for_zero_mean_isotropic_gaussian_prior(genotype_mat, expr_vec, prior_var, residual_var):
	posterior_sample_covariance = prior_var*np.dot(genotype_mat, np.transpose(genotype_mat)) + residual_var*np.eye(genotype_mat.shape[0])
	posterior_alpha = np.linalg.solve(posterior_sample_covariance, expr_vec)
	return prior_var*np.dot(np.transpose(genotype_mat), posterior_alpha)

def learn_standardized_gaussian_prior_variance_with_cross_validation(training_genotype_mat, training_expr_vec, prior_var_grid=None, n_folds=4):
	if prior_var_grid is None:
		prior_var_grid = np.power(10.0, np.linspace(-6.0, -1.0, 8))

	n_samples = training_genotype_mat.shape[0]
	n_folds = np.min([n_folds, n_samples])
	if n_folds < 2:
		return prior_var_grid[0]

	shuffled_indices = np.random.permutation(n_samples)
	fold_indices = np.array_split(shuffled_indices, n_folds)
	cv_errors = []
	all_indices = np.arange(n_samples)
	for prior_var in prior_var_grid:
		fold_errors = []
		for validation_indices in fold_indices:
			training_indices = np.setdiff1d(all_indices, validation_indices)
			if len(training_indices) < 2 or len(validation_indices) == 0:
				continue
			training_expr_mean = np.mean(training_expr_vec[training_indices])
			fold_training_expr_vec = training_expr_vec[training_indices] - training_expr_mean
			fold_validation_expr_vec = training_expr_vec[validation_indices] - training_expr_mean
			fold_training_genotype_mat = training_genotype_mat[training_indices, :]
			fold_validation_genotype_mat = training_genotype_mat[validation_indices, :]
			residual_var = np.maximum(np.var(fold_training_expr_vec, ddof=1), 1e-12)
			posterior_mean_betas = compute_posterior_mean_for_zero_mean_isotropic_gaussian_prior(fold_training_genotype_mat, fold_training_expr_vec, prior_var, residual_var)
			predicted_expr = np.dot(fold_validation_genotype_mat, posterior_mean_betas)
			fold_errors.append(np.mean(np.square(fold_validation_expr_vec - predicted_expr)))
		if len(fold_errors) == 0:
			cv_errors.append(np.inf)
		else:
			cv_errors.append(np.mean(fold_errors))

	return prior_var_grid[np.argmin(cv_errors)]

def compute_no_borzoi_standardized_gaussian_posterior_distribution(training_genotype_mat, training_expr_vec, snp_sdevs, full_covariance_snp_limit=2000):
	standardized_training_genotype_mat = training_genotype_mat/snp_sdevs[None, :]
	prior_var = learn_standardized_gaussian_prior_variance_with_cross_validation(standardized_training_genotype_mat, training_expr_vec)

	training_expr_vec_centered = training_expr_vec - np.mean(training_expr_vec)
	residual_var = np.maximum(np.var(training_expr_vec_centered, ddof=1), 1e-12)
	posterior_mean_standardized_betas = compute_posterior_mean_for_zero_mean_isotropic_gaussian_prior(standardized_training_genotype_mat, training_expr_vec_centered, prior_var, residual_var)
	training_residual_var = np.var(training_expr_vec_centered - np.dot(standardized_training_genotype_mat, posterior_mean_standardized_betas), ddof=1)
	training_residual_var = np.maximum(training_residual_var, 1e-12)

	posterior_sample_covariance = prior_var*np.dot(standardized_training_genotype_mat, np.transpose(standardized_training_genotype_mat)) + training_residual_var*np.eye(standardized_training_genotype_mat.shape[0])
	posterior_alpha = np.linalg.solve(posterior_sample_covariance, training_expr_vec_centered)
	posterior_mean_standardized_betas = prior_var*np.dot(np.transpose(standardized_training_genotype_mat), posterior_alpha)

	posterior_sample_covariance_inv_training_genotype = np.linalg.solve(posterior_sample_covariance, standardized_training_genotype_mat)
	posterior_var_standardized_betas = prior_var - np.square(prior_var)*np.sum(standardized_training_genotype_mat*posterior_sample_covariance_inv_training_genotype, axis=0)
	posterior_var_standardized_betas = np.maximum(posterior_var_standardized_betas, 0.0)

	posterior_mean_betas = posterior_mean_standardized_betas/snp_sdevs
	if training_genotype_mat.shape[1] > full_covariance_snp_limit:
		posterior_covariance_mat = np.diag(posterior_var_standardized_betas/np.square(snp_sdevs))
	else:
		posterior_covariance_mat = prior_var*np.eye(training_genotype_mat.shape[1]) - np.square(prior_var)*np.dot(np.transpose(standardized_training_genotype_mat), posterior_sample_covariance_inv_training_genotype)
		posterior_covariance_mat = posterior_covariance_mat/(snp_sdevs[:, None]*snp_sdevs[None, :])
		posterior_covariance_mat[np.diag_indices_from(posterior_covariance_mat)] = posterior_var_standardized_betas/np.square(snp_sdevs)

	return posterior_mean_betas, posterior_covariance_mat, prior_var

def compute_mean_causal_effects_from_distribution(borzoi_vec, anno, distribution_obj):
	if len(borzoi_vec) != anno.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	if distribution_obj.get('model_effect_scale') == 'allelic_grid_squared':
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

def estimate_snp_eqtl_effects_and_ses(genotype_mat, expr_vec):
	X = np.asarray(genotype_mat, dtype=float)
	y = np.asarray(expr_vec, dtype=float)
	n_samples = X.shape[0]
	n_snps = X.shape[1]
	if n_samples <= 2:
		return np.full(n_snps, np.nan), np.full(n_snps, np.nan)

	X = X - np.mean(X, axis=0)
	y = y - np.mean(y)
	sxx = np.sum(np.square(X), axis=0)
	valid_snps = np.isfinite(sxx) & (sxx > 0.0)

	eqtl_effect_sizes = np.full(n_snps, np.nan)
	eqtl_ses = np.full(n_snps, np.nan)
	eqtl_effect_sizes[valid_snps] = np.dot(y, X[:, valid_snps])/sxx[valid_snps]

	residuals = y[:, None] - X[:, valid_snps]*eqtl_effect_sizes[valid_snps][None, :]
	residual_var = np.sum(np.square(residuals), axis=0)/(n_samples - 2)
	eqtl_ses[valid_snps] = np.sqrt(residual_var/sxx[valid_snps])
	return eqtl_effect_sizes, eqtl_ses

def run_twas(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, gene_id_to_expression_vector, genotype_stem, distribution_obj, distribution, output_file, standardize_geno, training_expression_sample_indices, training_genotype_sample_indices, variant_to_gwas_sumstats):
	# Initialize output file
	t = open(output_file,'w')
	t.write('gene_id\tn_snps\tn_gwas_snps\tn_ld_imputed_gwas_snps\tmax_abs_gwas_z_score\thybrid_twas_z_score\thybrid_twas_se\thybrid_twas_pvalue\tno_borzoi_twas_z_score\tno_borzoi_twas_se\tno_borzoi_twas_pvalue\tborzoi_twas_z_score\tborzoi_twas_se\tborzoi_twas_pvalue\thybrid_twas_positive_probability\tborzoi_twas_positive_probability\thybrid_directional_fsr\tborzoi_directional_fsr\tcis_snp_h2\tcis_snp_h2_pvalue\n')


	# Loop through chromsomes
	for chrom_num in range(1,23):

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
			if distribution != 'ldscore_grid_squared' and gene_id not in gene_id_to_variant_gene_anno:
				continue

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])

			gwas_effect_sizes = []
			gwas_ses = []
			for variant in ordered_cis_variants:
				if variant in variant_to_gwas_sumstats:
					snp_info = rsid_to_snp_info[variant]
					geno_alleles = snp_info[:2]
					if geno_alleles[0] != variant.split('_')[3]:
						print('assumption erroror')
						pdb.set_trace()
					gwas_effect_sizes.append(variant_to_gwas_sumstats[variant][0])
					gwas_ses.append(variant_to_gwas_sumstats[variant][1])
				else:
					gwas_effect_sizes.append(np.nan)
					gwas_ses.append(np.nan)

			gwas_effect_sizes = np.asarray(gwas_effect_sizes)
			gwas_ses = np.asarray(gwas_ses)

			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 10:
				continue

			# Load in data for gene
			# Borzoi
			borzoi_effects_unstandardized, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])

			# Anno
			if distribution == 'ldscore_grid_squared':
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
				if distribution != 'ldscore_grid_squared' and borzoi_variant_alleles[var_index,:][0] != borzoi_anno_variant_alleles[var_index,:][0]:
					print('annotation alllele assumption erororo')
					pdb.set_trace()
			
			# Extract genotype
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = G[cis_genotype_indices,:].compute()
			training_geno_mat = geno_mat[:, training_genotype_sample_indices]
			# Mean impute using training samples only
			training_row_means = np.nanmean(training_geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(training_geno_mat))
			training_geno_mat[nan_rows, nan_cols] = training_row_means[nan_rows]

			training_geno_mat = 2.0 - training_geno_mat

			snp_means = np.mean(training_geno_mat, axis=1)
			snp_sdevs = np.std(training_geno_mat, axis=1)
			zero_var_snps = snp_sdevs == 0.0
			valid_training_snps = np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
			snp_sdevs[zero_var_snps] = 1.0
			if standardize_geno == 'True':
				print('assumption erroor')
				pdb.set_trace()
			else:
				training_geno_mat = (training_geno_mat - snp_means[:, None])
				borzoi_effects_standardized = borzoi_effects_unstandardized

			expr_vec = gene_id_to_expression_vector[gene_id]
			training_expr_vec = expr_vec[training_expression_sample_indices]


			if distribution != 'ldscore_grid_squared':
				print('not yet implemented')
				pdb.set_trace()

			if np.sum(valid_training_snps) < 10:
				continue

			cis_snp_h2, cis_snp_h2_pvalue = estimate_cis_snp_heritability_with_lrt(np.transpose(training_geno_mat[valid_training_snps, :]), training_expr_vec)

			snp_sdevs = snp_sdevs[valid_training_snps]
			borzoi_vec = borzoi_effects_standardized[valid_training_snps]
			training_genotype_mat = np.transpose(training_geno_mat[valid_training_snps, :])
			eqtl_effect_sizes, eqtl_ses = estimate_snp_eqtl_effects_and_ses(training_genotype_mat, training_expr_vec)

			gwas_effect_sizes = gwas_effect_sizes[valid_training_snps]
			gwas_ses = gwas_ses[valid_training_snps]
			valid_gwas_snps = np.isfinite(gwas_effect_sizes) & np.isfinite(gwas_ses) & (gwas_ses > 0.0)
			if np.sum(valid_gwas_snps) == 0:
				max_abs_gwas_z_score = np.nan
				n_ld_imputed_gwas_snps = 0
			else:
				max_abs_gwas_z_score = np.max(np.abs(gwas_effect_sizes[valid_gwas_snps]/gwas_ses[valid_gwas_snps]))
				n_ld_imputed_gwas_snps = len(gwas_effect_sizes) - np.sum(valid_gwas_snps)

			borzoi_mean_betas, borzoi_covariance_mat = compute_ldscore_grid_squared_prior_distribution(borzoi_vec, distribution_obj)
			borzoi_directional_fsr = compute_directional_fsr_from_genotype(training_genotype_mat, borzoi_mean_betas, borzoi_covariance_mat)

			hybrid_posterior_mean_betas, hybrid_posterior_covariance_mat = compute_ldscore_grid_squared_posterior_distribution(training_genotype_mat, training_expr_vec, borzoi_vec, distribution_obj)
			no_borzoi_posterior_mean_betas, no_borzoi_posterior_covariance_mat, optimal_prior_var = compute_no_borzoi_standardized_gaussian_posterior_distribution(training_genotype_mat, training_expr_vec, snp_sdevs)
			borzoi_mean_betas, borzoi_covariance_mat = compute_ldscore_grid_squared_prior_distribution(borzoi_vec, distribution_obj)
			hybrid_directional_fsr = compute_directional_fsr_from_genotype(training_genotype_mat, hybrid_posterior_mean_betas, hybrid_posterior_covariance_mat)

			if training_genotype_mat.shape[1] == 1:
				ld_mat = np.ones((1, 1))
			else:
				ld_mat = np.corrcoef(training_genotype_mat, rowvar=False)
			hybrid_twas_effect, hybrid_twas_se, hybrid_twas_pvalue = run_twas_from_per_allele_sumstats(gwas_effect_sizes, gwas_ses, hybrid_posterior_mean_betas, ld_mat, snp_sdevs)
			no_borzoi_twas_effect, no_borzoi_twas_se, no_borzoi_twas_pvalue = run_twas_from_per_allele_sumstats(gwas_effect_sizes, gwas_ses, no_borzoi_posterior_mean_betas, ld_mat, snp_sdevs)
			borzoi_twas_effect, borzoi_twas_se, borzoi_twas_pvalue = run_twas_from_per_allele_sumstats(gwas_effect_sizes, gwas_ses, borzoi_mean_betas, ld_mat, snp_sdevs)
			hybrid_twas_positive_probability = compute_twas_positive_probability_from_effect_distribution(gwas_effect_sizes, gwas_ses, hybrid_posterior_mean_betas, hybrid_posterior_covariance_mat, ld_mat, snp_sdevs)
			borzoi_twas_positive_probability = compute_twas_positive_probability_from_effect_distribution(gwas_effect_sizes, gwas_ses, borzoi_mean_betas, borzoi_covariance_mat, ld_mat, snp_sdevs)

			t.write(gene_id + '\t' + str(len(gwas_effect_sizes)) + '\t' + str(np.sum(valid_gwas_snps)) + '\t' + str(n_ld_imputed_gwas_snps) + '\t' + str(max_abs_gwas_z_score) + '\t' + str(hybrid_twas_effect) + '\t' + str(hybrid_twas_se) + '\t' + str(hybrid_twas_pvalue) + '\t' + str(no_borzoi_twas_effect) + '\t' + str(no_borzoi_twas_se) + '\t' + str(no_borzoi_twas_pvalue) + '\t' + str(borzoi_twas_effect) + '\t' + str(borzoi_twas_se) + '\t' + str(borzoi_twas_pvalue) + '\t' + str(hybrid_twas_positive_probability) + '\t' + str(borzoi_twas_positive_probability) + '\t' + str(hybrid_directional_fsr) + '\t' + str(borzoi_directional_fsr) + '\t' + str(cis_snp_h2) + '\t' + str(cis_snp_h2_pvalue) + '\n')
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


def load_in_ldscore_grid_squared_distribution_data(borzoi_based_prior_output_stem):
	summary_file = borzoi_based_prior_output_stem + '_ldscore_grid_squared_summary.txt'
	param_file = borzoi_based_prior_output_stem + '_ldscore_grid_squared_params.txt'

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
		print('ldscore_grid_squared prior loading error: could not find grid_knots in ' + summary_file)
		pdb.set_trace()
	if model_effect_scale != 'allelic_grid_squared':
		print('ldscore_grid_squared prior loading error: expected model_effect_scale=allelic_grid_squared')
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
		print('ldscore_grid_squared prior loading error: could not find a_prior in ' + param_file)
		pdb.set_trace()

	ordered_coef_indices = np.sort(np.asarray([*grid_coef_mapping]))
	if np.array_equal(ordered_coef_indices, np.arange(len(grid_knots) + 2)) == False:
		print('ldscore_grid_squared prior loading error: grid coefficient indices do not match expected basis size')
		pdb.set_trace()
	grid_coefs = np.asarray([grid_coef_mapping[index] for index in ordered_coef_indices])

	distribution_obj = {}
	distribution_obj['model_effect_scale'] = model_effect_scale
	distribution_obj['a_prior'] = a_prior
	distribution_obj['grid_knots'] = grid_knots
	distribution_obj['grid_coefs'] = grid_coefs
	return distribution_obj

def create_mapping_from_variant_id_to_gwas_ss(gwas_ss_file):
	f = gzip.open(gwas_ss_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data[4]) > 1 or len(data[5]) > 1:
			continue
		var_id1 = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5] + '_b38'
		var_id2 = 'chr' + data[1] + '_' + data[2] + '_' + data[5] + '_' + data[4] + '_b38'
		effect_size = float(data[10])
		effect_size_se = float(data[11])
		mapping[var_id1] = (effect_size, effect_size_se)
		mapping[var_id2] = (-1.0*effect_size, effect_size_se)

	f.close()

	return mapping

def create_temp_variant_to_gwas_sumstats_file(output_file):
	output_stem = os.path.splitext(output_file)[0]
	return output_stem + '_temp_variant_to_gwas_sumstats.pkl'



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
gwas_ss_file = sys.argv[10]
training_sample_size = int(sys.argv[11])


# Set seed
np.random.seed(1)

###########################
# Load in data
###########################



# Load in gwas data
variant_to_gwas_sumstats = create_mapping_from_variant_id_to_gwas_ss(gwas_ss_file)
'''
variant_to_gwas_sumstats_file = create_temp_variant_to_gwas_sumstats_file(output_file)
with open(variant_to_gwas_sumstats_file, 'wb') as f:
	pickle.dump(variant_to_gwas_sumstats, f, protocol=pickle.HIGHEST_PROTOCOL)
print('Saved variant_to_gwas_sumstats to ' + variant_to_gwas_sumstats_file)
'''
'''
with open(variant_to_gwas_sumstats_file, 'rb') as f:
	variant_to_gwas_sumstats = pickle.load(f)
print('Loaded variant_to_gwas_sumstats from ' + variant_to_gwas_sumstats_file)
'''

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)
expression_sample_indices = np.arange(len(genotype_sample_indices))

if training_sample_size < len(expression_sample_indices):
	downsample_indices = np.random.choice(np.arange(len(genotype_sample_indices)), size=training_sample_size, replace=False)
	genotype_sample_indices = genotype_sample_indices[downsample_indices]
	expression_sample_indices = expression_sample_indices[downsample_indices]
else:
	print('Requested training sample size ' + str(training_sample_size) + ' but only ' + str(len(expression_sample_indices)) + ' expression samples are available; using all available samples.')




if distribution == 'gaussian':
	distribution_obj = load_in_gaussian_distribution_data(borzoi_based_prior_output_stem)
elif distribution == 'ashr':
	distribution_obj = load_in_ashr_distribution_data(borzoi_based_prior_output_stem)
elif distribution == 'effect_size_grid':
	distribution_obj = load_in_effect_size_grid_distribution_data(borzoi_based_prior_output_stem)
elif distribution.startswith('uniform_prior_grid'):
	distribution_obj = load_in_uniform_prior_grid_distribution_data(borzoi_based_prior_output_stem)
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
if distribution == 'ldscore_grid_squared':
	gene_id_to_variant_gene_anno = None
	anno_names = np.asarray([])
else:
	gene_id_to_variant_gene_anno, anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(borzoi_annotation_file)


# Create mapping from gene id to expression vector
gene_id_to_expression_vector = create_mapping_from_gene_id_to_expression_vector(expr_file)


# run expression prediction
run_twas(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, gene_id_to_expression_vector, genotype_stem, distribution_obj, distribution, output_file, standardize_geno, expression_sample_indices, genotype_sample_indices, variant_to_gwas_sumstats)
