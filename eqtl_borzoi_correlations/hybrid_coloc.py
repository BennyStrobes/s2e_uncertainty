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

def convert_gwas_sumstats_to_ld_imputed_sumstats(gwas_effect_sizes, gwas_ses, ld_mat, snp_sdevs, ridge=0.1):
	gwas_effect_sizes = np.asarray(gwas_effect_sizes, dtype=float)
	gwas_ses = np.asarray(gwas_ses, dtype=float)
	snp_sdevs = np.asarray(snp_sdevs, dtype=float)
	observed_gwas_snps = np.isfinite(gwas_effect_sizes) & np.isfinite(gwas_ses) & (gwas_ses > 0.0)
	gwas_z_scores = np.full(len(gwas_effect_sizes), np.nan)
	gwas_z_scores[observed_gwas_snps] = gwas_effect_sizes[observed_gwas_snps]/gwas_ses[observed_gwas_snps]
	gwas_z_scores, ld_imputed_gwas_snps, gwas_z_imputation_r2 = ld_impute_missing_gwas_z_scores(gwas_z_scores, ld_mat, ridge=ridge)

	imputed_gwas_effect_sizes = np.copy(gwas_effect_sizes)
	imputed_gwas_ses = np.copy(gwas_ses)
	observed_with_sdev = observed_gwas_snps & np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
	if np.sum(observed_with_sdev) > 0:
		estimated_gwas_sample_sizes = 1.0/(np.square(gwas_ses[observed_with_sdev])*np.square(snp_sdevs[observed_with_sdev]))
		estimated_gwas_sample_size = np.median(estimated_gwas_sample_sizes[np.isfinite(estimated_gwas_sample_sizes) & (estimated_gwas_sample_sizes > 0.0)])
	else:
		estimated_gwas_sample_size = np.nan
	if np.isfinite(estimated_gwas_sample_size) and estimated_gwas_sample_size > 0.0:
		imputed_with_sdev = ld_imputed_gwas_snps & np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
		imputed_gwas_ses[imputed_with_sdev] = 1.0/(np.sqrt(estimated_gwas_sample_size)*snp_sdevs[imputed_with_sdev])
		imputed_gwas_effect_sizes[imputed_with_sdev] = gwas_z_scores[imputed_with_sdev]*imputed_gwas_ses[imputed_with_sdev]
	return imputed_gwas_effect_sizes, imputed_gwas_ses, observed_gwas_snps, ld_imputed_gwas_snps, gwas_z_imputation_r2

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

def logdiffexp(log_x, log_y):
	if log_y >= log_x:
		return -np.inf
	return log_x + np.log1p(-np.exp(log_y - log_x))

def compute_snp_causal_posteriors_from_log_abfs(study1_log_abf, study2_log_abf, p1=1e-4, p2=1e-4, p12=1e-5):
	study1_log_abf = np.asarray(study1_log_abf, dtype=float)
	study2_log_abf = np.asarray(study2_log_abf, dtype=float)
	valid = np.isfinite(study1_log_abf) & np.isfinite(study2_log_abf)
	study1_snp_pp = np.full(len(study1_log_abf), np.nan)
	study2_snp_pp = np.full(len(study2_log_abf), np.nan)
	if np.sum(valid) == 0:
		return study1_snp_pp, study2_snp_pp, np.full(5, np.nan)

	l1 = study1_log_abf[valid]
	l2 = study2_log_abf[valid]
	pph = compute_coloc_pph_from_log_abfs(l1, l2, p1=p1, p2=p2, p12=p12)
	log_sum_1 = logsumexp(l1)
	log_sum_2 = logsumexp(l2)
	log_sum_shared = logsumexp(l1 + l2)
	if log_sum_shared >= log_sum_1 + log_sum_2:
		log_sum_distinct = -np.inf
	else:
		log_sum_distinct = log_sum_1 + log_sum_2 + np.log1p(-np.exp(log_sum_shared - log_sum_1 - log_sum_2))

	valid_study1_pp = pph[1]*np.exp(l1 - log_sum_1) + pph[4]*np.exp(l1 + l2 - log_sum_shared)
	valid_study2_pp = pph[2]*np.exp(l2 - log_sum_2) + pph[4]*np.exp(l1 + l2 - log_sum_shared)
	if np.isfinite(log_sum_distinct):
		study1_h3_pp = np.zeros(len(l1))
		study2_h3_pp = np.zeros(len(l2))
		for snp_iter in range(len(l1)):
			log_study2_except_snp = logdiffexp(log_sum_2, l2[snp_iter])
			log_study1_except_snp = logdiffexp(log_sum_1, l1[snp_iter])
			if np.isfinite(log_study2_except_snp):
				study1_h3_pp[snp_iter] = pph[3]*np.exp(l1[snp_iter] + log_study2_except_snp - log_sum_distinct)
			if np.isfinite(log_study1_except_snp):
				study2_h3_pp[snp_iter] = pph[3]*np.exp(l2[snp_iter] + log_study1_except_snp - log_sum_distinct)
		valid_study1_pp = valid_study1_pp + study1_h3_pp
		valid_study2_pp = valid_study2_pp + study2_h3_pp

	study1_snp_pp[valid] = valid_study1_pp
	study2_snp_pp[valid] = valid_study2_pp
	return study1_snp_pp, study2_snp_pp, pph

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

def load_in_coloc_data(gene_id_to_est_borzoi_effects, gene_id_to_expression_vector, genotype_stem, standardize_geno, training_expression_sample_indices, training_genotype_sample_indices, variant_to_gwas_sumstats, borzoi_bins):
	# Initialize output file
	#t = open(output_file,'w')
	#t.write('gene_id\thybrid_correlation\tno_borzoi_correlation\tborzoi_correlation\thybrid_directional_FSR\tno_borzoi_directional_FSR\tborzoi_directional_FSR\thybrid_correlation_pvalue\tno_borzoi_correlation_pvalue\tborzoi_correlation_pvalue\thybrid_correlation_prob_negative\tno_borzoi_correlation_prob_negative\tborzoi_correlation_prob_negative\tcis_snp_h2\tcis_snp_h2_pvalue\n')

	coloc_data = []

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

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])

			gwas_effect_sizes = []
			gwas_ses = []
			for variant in ordered_cis_variants:
				if variant in variant_to_gwas_sumstats:
					gwas_effect_sizes.append(variant_to_gwas_sumstats[variant][0])
					gwas_ses.append(variant_to_gwas_sumstats[variant][1])
				else:
					gwas_effect_sizes.append(np.nan)
					gwas_ses.append(np.nan)

			gwas_effect_sizes = np.asarray(gwas_effect_sizes)
			gwas_ses = np.asarray(gwas_ses)

			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 30:
				continue

			# Load in data for gene
			# Borzoi
			borzoi_effects_unstandardized, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])


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

			training_geno_mat = (training_geno_mat - snp_means[:, None])

			expr_vec = gene_id_to_expression_vector[gene_id]
			training_expr_vec = expr_vec[training_expression_sample_indices]



			borzoi_vec = borzoi_effects_unstandardized[valid_training_snps]
			valid_snp_sdevs = snp_sdevs[valid_training_snps]
			training_genotype_mat = np.transpose(training_geno_mat[valid_training_snps, :])
			eqtl_effect_sizes, eqtl_ses = estimate_snp_eqtl_effects_and_ses(training_genotype_mat, training_expr_vec)

			gwas_effect_sizes = gwas_effect_sizes[valid_training_snps]
			gwas_ses = gwas_ses[valid_training_snps]

			if training_genotype_mat.shape[1] == 1:
				ld_mat = np.ones((1, 1))
			else:
				ld_mat = np.corrcoef(training_genotype_mat, rowvar=False)
			gwas_effect_sizes, gwas_ses, observed_gwas_snps, ld_imputed_gwas_snps, gwas_z_imputation_r2 = convert_gwas_sumstats_to_ld_imputed_sumstats(gwas_effect_sizes, gwas_ses, ld_mat, valid_snp_sdevs)
			if np.sum(np.isnan(gwas_effect_sizes)==False) < 10:
				continue 

			discretized_borzoi_effects = []
			for effect in borzoi_vec:
				bin_assignment = np.nan
				for bin_counter in range(1, len(borzoi_bins)):
					if np.abs(effect) >= borzoi_bins[bin_counter-1] and np.abs(effect) < borzoi_bins[bin_counter]:
						bin_assignment= bin_counter -1
				discretized_borzoi_effects.append((bin_assignment, np.sign(effect)))
				if np.isnan(bin_assignment):
					print('assumption eroorr')
					pdb.set_trace()

			coloc_data.append((gene_id, gwas_effect_sizes, gwas_ses, eqtl_effect_sizes, eqtl_ses, discretized_borzoi_effects))

	return coloc_data



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
		mapping[var_id2] = (effect_size, effect_size_se)

	f.close()

	return mapping

def create_temp_variant_to_gwas_sumstats_file(output_file):
	output_stem = os.path.splitext(output_file)[0]
	return output_stem + '_temp_variant_to_gwas_sumstats.pkl'

def create_temp_coloc_data_file(output_file):
	output_stem = os.path.splitext(output_file)[0]
	return output_stem + '_temp_coloc_data.pkl'

def extract_gene_posteriors(coloc_data, study1_prior_var, study2_prior_var):
	gene_posteriors = []
	for gene_coloc_data in coloc_data:
		gene_id, study1_effect_sizes, study1_ses, study2_effect_sizes, study2_ses, discretized_borzoi_effects = gene_coloc_data
		study1_log_abf = compute_coloc_log_abf(study1_effect_sizes, study1_ses, prior_var=study1_prior_var)
		study2_log_abf = compute_coloc_log_abf(study2_effect_sizes, study2_ses, prior_var=study2_prior_var)
		study1_snp_pp, study2_snp_pp, coloc_pph = compute_snp_causal_posteriors_from_log_abfs(study1_log_abf, study2_log_abf)
		gene_posteriors.append((
			gene_id,
			study1_snp_pp,
			study2_snp_pp,
			coloc_pph,
			discretized_borzoi_effects
		))
	return gene_posteriors

def get_n_borzoi_bins_from_coloc_data(coloc_data):
	max_bin = -1
	for gene_coloc_data in coloc_data:
		discretized_borzoi_effects = gene_coloc_data[5]
		for bin_assignment, _ in discretized_borzoi_effects:
			if np.isfinite(bin_assignment):
				max_bin = max(max_bin, int(bin_assignment))
	return max_bin + 1

def make_borzoi_informed_prior_vectors(discretized_borzoi_effects, bin_prior_means, bin_prior_vars):
	prior_means = []
	prior_vars = []
	for bin_assignment, borzoi_sign in discretized_borzoi_effects:
		if np.isfinite(bin_assignment) == False:
			prior_means.append(np.nan)
			prior_vars.append(np.nan)
			continue
		bin_index = int(bin_assignment)
		prior_means.append(float(borzoi_sign)*bin_prior_means[bin_index])
		prior_vars.append(bin_prior_vars[bin_index])
	return np.asarray(prior_means), np.asarray(prior_vars)

def extract_gene_posteriors_borzoi_informed_eqtl(coloc_data, study1_prior_var, eqtl_bin_prior_means, eqtl_bin_prior_vars):
	gene_posteriors = []
	for gene_coloc_data in coloc_data:
		gene_id, study1_effect_sizes, study1_ses, study2_effect_sizes, study2_ses, discretized_borzoi_effects = gene_coloc_data
		eqtl_prior_means, eqtl_prior_vars = make_borzoi_informed_prior_vectors(discretized_borzoi_effects, eqtl_bin_prior_means, eqtl_bin_prior_vars)
		study1_log_abf = compute_coloc_log_abf(study1_effect_sizes, study1_ses, prior_var=study1_prior_var)
		study2_log_abf = compute_coloc_log_abf(study2_effect_sizes, study2_ses, prior_mean=eqtl_prior_means, prior_var=eqtl_prior_vars)
		study1_snp_pp, study2_snp_pp, coloc_pph = compute_snp_causal_posteriors_from_log_abfs(study1_log_abf, study2_log_abf)
		gene_posteriors.append((
			gene_id,
			study1_snp_pp,
			study2_snp_pp,
			coloc_pph,
			discretized_borzoi_effects
		))
	return gene_posteriors

def compute_prior_variance_update(effect_sizes, effect_ses, causal_snp_pp, current_prior_var):
	effect_sizes = np.asarray(effect_sizes, dtype=float)
	effect_vars = np.square(np.asarray(effect_ses, dtype=float))
	causal_snp_pp = np.asarray(causal_snp_pp, dtype=float)
	valid = np.isfinite(effect_sizes) & np.isfinite(effect_vars) & np.isfinite(causal_snp_pp) & (effect_vars > 0.0)
	if np.sum(valid) == 0 or np.sum(causal_snp_pp[valid]) <= 0.0:
		return 0.0, 0.0

	beta_hat = effect_sizes[valid]
	V = effect_vars[valid]
	posterior_var = (current_prior_var*V)/(current_prior_var + V)
	posterior_mean = (current_prior_var/(current_prior_var + V))*beta_hat
	posterior_second_moment = posterior_var + np.square(posterior_mean)
	weighted_second_moment = np.sum(causal_snp_pp[valid]*posterior_second_moment)
	expected_n_causal = np.sum(causal_snp_pp[valid])
	return weighted_second_moment, expected_n_causal

def update_study_prior_variances(coloc_data, gene_posteriors, study1_prior_var, study2_prior_var, min_prior_var=1e-12):
	study1_weighted_second_moment = 0.0
	study1_expected_n_causal = 0.0
	study2_weighted_second_moment = 0.0
	study2_expected_n_causal = 0.0

	for gene_coloc_data, gene_posterior in zip(coloc_data, gene_posteriors):
		_, study1_effect_sizes, study1_ses, study2_effect_sizes, study2_ses, _ = gene_coloc_data
		_, study1_snp_pp, study2_snp_pp, _, _ = gene_posterior

		weighted_second_moment, expected_n_causal = compute_prior_variance_update(
			study1_effect_sizes,
			study1_ses,
			study1_snp_pp,
			study1_prior_var
		)
		study1_weighted_second_moment = study1_weighted_second_moment + weighted_second_moment
		study1_expected_n_causal = study1_expected_n_causal + expected_n_causal

		weighted_second_moment, expected_n_causal = compute_prior_variance_update(
			study2_effect_sizes,
			study2_ses,
			study2_snp_pp,
			study2_prior_var
		)
		study2_weighted_second_moment = study2_weighted_second_moment + weighted_second_moment
		study2_expected_n_causal = study2_expected_n_causal + expected_n_causal

	new_study1_prior_var = study1_prior_var
	new_study2_prior_var = study2_prior_var
	if study1_expected_n_causal > 0.0:
		new_study1_prior_var = max(min_prior_var, study1_weighted_second_moment/study1_expected_n_causal)
	if study2_expected_n_causal > 0.0:
		new_study2_prior_var = max(min_prior_var, study2_weighted_second_moment/study2_expected_n_causal)
	return new_study1_prior_var, new_study2_prior_var

def update_borzoi_informed_eqtl_prior_params(coloc_data, gene_posteriors, eqtl_bin_prior_means, eqtl_bin_prior_vars, mean_prior_var=1.0, min_prior_var=1e-12):
	n_bins = len(eqtl_bin_prior_means)
	weighted_first_moments = np.zeros(n_bins)
	weighted_second_moments = np.zeros(n_bins)
	expected_n_causal = np.zeros(n_bins)

	for gene_coloc_data, gene_posterior in zip(coloc_data, gene_posteriors):
		_, _, _, eqtl_effect_sizes, eqtl_ses, discretized_borzoi_effects = gene_coloc_data
		_, _, eqtl_snp_pp, _, _ = gene_posterior
		eqtl_effect_sizes = np.asarray(eqtl_effect_sizes, dtype=float)
		eqtl_vars = np.square(np.asarray(eqtl_ses, dtype=float))
		eqtl_snp_pp = np.asarray(eqtl_snp_pp, dtype=float)

		for snp_iter, (bin_assignment, borzoi_sign) in enumerate(discretized_borzoi_effects):
			if np.isfinite(bin_assignment) == False or np.isfinite(borzoi_sign) == False or borzoi_sign == 0.0:
				continue
			if np.isfinite(eqtl_effect_sizes[snp_iter]) == False or np.isfinite(eqtl_vars[snp_iter]) == False or np.isfinite(eqtl_snp_pp[snp_iter]) == False or eqtl_vars[snp_iter] <= 0.0:
				continue
			bin_index = int(bin_assignment)
			prior_mean = eqtl_bin_prior_means[bin_index]
			prior_var = eqtl_bin_prior_vars[bin_index]
			signed_effect = float(borzoi_sign)*eqtl_effect_sizes[snp_iter]
			V = eqtl_vars[snp_iter]
			posterior_var = (prior_var*V)/(prior_var + V)
			posterior_mean = (prior_var/(prior_var + V))*signed_effect + (V/(prior_var + V))*prior_mean
			posterior_second_moment = posterior_var + np.square(posterior_mean)
			weighted_first_moments[bin_index] = weighted_first_moments[bin_index] + eqtl_snp_pp[snp_iter]*posterior_mean
			weighted_second_moments[bin_index] = weighted_second_moments[bin_index] + eqtl_snp_pp[snp_iter]*posterior_second_moment
			expected_n_causal[bin_index] = expected_n_causal[bin_index] + eqtl_snp_pp[snp_iter]

	new_bin_prior_means = np.copy(eqtl_bin_prior_means)
	new_bin_prior_vars = np.copy(eqtl_bin_prior_vars)
	for bin_index in range(n_bins):
		if expected_n_causal[bin_index] <= 0.0:
			continue
		mean_prior_precision = eqtl_bin_prior_vars[bin_index]/mean_prior_var
		new_bin_prior_means[bin_index] = weighted_first_moments[bin_index]/(expected_n_causal[bin_index] + mean_prior_precision)
		mean_centered_second_moment = (
			weighted_second_moments[bin_index]/expected_n_causal[bin_index] -
			2.0*new_bin_prior_means[bin_index]*weighted_first_moments[bin_index]/expected_n_causal[bin_index] +
			np.square(new_bin_prior_means[bin_index])
		)
		new_bin_prior_vars[bin_index] = max(
			min_prior_var,
			mean_centered_second_moment
		)
	return new_bin_prior_means, new_bin_prior_vars


def bayesian_coloc(coloc_data, version='gaussian_priors', max_iter=100, convergence_thresh=1e-6):
	study1_prior_var = .15*.15
	study2_prior_var = .15*.15

	if version == 'borzoi_informed_eqtl_priors':
		n_borzoi_bins = get_n_borzoi_bins_from_coloc_data(coloc_data)
		eqtl_bin_prior_means = np.zeros(n_borzoi_bins)
		eqtl_bin_prior_vars = np.full(n_borzoi_bins, study2_prior_var)
		for itera in range(max_iter):
			gene_posteriors = extract_gene_posteriors_borzoi_informed_eqtl(
				coloc_data,
				study1_prior_var,
				eqtl_bin_prior_means,
				eqtl_bin_prior_vars
			)
			new_study1_prior_var = study1_prior_var
			study1_weighted_second_moment = 0.0
			study1_expected_n_causal = 0.0
			for gene_coloc_data, gene_posterior in zip(coloc_data, gene_posteriors):
				_, study1_effect_sizes, study1_ses, _, _, _ = gene_coloc_data
				_, study1_snp_pp, _, _, _ = gene_posterior
				weighted_second_moment, expected_n_causal = compute_prior_variance_update(
					study1_effect_sizes,
					study1_ses,
					study1_snp_pp,
					study1_prior_var
				)
				study1_weighted_second_moment = study1_weighted_second_moment + weighted_second_moment
				study1_expected_n_causal = study1_expected_n_causal + expected_n_causal
			if study1_expected_n_causal > 0.0:
				new_study1_prior_var = max(1e-12, study1_weighted_second_moment/study1_expected_n_causal)
			new_eqtl_bin_prior_means, new_eqtl_bin_prior_vars = update_borzoi_informed_eqtl_prior_params(
				coloc_data,
				gene_posteriors,
				eqtl_bin_prior_means,
				eqtl_bin_prior_vars
			)
			prior_delta = max(
				np.max(np.abs(new_eqtl_bin_prior_means - eqtl_bin_prior_means)),
				np.max(np.abs(new_eqtl_bin_prior_vars - eqtl_bin_prior_vars)),
				np.abs(new_study1_prior_var - study1_prior_var)
			)
			study1_prior_var = new_study1_prior_var
			eqtl_bin_prior_means = new_eqtl_bin_prior_means
			eqtl_bin_prior_vars = new_eqtl_bin_prior_vars
			print('iter ' + str(itera) + '\tstudy1_prior_var=' + str(study1_prior_var) + '\teqtl_bin_prior_means=' + ','.join(eqtl_bin_prior_means.astype(str)) + '\teqtl_bin_prior_vars=' + ','.join(eqtl_bin_prior_vars.astype(str)))
			if prior_delta < convergence_thresh:
				break
		return gene_posteriors, study1_prior_var, eqtl_bin_prior_means, eqtl_bin_prior_vars

	if version != 'gaussian_priors':
		print('bayesian_coloc version not implemented: ' + version)
		pdb.set_trace()

	for itera in range(max_iter):

		gene_posteriors = extract_gene_posteriors(coloc_data, study1_prior_var, study2_prior_var)
		new_study1_prior_var, new_study2_prior_var = update_study_prior_variances(
			coloc_data,
			gene_posteriors,
			study1_prior_var,
			study2_prior_var
		)
		prior_var_delta = max(
			np.abs(new_study1_prior_var - study1_prior_var),
			np.abs(new_study2_prior_var - study2_prior_var)
		)
		study1_prior_var = new_study1_prior_var
		study2_prior_var = new_study2_prior_var
		print('iter ' + str(itera) + '\tstudy1_prior_var=' + str(study1_prior_var) + '\tstudy2_prior_var=' + str(study2_prior_var))
		if prior_var_delta < convergence_thresh:
			break
	return gene_posteriors, study1_prior_var, study2_prior_var

def extract_pph4_by_gene(gene_posteriors):
	pph4_by_gene = {}
	for gene_posterior in gene_posteriors:
		gene_id = gene_posterior[0]
		coloc_pph = gene_posterior[3]
		if len(coloc_pph) != 5:
			pph4_by_gene[gene_id] = np.nan
		else:
			pph4_by_gene[gene_id] = coloc_pph[4]
	return pph4_by_gene

def write_coloc_pph4_summary(output_file, coloc_data, standard_gene_posteriors, gaussian_gene_posteriors, borzoi_gene_posteriors):
	standard_pph4 = extract_pph4_by_gene(standard_gene_posteriors)
	gaussian_pph4 = extract_pph4_by_gene(gaussian_gene_posteriors)
	borzoi_pph4 = extract_pph4_by_gene(borzoi_gene_posteriors)

	t = open(output_file, 'w')
	t.write('gene_id\tstandard_default_prior_pph4\tgaussian_empirical_prior_pph4\tborzoi_informed_eqtl_prior_pph4\n')
	for gene_coloc_data in coloc_data:
		gene_id = gene_coloc_data[0]
		t.write(
			gene_id + '\t' +
			str(standard_pph4.get(gene_id, np.nan)) + '\t' +
			str(gaussian_pph4.get(gene_id, np.nan)) + '\t' +
			str(borzoi_pph4.get(gene_id, np.nan)) + '\n'
		)
	t.close()
	print('Saved coloc PPH4 summary to ' + output_file)



########################
# Command line args
########################
borzoi_effect_file = sys.argv[1]
borzoi_annotation_file = sys.argv[2]
genotype_stem = sys.argv[3]
genotype_sample_mapping_file = sys.argv[4]
expr_file = sys.argv[5]
output_stem = sys.argv[6]
standardize_geno = sys.argv[7]
training_sample_size = int(sys.argv[8])
gwas_ss_file = sys.argv[9]

# Set seed
np.random.seed(1)

###########################
# Load in data
###########################
borzoi_bins = np.asarray([0, .01, .075, .2, 100000])
# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)

# Load in gwas data
variant_to_gwas_sumstats = create_mapping_from_variant_id_to_gwas_ss(gwas_ss_file)


# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)
expression_sample_indices = np.arange(len(genotype_sample_indices))



training_expression_sample_indices = np.random.choice(expression_sample_indices, size=training_sample_size, replace=False)
training_genotype_sample_indices = genotype_sample_indices[training_expression_sample_indices]



# Create mapping from gene id to expression vector
gene_id_to_expression_vector = create_mapping_from_gene_id_to_expression_vector(expr_file)


# run expression prediction
coloc_data = load_in_coloc_data(gene_id_to_est_borzoi_effects, gene_id_to_expression_vector, genotype_stem, standardize_geno, training_expression_sample_indices, training_genotype_sample_indices, variant_to_gwas_sumstats, borzoi_bins)


'''
coloc_data_file = create_temp_coloc_data_file(output_stem + '.txt')

with open(coloc_data_file, 'wb') as f:
	pickle.dump(coloc_data, f, protocol=pickle.HIGHEST_PROTOCOL)
print('Saved coloc_data to ' + coloc_data_file)

with open(coloc_data_file, 'rb') as f:
	coloc_data = pickle.load(f)
print('Loaded coloc_data from ' + coloc_data_file)

'''

# bayesian coloc
standard_gene_posteriors = extract_gene_posteriors(coloc_data, .15*.15, .15*.15)
borzoi_gene_posteriors, borzoi_gwas_prior_var, borzoi_eqtl_bin_prior_means, borzoi_eqtl_bin_prior_vars = bayesian_coloc(coloc_data, version='borzoi_informed_eqtl_priors')
gaussian_gene_posteriors, gaussian_study1_prior_var, gaussian_study2_prior_var = bayesian_coloc(coloc_data, version='gaussian_priors')
write_coloc_pph4_summary(output_stem, coloc_data, standard_gene_posteriors, gaussian_gene_posteriors, borzoi_gene_posteriors)




