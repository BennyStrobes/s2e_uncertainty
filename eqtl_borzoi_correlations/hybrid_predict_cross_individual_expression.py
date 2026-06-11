import numpy as np
import os
import sys
import pdb
import gzip
import pickle
import warnings
from pandas_plink import read_plink
import time
from scipy.optimize import minimize_scalar
from scipy import stats
from sklearn.linear_model import ElasticNet
from sklearn.exceptions import ConvergenceWarning

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

def compute_prediction_correlation_from_effect_mean(validation_genotype_mat, validation_expr_vec, effect_mean):
	predicted_expr = np.dot(validation_genotype_mat, effect_mean)
	prediction_correlation = np.corrcoef(validation_expr_vec, predicted_expr)[0,1]
	_, correlation_pvalue = stats.pearsonr(validation_expr_vec, predicted_expr)
	fisher_z = np.arctanh(prediction_correlation)
	fisher_z_se = 1.0/np.sqrt(len(validation_expr_vec) - 3)
	correlation_prob_negative = stats.norm.cdf(0.0, loc=fisher_z, scale=fisher_z_se)

	return prediction_correlation, correlation_pvalue, correlation_prob_negative

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

def compute_elastic_net_alpha_grid(genotype_mat, expr_vec, l1_ratio, n_alphas=5, eps=1e-4):
	n_samples = genotype_mat.shape[0]
	if n_samples == 0:
		return np.asarray([1.0])

	alpha_max = np.max(np.abs(np.dot(np.transpose(genotype_mat), expr_vec)))/(n_samples*np.maximum(l1_ratio, 1e-12))
	if np.isfinite(alpha_max) == False or alpha_max <= 0.0:
		alpha_max = 1.0
	return alpha_max*np.power(eps, np.linspace(0.0, 1.0, n_alphas))

def fit_elastic_net_effect_mean(genotype_mat, expr_vec, alpha, l1_ratio, fail_on_convergence_warning=False):
	model = ElasticNet(
		alpha=alpha,
		l1_ratio=l1_ratio,
		fit_intercept=False,
		max_iter=200000,
		tol=1e-3,
		selection='random',
		random_state=1
	)
	if fail_on_convergence_warning:
		with warnings.catch_warnings():
			warnings.filterwarnings('error', category=ConvergenceWarning)
			model.fit(genotype_mat, expr_vec)
	else:
		model.fit(genotype_mat, expr_vec)
	return model.coef_

def learn_standardized_elastic_net_hyperparameters_with_cross_validation(training_genotype_mat, training_expr_vec, l1_ratio_grid=None, n_alphas=5, eps=1e-4, n_folds=4):
	if l1_ratio_grid is None:
		l1_ratio_grid = np.asarray([0.01, 0.1, 0.5, 0.9])

	n_samples = training_genotype_mat.shape[0]
	n_folds = np.min([n_folds, n_samples])
	training_expr_vec_centered = training_expr_vec - np.mean(training_expr_vec)
	if n_folds < 2:
		return compute_elastic_net_alpha_grid(training_genotype_mat, training_expr_vec_centered, l1_ratio_grid[0], n_alphas=n_alphas, eps=eps)[0], l1_ratio_grid[0]

	shuffled_indices = np.random.permutation(n_samples)
	fold_indices = np.array_split(shuffled_indices, n_folds)
	all_indices = np.arange(n_samples)
	best_cv_error = np.inf
	best_alpha = None
	best_l1_ratio = None
	for l1_ratio in l1_ratio_grid:
		alpha_grid = compute_elastic_net_alpha_grid(training_genotype_mat, training_expr_vec_centered, l1_ratio, n_alphas=n_alphas, eps=eps)
		for alpha in alpha_grid:
			fold_errors = []
			failed_fit = False
			for validation_indices in fold_indices:
				training_indices = np.setdiff1d(all_indices, validation_indices)
				if len(training_indices) < 2 or len(validation_indices) == 0:
					continue
				training_expr_mean = np.mean(training_expr_vec[training_indices])
				fold_training_expr_vec = training_expr_vec[training_indices] - training_expr_mean
				fold_validation_expr_vec = training_expr_vec[validation_indices] - training_expr_mean
				fold_training_genotype_mat = training_genotype_mat[training_indices, :]
				fold_validation_genotype_mat = training_genotype_mat[validation_indices, :]
				try:
					effect_mean = fit_elastic_net_effect_mean(fold_training_genotype_mat, fold_training_expr_vec, alpha, l1_ratio, fail_on_convergence_warning=True)
					predicted_expr = np.dot(fold_validation_genotype_mat, effect_mean)
					fold_errors.append(np.mean(np.square(fold_validation_expr_vec - predicted_expr)))
				except (Exception, ConvergenceWarning):
					failed_fit = True
					break
			if failed_fit or len(fold_errors) == 0:
				continue
			cv_error = np.mean(fold_errors)
			if cv_error < best_cv_error:
				best_cv_error = cv_error
				best_alpha = alpha
				best_l1_ratio = l1_ratio

	if best_alpha is None:
		best_l1_ratio = l1_ratio_grid[0]
		best_alpha = compute_elastic_net_alpha_grid(training_genotype_mat, training_expr_vec_centered, best_l1_ratio, n_alphas=n_alphas, eps=eps)[0]

	return best_alpha, best_l1_ratio

def compute_no_borzoi_standardized_elastic_net_effect_mean(training_genotype_mat, training_expr_vec, snp_sdevs):
	standardized_training_genotype_mat = training_genotype_mat/snp_sdevs[None, :]
	alpha, l1_ratio = learn_standardized_elastic_net_hyperparameters_with_cross_validation(standardized_training_genotype_mat, training_expr_vec)

	training_expr_vec_centered = training_expr_vec - np.mean(training_expr_vec)
	try:
		effect_mean_standardized = fit_elastic_net_effect_mean(standardized_training_genotype_mat, training_expr_vec_centered, alpha, l1_ratio)
	except Exception as e:
		print('ElasticNet failed: ' + str(e), flush=True)
		effect_mean_standardized = np.zeros(standardized_training_genotype_mat.shape[1])

	effect_mean = effect_mean_standardized/snp_sdevs
	return effect_mean, alpha, l1_ratio

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

def run_expression_prediction(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, gene_id_to_expression_vector, genotype_stem, distribution_obj, distribution, output_file, standardize_geno, training_expression_sample_indices, training_genotype_sample_indices, validation_expression_sample_indices, validation_genotype_sample_indices):
	# Initialize output file
	t = open(output_file,'w')
	t.write('gene_id\thybrid_correlation\tno_borzoi_correlation\telastic_net_correlation\tborzoi_correlation\thybrid_directional_FSR\tno_borzoi_directional_FSR\telastic_net_directional_FSR\tborzoi_directional_FSR\thybrid_correlation_pvalue\tno_borzoi_correlation_pvalue\telastic_net_correlation_pvalue\tborzoi_correlation_pvalue\thybrid_correlation_prob_negative\tno_borzoi_correlation_prob_negative\telastic_net_correlation_prob_negative\tborzoi_correlation_prob_negative\tcis_snp_h2\tcis_snp_h2_pvalue\n')


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
			if distribution != 'ldscore_grid_squared' and gene_id not in gene_id_to_variant_gene_anno:
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
			validation_geno_mat = geno_mat[:, validation_genotype_sample_indices]
			# Mean impute using training samples only
			training_row_means = np.nanmean(training_geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(training_geno_mat))
			training_geno_mat[nan_rows, nan_cols] = training_row_means[nan_rows]

			nan_rows, nan_cols = np.where(np.isnan(validation_geno_mat))
			validation_geno_mat[nan_rows, nan_cols] = training_row_means[nan_rows]

			training_geno_mat = 2.0 - training_geno_mat
			validation_geno_mat = 2.0 - validation_geno_mat

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
				validation_geno_mat = (validation_geno_mat - snp_means[:, None])
				borzoi_effects_standardized = borzoi_effects_unstandardized

			expr_vec = gene_id_to_expression_vector[gene_id]
			training_expr_vec = expr_vec[training_expression_sample_indices]
			validation_expr_vec = expr_vec[validation_expression_sample_indices]


			if distribution != 'ldscore_grid_squared':
				print('not yet implemented')
				pdb.set_trace()


			cis_snp_h2, cis_snp_h2_pvalue = estimate_cis_snp_heritability_with_lrt(np.transpose(training_geno_mat[valid_training_snps, :]), training_expr_vec)


			borzoi_vec = borzoi_effects_standardized[valid_training_snps]
			training_genotype_mat = np.transpose(training_geno_mat[valid_training_snps, :])
			validation_genotype_mat = np.transpose(validation_geno_mat[valid_training_snps, :])

			hybrid_posterior_mean_betas, hybrid_posterior_covariance_mat = compute_ldscore_grid_squared_posterior_distribution(training_genotype_mat, training_expr_vec, borzoi_vec, distribution_obj)
			no_borzoi_posterior_mean_betas, no_borzoi_posterior_covariance_mat, optimal_prior_var = compute_no_borzoi_standardized_gaussian_posterior_distribution(training_genotype_mat, training_expr_vec, snp_sdevs[valid_training_snps])
			elastic_net_effect_mean_betas, _, _ = compute_no_borzoi_standardized_elastic_net_effect_mean(training_genotype_mat, training_expr_vec, snp_sdevs[valid_training_snps])
			borzoi_mean_betas, borzoi_covariance_mat = compute_ldscore_grid_squared_prior_distribution(borzoi_vec, distribution_obj)

			

			hybrid_correlation, hybrid_directional_fsr, hybrid_correlation_pvalue, hybrid_correlation_prob_negative = compute_prediction_correlation_and_directional_fsr(validation_genotype_mat, validation_expr_vec, hybrid_posterior_mean_betas, hybrid_posterior_covariance_mat)
			no_borzoi_correlation, no_borzoi_directional_fsr, no_borzoi_correlation_pvalue, no_borzoi_correlation_prob_negative = compute_prediction_correlation_and_directional_fsr(validation_genotype_mat, validation_expr_vec, no_borzoi_posterior_mean_betas, no_borzoi_posterior_covariance_mat)
			elastic_net_correlation, elastic_net_correlation_pvalue, elastic_net_correlation_prob_negative = compute_prediction_correlation_from_effect_mean(validation_genotype_mat, validation_expr_vec, elastic_net_effect_mean_betas)
			elastic_net_directional_fsr = np.nan
			borzoi_correlation, borzoi_directional_fsr, borzoi_correlation_pvalue, borzoi_correlation_prob_negative = compute_prediction_correlation_and_directional_fsr(validation_genotype_mat, validation_expr_vec, borzoi_mean_betas, borzoi_covariance_mat)

			t.write(gene_id + '\t' + str(hybrid_correlation) + '\t' + str(no_borzoi_correlation) + '\t' + str(elastic_net_correlation) + '\t' + str(borzoi_correlation) + '\t' + str(hybrid_directional_fsr) + '\t' + str(no_borzoi_directional_fsr) + '\t' + str(elastic_net_directional_fsr) + '\t' + str(borzoi_directional_fsr) + '\t' + str(hybrid_correlation_pvalue) + '\t' + str(no_borzoi_correlation_pvalue) + '\t' + str(elastic_net_correlation_pvalue) + '\t' + str(borzoi_correlation_pvalue) + '\t' + str(hybrid_correlation_prob_negative) + '\t' + str(no_borzoi_correlation_prob_negative) + '\t' + str(elastic_net_correlation_prob_negative) + '\t' + str(borzoi_correlation_prob_negative) + '\t' + str(cis_snp_h2) + '\t' + str(cis_snp_h2_pvalue) + '\n')
			t.flush()
			print(time.time())

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
training_sample_size = int(sys.argv[10])
max_training_sample_size = int(sys.argv[11])


# Set seed
np.random.seed(1)


###########################
# Load in data
###########################


# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)
expression_sample_indices = np.arange(len(genotype_sample_indices))



if max_training_sample_size > len(genotype_sample_indices):
	print('max_training_sample_size exceeds number of available samples')
	pdb.set_trace()


validation_sample_size = len(genotype_sample_indices) - max_training_sample_size

validation_expression_sample_indices = np.random.choice(expression_sample_indices, size=validation_sample_size, replace=False)
validation_genotype_sample_indices = genotype_sample_indices[validation_expression_sample_indices]

training_expression_sample_pool = np.setdiff1d(expression_sample_indices, validation_expression_sample_indices)
if training_sample_size > len(training_expression_sample_pool):
	print('training_sample_size exceeds number of non-validation samples')
	pdb.set_trace()
training_expression_sample_indices = np.random.choice(training_expression_sample_pool, size=training_sample_size, replace=False)
training_genotype_sample_indices = genotype_sample_indices[training_expression_sample_indices]


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
run_expression_prediction(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, gene_id_to_expression_vector, genotype_stem, distribution_obj, distribution, output_file, standardize_geno, training_expression_sample_indices, training_genotype_sample_indices, validation_expression_sample_indices, validation_genotype_sample_indices)
