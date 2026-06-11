import numpy as np
import os
import sys
import pdb
from statistics import NormalDist
from pandas_plink import read_plink


def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	mapping = {}

	for snp_iter, snp_name in enumerate(ordered_snps):
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


def infer_genotype_stem(processed_genotype_data_dir):
	processed_genotype_data_dir = processed_genotype_data_dir.rstrip('/')
	if processed_genotype_data_dir.endswith('gtex_v9_eqtl_chr'):
		return processed_genotype_data_dir
	return os.path.join(processed_genotype_data_dir, 'gtex_v9_eqtl_chr')


def mean_impute_genotype_matrix(geno_mat):
	row_means = np.nanmean(geno_mat, axis=1)
	nan_rows, nan_cols = np.where(np.isnan(geno_mat))
	geno_mat[nan_rows, nan_cols] = row_means[nan_rows]
	return geno_mat


def select_random_gene_tss_from_chromosome(bim, cis_window=100000, min_variants=10, rng=None):
	if rng is None:
		rng = np.random.default_rng()

	positions = np.asarray(bim['pos']).astype(int)
	if len(positions) == 0:
		raise ValueError('No variants found on chromosome')

	sort_indices = np.argsort(positions)
	sorted_positions = positions[sort_indices]
	left_edges = np.searchsorted(sorted_positions, sorted_positions - cis_window, side='left')
	right_edges = np.searchsorted(sorted_positions, sorted_positions + cis_window, side='right')
	n_cis_variants = right_edges - left_edges
	candidate_sorted_indices = np.where(n_cis_variants >= min_variants)[0]
	if len(candidate_sorted_indices) == 0:
		raise ValueError('No chromosome positions have at least ' + str(min_variants) + ' variants in the cis window')

	chosen_sorted_index = rng.choice(candidate_sorted_indices)
	chosen_bim_index = sort_indices[chosen_sorted_index]
	return positions[chosen_bim_index]


def load_genotype_data_for_gene_window(processed_genotype_data_dir, chrom_num, gene_tss, cis_window=100000, genotype_sample_indices=None, standardize_genotype=True):
	genotype_stem = infer_genotype_stem(processed_genotype_data_dir)
	bim, fam, G = read_plink(genotype_stem + str(chrom_num))

	ordered_snps = np.asarray(bim['snp'])
	snp_positions = np.asarray(bim['pos']).astype(int)
	rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(ordered_snps)
	rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(
		ordered_snps,
		np.asarray(bim['a0']),
		np.asarray(bim['a1']),
		np.asarray(bim['chrom']),
		snp_positions,
	)

	ordered_cis_variants = ordered_snps[np.abs(snp_positions - int(gene_tss)) <= cis_window]
	if len(ordered_cis_variants) == 0:
		raise ValueError('No cis variants found for chromosome ' + str(chrom_num) + ' TSS ' + str(gene_tss))

	cis_genotype_indices = np.asarray([rsid_to_genotype_index[variant_id] for variant_id in ordered_cis_variants])
	geno_mat = G[cis_genotype_indices, :].compute()
	if genotype_sample_indices is not None:
		genotype_sample_indices = np.asarray(genotype_sample_indices).astype(int)
		geno_mat = geno_mat[:, genotype_sample_indices]

	geno_mat = mean_impute_genotype_matrix(geno_mat)
	geno_mat = 2.0 - geno_mat
	unstandardized_genotype_mat = geno_mat.copy()

	snp_means = np.mean(geno_mat, axis=1)
	snp_sdevs = np.std(geno_mat, axis=1)
	valid_snps = np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
	if np.sum(valid_snps) == 0:
		raise ValueError('No non-constant cis variants found for chromosome ' + str(chrom_num) + ' TSS ' + str(gene_tss))

	if standardize_genotype:
		geno_mat = (geno_mat - snp_means[:, None])/snp_sdevs[:, None]
	else:
		geno_mat = geno_mat - snp_means[:, None]

	if 'iid' in fam:
		sample_ids = np.asarray(fam['iid'])
	else:
		sample_ids = np.asarray(fam.iloc[:, 1])
	if genotype_sample_indices is not None:
		sample_ids = sample_ids[genotype_sample_indices]

	return {
		'gene_id': 'chr' + str(chrom_num) + '_random_gene_' + str(int(gene_tss)),
		'chrom_num': chrom_num,
		'gene_tss': int(gene_tss),
		'ordered_variants': ordered_cis_variants[valid_snps],
		'variant_positions': np.asarray([rsid_to_snp_info[variant_id][3] for variant_id in ordered_cis_variants])[valid_snps].astype(int),
		'variant_alleles': np.asarray([(rsid_to_snp_info[variant_id][0], rsid_to_snp_info[variant_id][1]) for variant_id in ordered_cis_variants])[valid_snps],
		'genotype_mat': np.transpose(geno_mat[valid_snps, :]),
		'unstandardized_genotype_mat': np.transpose(unstandardized_genotype_mat[valid_snps, :]),
		'snp_means': snp_means[valid_snps],
		'snp_sdevs': snp_sdevs[valid_snps],
		'sample_ids': sample_ids,
	}


def load_genotype_data_for_random_gene_on_chrom_10(processed_genotype_data_dir, cis_window=100000, min_variants=10, genotype_sample_indices=None, standardize_genotype=True, random_seed=None):
	rng = np.random.default_rng(random_seed)
	genotype_stem = infer_genotype_stem(processed_genotype_data_dir)
	bim, _, _ = read_plink(genotype_stem + '10')
	gene_tss = select_random_gene_tss_from_chromosome(bim, cis_window=cis_window, min_variants=min_variants, rng=rng)
	return load_genotype_data_for_gene_window(
		processed_genotype_data_dir,
		10,
		gene_tss,
		cis_window=cis_window,
		genotype_sample_indices=genotype_sample_indices,
		standardize_genotype=standardize_genotype,
	)

def simulate_noisy_estimated_eqtl_effects_with_target_fsr(std_genotype, causal_effects, target_fsr, rng=None):
	if rng is None:
		rng = np.random.default_rng()
	if target_fsr < 0.0 or target_fsr >= 0.5:
		raise ValueError('target_fsr must be in [0, 0.5)')

	true_pred_expr = np.dot(std_genotype, causal_effects)
	true_signal = np.dot(true_pred_expr, true_pred_expr)
	if true_signal <= 0.0:
		raise ValueError('Cannot calibrate noisy effects when true predicted expression has zero variance')

	if target_fsr == 0.0:
		return causal_effects.copy(), 0.0, 0.0

	noise_projection_sd_per_unit = np.sqrt(np.sum(np.square(np.dot(np.transpose(std_genotype), true_pred_expr))))
	if noise_projection_sd_per_unit <= 0.0:
		raise ValueError('Cannot calibrate noisy effects because genotype-noise projection has zero variance')

	z_score = NormalDist().inv_cdf(1.0 - target_fsr)
	noise_sd = true_signal/(z_score*noise_projection_sd_per_unit)
	noisy_estimated_effects = causal_effects + rng.normal(loc=0.0, scale=noise_sd, size=len(causal_effects))
	expected_fsr = NormalDist().cdf(-true_signal/(noise_sd*noise_projection_sd_per_unit))

	return noisy_estimated_effects, noise_sd, expected_fsr

def compute_ld_matrix(std_genotype):
	ld_mat = np.corrcoef(std_genotype, rowvar=False)
	if np.any(np.isfinite(ld_mat) == False):
		raise ValueError('LD matrix contains non-finite values')
	return ld_mat

def simulate_multivariate_normal_noise(covariance_mat, rng):
	eigenvalues, eigenvectors = np.linalg.eigh(covariance_mat)
	eigenvalues[eigenvalues < 0.0] = 0.0
	return np.dot(eigenvectors, np.sqrt(eigenvalues)*rng.normal(size=len(eigenvalues)))

def simulate_gwas_z_scores_using_ld(std_genotype, causal_effects, sim_alpha, gwas_ss, rng=None):
	if rng is None:
		rng = np.random.default_rng()
	if gwas_ss <= 0.0:
		raise ValueError('gwas_ss must be positive')

	ld_mat = compute_ld_matrix(std_genotype)
	genetic_component = np.dot(std_genotype, causal_effects)
	genetic_component_mean = np.mean(genetic_component)
	genetic_component_sdev = np.std(genetic_component)
	if genetic_component_sdev <= 0.0 or np.isfinite(genetic_component_sdev) == False:
		raise ValueError('Cannot simulate GWAS effects when genetic expression component has zero variance')
	standardized_genetic_component = (genetic_component - genetic_component_mean)/genetic_component_sdev
	gwas_causal_effects = sim_alpha*causal_effects/genetic_component_sdev
	marginal_effect_mean = np.dot(ld_mat, gwas_causal_effects)
	marginal_effects = marginal_effect_mean + simulate_multivariate_normal_noise(ld_mat/gwas_ss, rng)
	gwas_z_scores = np.sqrt(gwas_ss)*marginal_effects

	return {
		'ld_mat': ld_mat,
		'standardized_genetic_component': standardized_genetic_component,
		'gwas_causal_effects': gwas_causal_effects,
		'gwas_marginal_effect_mean': marginal_effect_mean,
		'gwas_marginal_effects': marginal_effects,
		'gwas_z_scores': gwas_z_scores,
	}

def compute_twas_association_statistics(noisy_estimated_effects, gwas_marginal_effects, gwas_z_scores, ld_mat):
	pred_expr_var = np.dot(noisy_estimated_effects, np.dot(ld_mat, noisy_estimated_effects))
	if pred_expr_var <= 0.0 or np.isfinite(pred_expr_var) == False:
		raise ValueError('Cannot compute TWAS statistic with non-positive predicted expression variance')

	twas_effect_size = np.dot(noisy_estimated_effects, gwas_marginal_effects)/pred_expr_var
	twas_z_score = np.dot(noisy_estimated_effects, gwas_z_scores)/np.sqrt(pred_expr_var)
	twas_pvalue = 2.0*(1.0 - NormalDist().cdf(np.abs(twas_z_score)))

	return twas_effect_size, twas_z_score, twas_pvalue

def logistic_function(x):
	return 1.0/(1.0 + np.exp(-x))

def log_normal_density(x, mean, var):
	return -0.5*(np.log(2.0*np.pi*var) + (np.square(x - mean)/var))

def simulate_causal_eqtl_effects_from_delta_prior(n_snps, baseline_inclusion_prob=0.05, delta_inclusion_slope=1.0, delta_prior_slope=0.02, delta_prior_resid_sd=0.01, rng=None):
	if rng is None:
		rng = np.random.default_rng()
	delta = rng.normal(loc=0.0, scale=1.0, size=n_snps)
	baseline_log_odds = np.log(baseline_inclusion_prob/(1.0 - baseline_inclusion_prob))
	inclusion_probs = logistic_function(baseline_log_odds + delta_inclusion_slope*delta)
	causal_indicators = rng.uniform(size=n_snps) < inclusion_probs
	nonzero_effects = rng.normal(loc=delta_prior_slope*delta, scale=delta_prior_resid_sd, size=n_snps)
	causal_effects = np.zeros(n_snps)
	causal_effects[causal_indicators] = nonzero_effects[causal_indicators]
	return causal_effects, delta, causal_indicators, inclusion_probs

def compute_posterior_beta_moments_given_delta(noisy_estimated_effects, eqtl_effect_noise_sd, delta, baseline_inclusion_prob=0.05, delta_inclusion_slope=1.0, delta_prior_slope=0.02, delta_prior_resid_sd=0.01):
	baseline_log_odds = np.log(baseline_inclusion_prob/(1.0 - baseline_inclusion_prob))
	prior_inclusion_probs = logistic_function(baseline_log_odds + delta_inclusion_slope*delta)
	slab_prior_mean = delta_prior_slope*delta
	slab_prior_var = np.square(delta_prior_resid_sd)
	if slab_prior_var <= 0.0:
		raise ValueError('delta_prior_resid_sd must be positive')
	if eqtl_effect_noise_sd == 0.0:
		return noisy_estimated_effects.copy(), np.zeros(len(noisy_estimated_effects))

	noise_var = np.square(eqtl_effect_noise_sd)
	slab_posterior_var = 1.0/((1.0/slab_prior_var) + (1.0/noise_var))
	slab_posterior_mean = slab_posterior_var*((slab_prior_mean/slab_prior_var) + (noisy_estimated_effects/noise_var))

	log_null_weight = np.log(1.0 - prior_inclusion_probs) + log_normal_density(noisy_estimated_effects, 0.0, noise_var)
	log_slab_weight = np.log(prior_inclusion_probs) + log_normal_density(noisy_estimated_effects, slab_prior_mean, slab_prior_var + noise_var)
	max_log_weight = np.maximum(log_null_weight, log_slab_weight)
	null_weight = np.exp(log_null_weight - max_log_weight)
	slab_weight = np.exp(log_slab_weight - max_log_weight)
	posterior_inclusion_probs = slab_weight/(null_weight + slab_weight)

	posterior_mean = posterior_inclusion_probs*slab_posterior_mean
	posterior_second_moment = posterior_inclusion_probs*(slab_posterior_var + np.square(slab_posterior_mean))
	posterior_var = posterior_second_moment - np.square(posterior_mean)
	return posterior_mean, posterior_var

def standardize_vector(vec):
	vec_mean = np.mean(vec)
	vec_sdev = np.std(vec)
	if vec_sdev <= 0.0 or np.isfinite(vec_sdev) == False:
		raise ValueError('Cannot standardize vector with non-positive standard deviation')
	return (vec - vec_mean)/vec_sdev, vec_mean, vec_sdev

def simulate_individual_trait(std_genotype, causal_effects, sim_alpha, sigma_y_sq=1.0, rng=None):
	if rng is None:
		rng = np.random.default_rng()
	true_genetic_component, _, _ = standardize_vector(np.dot(std_genotype, causal_effects))
	trait = sim_alpha*true_genetic_component + rng.normal(loc=0.0, scale=np.sqrt(sigma_y_sq), size=std_genotype.shape[0])
	return trait, true_genetic_component

def construct_genetic_expression_uncertainty_model(std_genotype, posterior_beta_mean, posterior_beta_var):
	raw_mu_g = np.dot(std_genotype, posterior_beta_mean)
	mu_g, _, mu_g_sdev = standardize_vector(raw_mu_g)
	if np.all(posterior_beta_var == 0.0):
		k_g = np.zeros((std_genotype.shape[0], std_genotype.shape[0]))
	else:
		k_g = np.dot(std_genotype*posterior_beta_var[None, :], np.transpose(std_genotype))/np.square(mu_g_sdev)
	return mu_g, k_g

def compute_uncertainty_aware_alpha_log_likelihood(alpha, transformed_trait, transformed_mu_g, k_g_eigenvalues, sigma_y_sq):
	conditional_variance = sigma_y_sq + np.square(alpha)*k_g_eigenvalues
	if np.any(conditional_variance <= 0.0):
		return -np.inf
	resid = transformed_trait - alpha*transformed_mu_g
	return -0.5*(np.sum(np.log(conditional_variance)) + np.sum(np.square(resid)/conditional_variance))

def fit_uncertainty_aware_twas_model(trait, mu_g, k_g, sigma_y_sq=1.0, alpha_grid=None):
	if alpha_grid is None:
		alpha_grid = np.linspace(-0.25, 0.25, 5001)
	k_g_eigenvalues, k_g_eigenvectors = np.linalg.eigh(k_g)
	k_g_eigenvalues[k_g_eigenvalues < 0.0] = 0.0
	transformed_trait = np.dot(np.transpose(k_g_eigenvectors), trait)
	transformed_mu_g = np.dot(np.transpose(k_g_eigenvectors), mu_g)

	log_likelihoods = []
	for alpha in alpha_grid:
		log_likelihoods.append(compute_uncertainty_aware_alpha_log_likelihood(alpha, transformed_trait, transformed_mu_g, k_g_eigenvalues, sigma_y_sq))
	log_likelihoods = np.asarray(log_likelihoods)

	best_index = np.argmax(log_likelihoods)
	alpha_hat = alpha_grid[best_index]
	null_log_likelihood = compute_uncertainty_aware_alpha_log_likelihood(0.0, transformed_trait, transformed_mu_g, k_g_eigenvalues, sigma_y_sq)
	lrt_stat = np.maximum(2.0*(log_likelihoods[best_index] - null_log_likelihood), 0.0)
	lrt_pvalue = 2.0*(1.0 - NormalDist().cdf(np.sqrt(lrt_stat)))

	posterior_weights = np.exp(log_likelihoods - np.max(log_likelihoods))
	p_alpha_gt_zero = np.sum(posterior_weights[alpha_grid > 0.0])/np.sum(posterior_weights)
	directional_fsr = np.minimum(p_alpha_gt_zero, 1.0 - p_alpha_gt_zero)

	return {
		'alpha_hat': alpha_hat,
		'lrt_stat': lrt_stat,
		'lrt_pvalue': lrt_pvalue,
		'p_alpha_gt_zero': p_alpha_gt_zero,
		'directional_fsr': directional_fsr,
		'alpha_grid': alpha_grid,
		'log_likelihoods': log_likelihoods,
	}

def run_twas_simulation(std_genotype, sim_fsr, sim_alpha_var, gwas_ss):
	sign_concordances = []
	posterior_mean_sign_concordances = []
	uncertainty_sign_concordances = []
	uncertainty_directional_fsrs = []
	hit_pvalue_threshold = .05/10000
	effective_trait_noise_var = std_genotype.shape[0]/gwas_ss
	for itera in range(100):
		print(itera)
		# Simulate causal effects
		n_snps = std_genotype.shape[1]
		causal_effects, delta, causal_indicators, prior_inclusion_probs = simulate_causal_eqtl_effects_from_delta_prior(n_snps)

		noisy_estimated_effects, noisy_effect_noise_sd, expected_fsr = simulate_noisy_estimated_eqtl_effects_with_target_fsr(std_genotype, causal_effects, sim_fsr)
		posterior_beta_mean, posterior_beta_var = compute_posterior_beta_moments_given_delta(noisy_estimated_effects, noisy_effect_noise_sd, delta)

		corry = np.corrcoef(np.dot(std_genotype, causal_effects), np.dot(std_genotype, noisy_estimated_effects))[0,1]

		sim_alpha = np.random.normal(loc=0, scale=np.sqrt(sim_alpha_var))
		gwas_simulation = simulate_gwas_z_scores_using_ld(std_genotype, causal_effects, sim_alpha, gwas_ss)
		gwas_z_scores = gwas_simulation['gwas_z_scores']
		twas_effect_size, twas_z_score, twas_pvalue = compute_twas_association_statistics(
			noisy_estimated_effects,
			gwas_simulation['gwas_marginal_effects'],
			gwas_z_scores,
			gwas_simulation['ld_mat'],
		)
		posterior_mean_twas_effect_size, posterior_mean_twas_z_score, posterior_mean_twas_pvalue = compute_twas_association_statistics(
			posterior_beta_mean,
			gwas_simulation['gwas_marginal_effects'],
			gwas_z_scores,
			gwas_simulation['ld_mat'],
		)

		trait, true_genetic_component = simulate_individual_trait(std_genotype, causal_effects, sim_alpha, sigma_y_sq=effective_trait_noise_var)
		mu_g, k_g = construct_genetic_expression_uncertainty_model(std_genotype, posterior_beta_mean, posterior_beta_var)
		uncertainty_twas_result = fit_uncertainty_aware_twas_model(trait, mu_g, k_g, sigma_y_sq=effective_trait_noise_var)

		if twas_pvalue < hit_pvalue_threshold:
			if twas_effect_size*sim_alpha > 0:
				sign_concordances.append(1.0)
			else:
				sign_concordances.append(-1.0)
		if posterior_mean_twas_pvalue < hit_pvalue_threshold:
			if posterior_mean_twas_effect_size*sim_alpha > 0:
				posterior_mean_sign_concordances.append(1.0)
			else:
				posterior_mean_sign_concordances.append(-1.0)
		if uncertainty_twas_result['lrt_pvalue'] < hit_pvalue_threshold:
			if uncertainty_twas_result['alpha_hat']*sim_alpha > 0:
				uncertainty_sign_concordances.append(1.0)
			else:
				uncertainty_sign_concordances.append(-1.0)
			uncertainty_directional_fsrs.append(uncertainty_twas_result['directional_fsr'])
	sign_concordances = np.asarray(sign_concordances)
	posterior_mean_sign_concordances = np.asarray(posterior_mean_sign_concordances)
	uncertainty_sign_concordances = np.asarray(uncertainty_sign_concordances)
	uncertainty_directional_fsrs = np.asarray(uncertainty_directional_fsrs)
	if len(sign_concordances) > 0:
		print('plugin_twas_hits\t' + str(len(sign_concordances)) + '\tplugin_false_sign_rate\t' + str(np.sum(sign_concordances < 0.0)/len(sign_concordances)))
	else:
		print('plugin_twas_hits\t0\tplugin_false_sign_rate\tnan')
	if len(posterior_mean_sign_concordances) > 0:
		print('posterior_mean_twas_hits\t' + str(len(posterior_mean_sign_concordances)) + '\tposterior_mean_false_sign_rate\t' + str(np.sum(posterior_mean_sign_concordances < 0.0)/len(posterior_mean_sign_concordances)))
	else:
		print('posterior_mean_twas_hits\t0\tposterior_mean_false_sign_rate\tnan')
	if len(uncertainty_sign_concordances) > 0:
		print('uncertainty_twas_hits\t' + str(len(uncertainty_sign_concordances)) + '\tuncertainty_false_sign_rate\t' + str(np.sum(uncertainty_sign_concordances < 0.0)/len(uncertainty_sign_concordances)) + '\tavg_model_directional_fsr\t' + str(np.mean(uncertainty_directional_fsrs)))
	else:
		print('uncertainty_twas_hits\t0\tuncertainty_false_sign_rate\tnan\tavg_model_directional_fsr\tnan')
	pdb.set_trace()





def main():
	######################
	# Command line args
	######################
	processed_genotype_data_dir = sys.argv[1]
	twas_simulation_output_stem = sys.argv[2]
	sim_fsr = float(sys.argv[3])
	sim_alpha_var = float(sys.argv[4])
	gwas_ss = float(sys.argv[5])

	random_gene_genotype_data = load_genotype_data_for_random_gene_on_chrom_10(processed_genotype_data_dir)
	print(
		random_gene_genotype_data['gene_id'] + '\t' +
		str(random_gene_genotype_data['genotype_mat'].shape[0]) + '\t' +
		str(random_gene_genotype_data['genotype_mat'].shape[1])
	)

	run_twas_simulation(random_gene_genotype_data['genotype_mat'], sim_fsr, sim_alpha_var, gwas_ss)




if __name__ == '__main__':
	main()
