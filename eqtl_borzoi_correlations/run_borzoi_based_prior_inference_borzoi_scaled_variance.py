import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time
from scipy.optimize import minimize

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


def generate_input_data_object_for_bayesian_inf(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, standardize_geno):
	# Initialize output object
	gene_to_data = {}
	
	# Loop through chromsomes
	for chrom_num in range(18,23):
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

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])
			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 10:
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
						print('need to flix flip. should flip genotypes not snps')
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
				borzoi_effects_for_model = borzoi_effects_unstandardized*snp_sdevs
			elif standardize_geno == 'False':
				geno_mat = (geno_mat - snp_means[:, None])
				borzoi_effects_for_model = np.copy(borzoi_effects_unstandardized)

			gene_to_data[gene_id] = {}
			gene_to_data[gene_id]['expression'] = gene_id_to_expression_vector[gene_id]
			gene_to_data[gene_id]['borzoi'] = borzoi_effects_for_model[valid_snps]
			gene_to_data[gene_id]['genotype'] = np.ascontiguousarray(geno_mat[valid_snps, :])

	return gene_to_data



def initialize_data(inf_input_data_obj):
	gene_names = np.asarray([*inf_input_data_obj])
	n_genes = len(gene_names)

	# Initialize priors
	a_prior = 1.0
	b_prior = 1.0/300.0
	c_prior = 0.0

	resid_vars = np.ones(n_genes)

	# Initialize causal effects
	betas = []
	sig_sq_priors = []

	expression_resids = []

	for gene_name in gene_names:

		borzoi_preds = inf_input_data_obj[gene_name]['borzoi']

		n_snps = len(borzoi_preds)

		if n_snps < 10:
			print('asusmption error: too few snps')
			pdb.set_trace()

		init_gene_causal_effects = np.zeros(len(borzoi_preds))
		betas.append(init_gene_causal_effects)

		sig_sq_priors.append(b_prior + (c_prior*np.square(borzoi_preds)))

		resid_E = np.copy(inf_input_data_obj[gene_name]['expression']) - np.dot(init_gene_causal_effects, inf_input_data_obj[gene_name]['genotype'])
		expression_resids.append(resid_E)

	return betas, a_prior, b_prior, c_prior, sig_sq_priors, resid_vars, expression_resids, gene_names


@njit(cache=True)
def update_causal_effect_for_single_gene_shell(cur_gene_beta, prior_means, per_snp_sig_sq_prior, gene_resid_var, resid_expr_vec, genotype_mat):
	# Number of snps per gene
	n_snps = len(prior_means)
	n_samples = len(resid_expr_vec)

	# Looop through snps in gene (in random order)
	for snp_iter in np.random.permutation(n_snps):
		# Add the current SNP contribution back into the residual before sampling it
		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] + (geno_val*cur_gene_beta[snp_iter])

		prior_var = per_snp_sig_sq_prior[snp_iter]
		prior_mean = prior_means[snp_iter]

		sum_sq_geno = 0.0
		geno_resid_dot = 0.0
		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			sum_sq_geno = sum_sq_geno + (geno_val*geno_val)
			geno_resid_dot = geno_resid_dot + (geno_val*resid_expr_vec[sample_iter])

		likelihood_precision = sum_sq_geno/gene_resid_var
		prior_precision = 1.0/prior_var

		posterior_var = 1.0/(prior_precision + likelihood_precision)
		posterior_mean = posterior_var*((prior_mean/prior_var) + (geno_resid_dot/gene_resid_var))

		new_beta = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))

		# Remove the updated SNP contribution from the residual
		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] - (geno_val*new_beta)
		cur_gene_beta[snp_iter] = new_beta

	return cur_gene_beta, resid_expr_vec

def update_causal_effect_for_single_gene(cur_gene_beta, a_prior, per_snp_sig_sq_prior, gene_resid_var, resid_expr_vec, borzoi_vec, genotype_mat):
	prior_means = a_prior*borzoi_vec
	return update_causal_effect_for_single_gene_shell(cur_gene_beta, prior_means, per_snp_sig_sq_prior, gene_resid_var, resid_expr_vec, genotype_mat)

def update_resid_var_for_single_gene(resid_expr_vec, v0=1e-10, s_sq0=1e-10):
	n_samples = len(resid_expr_vec)
	shape_post = v0 + (n_samples/2.0)
	scale_post = s_sq0 + (np.sum(np.square(resid_expr_vec))/2.0)
	return 1.0/np.random.gamma(shape_post, 1.0/scale_post)

@njit(cache=True)
def compute_posterior_beta_moments_for_single_gene_shell(cur_beta_mean, prior_means, per_snp_sig_sq_prior, gene_resid_var, expression_vec, genotype_mat, n_coordinate_sweeps):
	n_snps = len(prior_means)
	n_samples = len(expression_vec)
	beta_means = np.copy(cur_beta_mean)
	posterior_vars = np.zeros(n_snps)
	sum_sq_geno_arr = np.zeros(n_snps)
	resid_expr_vec = np.copy(expression_vec)

	# Initialize residual to y - X^T beta_mean and cache per-SNP genotype sums of squares.
	for snp_iter in range(n_snps):
		beta_val = beta_means[snp_iter]
		sum_sq_geno = 0.0
		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			sum_sq_geno = sum_sq_geno + (geno_val*geno_val)
			resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] - (geno_val*beta_val)
		sum_sq_geno_arr[snp_iter] = sum_sq_geno

	for _ in range(n_coordinate_sweeps):
		for snp_iter in range(n_snps):
			beta_old = beta_means[snp_iter]
			for sample_iter in range(n_samples):
				geno_val = genotype_mat[snp_iter, sample_iter]
				resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] + (geno_val*beta_old)

			prior_var = per_snp_sig_sq_prior[snp_iter]
			if prior_var < 1e-12:
				prior_var = 1e-12
			prior_mean = prior_means[snp_iter]
			sum_sq_geno = sum_sq_geno_arr[snp_iter]
			geno_resid_dot = 0.0
			for sample_iter in range(n_samples):
				geno_val = genotype_mat[snp_iter, sample_iter]
				geno_resid_dot = geno_resid_dot + (geno_val*resid_expr_vec[sample_iter])

			likelihood_precision = sum_sq_geno/gene_resid_var
			prior_precision = 1.0/prior_var
			posterior_var = 1.0/(prior_precision + likelihood_precision)
			posterior_mean = posterior_var*((prior_mean/prior_var) + (geno_resid_dot/gene_resid_var))
			posterior_vars[snp_iter] = posterior_var
			beta_means[snp_iter] = posterior_mean

			for sample_iter in range(n_samples):
				geno_val = genotype_mat[snp_iter, sample_iter]
				resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] - (geno_val*posterior_mean)

	expected_resid_ss = 0.0
	for sample_iter in range(n_samples):
		expected_resid_ss = expected_resid_ss + (resid_expr_vec[sample_iter]*resid_expr_vec[sample_iter])
	for snp_iter in range(n_snps):
		expected_resid_ss = expected_resid_ss + (posterior_vars[snp_iter]*sum_sq_geno_arr[snp_iter])

	return beta_means, posterior_vars, resid_expr_vec, expected_resid_ss

def compute_posterior_beta_moments_for_single_gene(cur_beta_mean, prior_means, per_snp_sig_sq_prior, gene_resid_var, expression_vec, genotype_mat, v0=1e-10, s_sq0=1e-10, n_coordinate_sweeps=3):
	posterior_mean, posterior_var, resid_mean, expected_resid_ss = compute_posterior_beta_moments_for_single_gene_shell(
		cur_beta_mean,
		prior_means,
		per_snp_sig_sq_prior,
		gene_resid_var,
		expression_vec,
		genotype_mat,
		n_coordinate_sweeps
	)
	shape_post = v0 + (len(expression_vec)/2.0)
	scale_post = s_sq0 + (expected_resid_ss/2.0)
	if shape_post <= 1.0:
		gene_resid_var_est = scale_post
	else:
		gene_resid_var_est = scale_post/(shape_post - 1.0)
	return posterior_mean, posterior_var, resid_mean, gene_resid_var_est

def update_causal_effect_estimates(betas, a_prior, sig_sq_priors, resid_vars, expression_resids, gene_names, inf_input_data_obj, n_coordinate_sweeps=3):
	# Loop through genes (perform update seperately for each gene)
	posterior_beta_vars = []
	for gene_iter, gene_name in enumerate(gene_names):

		#####################
		# Load in data for this gene

		prior_means = a_prior*inf_input_data_obj[gene_name]['borzoi']
		gene_beta_mean, gene_beta_var, gene_expression_resid, gene_resid_var_est = compute_posterior_beta_moments_for_single_gene(
			betas[gene_iter],
			prior_means,
			sig_sq_priors[gene_iter],
			resid_vars[gene_iter],
			inf_input_data_obj[gene_name]['expression'],
			inf_input_data_obj[gene_name]['genotype'],
			n_coordinate_sweeps=n_coordinate_sweeps
		)
		betas[gene_iter] = gene_beta_mean
		posterior_beta_vars.append(gene_beta_var)
		expression_resids[gene_iter] = gene_expression_resid
		resid_vars[gene_iter] = gene_resid_var_est

	return betas, posterior_beta_vars, expression_resids, resid_vars


def scaled_variance_objective_and_gradient(params, delta_sq, expected_sq_resid):
	b_prior = params[0]
	c_prior = params[1]
	sig_sq = b_prior + (c_prior*delta_sq)

	if b_prior < 1e-12 or c_prior < 0.0 or np.any(sig_sq <= 0.0):
		return np.inf, np.zeros(len(params))

	objective = 0.5*np.sum(np.log(sig_sq) + (expected_sq_resid/sig_sq))
	common_grad = 0.5*((1.0/sig_sq) - (expected_sq_resid/np.square(sig_sq)))
	gradient = np.asarray([
		np.sum(common_grad),
		np.sum(common_grad*delta_sq)
	])
	return objective, gradient

def optimize_scaled_variance_params(delta_sq, expected_sq_resid, b_prior, c_prior):
	x0 = np.asarray([max(b_prior, 1e-12), max(c_prior, 0.0)])

	def objective_wrapper(params):
		return scaled_variance_objective_and_gradient(params, delta_sq, expected_sq_resid)

	fit_res = minimize(
		fun=objective_wrapper,
		x0=x0,
		method='L-BFGS-B',
		jac=True,
		bounds=[(1e-12, None), (0.0, None)],
		options={'maxiter': 200}
	)
	if fit_res.success == False:
		print('Warning: b,c optimization did not fully converge: ' + str(fit_res.message), flush=True)
	return fit_res.x[0], fit_res.x[1]

def update_priors(beta_means, beta_vars, gene_names, inf_input_data_obj, b_prior, c_prior, tau_sq=100.0, max_m_step_iters=5):
	all_beta_mean = []
	all_beta_var = []
	all_borzoi = []
	gene_lengths = []
	for gene_iter, gene_name in enumerate(gene_names):
		gene_beta_mean = beta_means[gene_iter]
		gene_beta_var = beta_vars[gene_iter]
		all_beta_mean.append(gene_beta_mean)
		all_beta_var.append(gene_beta_var)
		all_borzoi.append(inf_input_data_obj[gene_name]['borzoi'])
		gene_lengths.append(len(gene_beta_mean))

	all_beta_mean = np.hstack(all_beta_mean)
	all_beta_var = np.hstack(all_beta_var)
	all_borzoi = np.hstack(all_borzoi)
	all_borzoi_sq = np.square(all_borzoi)

	a_prior = 0.0
	for _ in range(max_m_step_iters):
		sig_sq = b_prior + (c_prior*all_borzoi_sq)
		sig_sq[sig_sq < 1e-12] = 1e-12

		a_denom = np.sum(all_borzoi_sq/sig_sq) + (1.0/tau_sq)
		a_numer = np.sum((all_borzoi*all_beta_mean)/sig_sq)
		a_prior = a_numer/a_denom

		expected_sq_resid = all_beta_var + np.square(all_beta_mean - (a_prior*all_borzoi))
		b_prior, c_prior = optimize_scaled_variance_params(all_borzoi_sq, expected_sq_resid, b_prior, c_prior)

	sig_sq = b_prior + (c_prior*all_borzoi_sq)
	sig_sq[sig_sq < 1e-12] = 1e-12
	a_denom = np.sum(all_borzoi_sq/sig_sq) + (1.0/tau_sq)
	a_numer = np.sum((all_borzoi*all_beta_mean)/sig_sq)
	a_prior = a_numer/a_denom

	sig_sq_priors = []
	cur_index = 0
	for gene_length in gene_lengths:
		sig_sq_priors.append(sig_sq[cur_index:(cur_index + gene_length)])
		cur_index = cur_index + gene_length

	signal = np.var(a_prior*all_borzoi)
	avg_prior_var = np.mean(sig_sq)
	r_squares = np.asarray([signal/(signal + avg_prior_var)])

	return a_prior, b_prior, c_prior, sig_sq_priors, r_squares

def save_posterior_samples(output_stem, posterior_samples, param_names):
	sample_output_file = output_stem + '_posterior_samples.txt'

	t = open(sample_output_file,'w')
	t.write('sample_iter\t' + '\t'.join(param_names) + '\n')
	for sample_iter in range(posterior_samples.shape[0]):
		t.write(str(sample_iter) + '\t' + '\t'.join(posterior_samples[sample_iter, :].astype(str)) + '\n')
	t.close()

	summary_output_file = output_stem + '_posterior_summary.txt'
	t = open(summary_output_file,'w')
	t.write('parameter\tmean\tstandard_error\tci_2.5\tci_97.5\n')
	for param_iter, param_name in enumerate(param_names):
		param_samples = posterior_samples[:, param_iter]
		param_mean = np.mean(param_samples)
		if len(param_samples) > 1:
			param_se = np.std(param_samples, ddof=1)/np.sqrt(len(param_samples))
		else:
			param_se = 0.0
		param_lb = np.percentile(param_samples, 2.5)
		param_ub = np.percentile(param_samples, 97.5)
		t.write(param_name + '\t' + str(param_mean) + '\t' + str(param_se) + '\t' + str(param_lb) + '\t' + str(param_ub) + '\n')
	t.close()

def create_posterior_param_names():
	param_names = []
	param_names.append('a_prior')
	param_names.append('b_prior')
	param_names.append('c_prior')
	param_names.append('r_square')
	return param_names

def run_inference(inf_input_data_obj, output_stem, burn_in_iters=50, total_iters = 100):
	# Initialize model parameters
	betas, a_prior, b_prior, c_prior, sig_sq_priors, resid_vars, expression_resids, gene_names = initialize_data(inf_input_data_obj)

	print(len(gene_names))
	posterior_samples = []

	# Begin sampling
	for itera in range(total_iters):

		t1 = time.time()
		# E-step: compute posterior means and variances of beta under current hyperparameters.
		betas, posterior_beta_vars, expression_resids, resid_vars = update_causal_effect_estimates(betas, a_prior, sig_sq_priors, resid_vars, expression_resids, gene_names, inf_input_data_obj)

		print('update priors')
		# M-step: update a, b, and c using posterior second moments.
		a_prior, b_prior, c_prior, sig_sq_priors, r_squares = update_priors(betas, posterior_beta_vars, gene_names, inf_input_data_obj, b_prior, c_prior)
		t2 = time.time()

		print('###################', flush=True)
		print('Iteration ' + str(itera), flush=True)
		print(a_prior, flush=True)
		print(b_prior, flush=True)
		print(c_prior, flush=True)
		print(r_squares, flush=True)
		print(t2-t1, flush=True)

		if itera >= burn_in_iters:
			posterior_samples.append(np.hstack((np.asarray([a_prior, b_prior, c_prior]), r_squares)))

	if len(posterior_samples) > 0:
		param_names = create_posterior_param_names()
		save_posterior_samples(output_stem, np.vstack(posterior_samples), param_names)


########################
# Command line args
########################
borzoi_effect_file = sys.argv[1]
genotype_stem = sys.argv[2]
genotype_sample_mapping_file = sys.argv[3]
expr_file = sys.argv[4]
distribution = sys.argv[5]
output_stem = sys.argv[6]
standardize_geno = sys.argv[7]

###########################
# Load in data
###########################
# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

# Create mapping from gene id to expression vector
gene_id_to_expression_vector = create_mapping_from_gene_id_to_expression_vector(expr_file)


# Generate input data Object
inf_input_data_obj = generate_input_data_object_for_bayesian_inf(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, standardize_geno)
del gene_id_to_est_borzoi_effects
del genotype_sample_indices
del gene_id_to_expression_vector

'''
pickle_output_file = output_stem + '_inf_input_data_obj.pkl'
with open(pickle_output_file, 'wb') as f:
	pickle.dump(inf_input_data_obj, f)


pickle_output_file = output_stem + '_inf_input_data_obj.pkl'
with open(pickle_output_file, 'rb') as f:
	inf_input_data_obj = pickle.load(f)
'''

run_inference(inf_input_data_obj, output_stem, burn_in_iters=500, total_iters = 800)
