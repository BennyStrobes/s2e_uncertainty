import numpy as np
import sys
import pdb
import gzip
import pickle
import time
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize


def create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file):
	f = gzip.open(est_borzoi_effect_size_file, 'rt')
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
		if passing is False:
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


def generate_input_data_object_for_bayesian_inf(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, standardize_geno):
	from pandas_plink import read_plink
	gene_to_data = {}

	for chrom_num in range(19, 23):
		print(chrom_num)
		(bim, fam, G) = read_plink(genotype_stem + str(chrom_num))
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(
			np.asarray(bim['snp']),
			np.asarray(bim['a0']),
			np.asarray(bim['a1']),
			np.asarray(bim['chrom']),
			np.asarray(bim['pos'])
		)

		for gene_id in [*gene_id_to_est_borzoi_effects]:
			gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
			if str(gene_chrom_num) != str(chrom_num):
				continue
			if gene_id not in gene_id_to_expression_vector:
				continue

			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(
				rsid_to_genotype_index,
				rsid_to_snp_info,
				gene_id_to_est_borzoi_effects[gene_id]
			)
			if len(ordered_cis_variants) < 10:
				continue

			borzoi_effects_unstandardized, borzoi_variant_alleles = load_in_snp_gene_data(
				ordered_cis_variants,
				gene_id_to_est_borzoi_effects[gene_id]
			)

			cis_genotype_indices = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				if np.isnan(borzoi_effects_unstandardized[var_index]) is False:
					if borzoi_variant_alleles[var_index, :][0] == geno_alleles[0]:
						borzoi_effects_unstandardized[var_index] = -1.0*borzoi_effects_unstandardized[var_index]
						print('need to flix flip. should flip genotypes not snps')
						pdb.set_trace()

			cis_genotype_indices = np.asarray(cis_genotype_indices)
			geno_mat = (G[cis_genotype_indices, :].compute())[:, genotype_sample_indices]
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
			else:
				print('standardize_geno assumption error')
				pdb.set_trace()

			abs_borzoi = np.abs(borzoi_effects_for_model)
			anno_mat = np.hstack((
				np.ones((len(borzoi_effects_for_model), 1)),
				abs_borzoi.reshape((len(borzoi_effects_for_model), 1)),
				np.square(abs_borzoi).reshape((len(borzoi_effects_for_model), 1)),
				np.power(abs_borzoi, 3).reshape((len(borzoi_effects_for_model), 1)),
				np.power(abs_borzoi, 4).reshape((len(borzoi_effects_for_model), 1))
			))

			gene_to_data[gene_id] = {}
			gene_to_data[gene_id]['expression'] = gene_id_to_expression_vector[gene_id]
			gene_to_data[gene_id]['anno'] = anno_mat[valid_snps, :]
			gene_to_data[gene_id]['borzoi'] = borzoi_effects_for_model[valid_snps]
			gene_to_data[gene_id]['genotype'] = np.ascontiguousarray(geno_mat[valid_snps, :])

	if len(gene_to_data) > 0:
		ordered_gene_names = np.asarray([*gene_to_data])
		all_anno = np.vstack([gene_to_data[gene_name]['anno'] for gene_name in ordered_gene_names])
		if all_anno.shape[1] > 1:
			anno_means = np.mean(all_anno[:, 1:], axis=0)
			anno_sdevs = np.std(all_anno[:, 1:], axis=0)
			invalid_sdevs = np.isnan(anno_sdevs) | (anno_sdevs == 0.0)
			anno_sdevs[invalid_sdevs] = 1.0
			for gene_name in ordered_gene_names:
				gene_to_data[gene_name]['anno'][:, 1:] = (gene_to_data[gene_name]['anno'][:, 1:] - anno_means)/anno_sdevs

	return gene_to_data


def create_optimization_data(gene_to_data):
	opt_data = []
	for gene_name in [*gene_to_data]:
		geno = gene_to_data[gene_name]['genotype']
		expr = gene_to_data[gene_name]['expression']
		xtx = np.dot(geno, geno.T)
		xty = np.dot(geno, expr)
		yty = np.dot(expr, expr)
		expr_var = np.var(expr)
		if expr_var <= 1e-8:
			expr_var = 1.0
		opt_data.append({
			'gene_name': gene_name,
			'anno': gene_to_data[gene_name]['anno'],
			'borzoi': gene_to_data[gene_name]['borzoi'],
			'xtx': xtx,
			'xty': xty,
			'yty': yty,
			'n_samples': len(expr),
			'init_log_resid_var': np.log(expr_var)
		})
	return opt_data


def unpack_params(params, n_gamma, n_genes):
	mu = params[0]
	gamma = params[1:(1 + n_gamma)]
	log_resid_vars = params[(1 + n_gamma):(1 + n_gamma + n_genes)]
	return mu, gamma, log_resid_vars


def robust_cholesky(mat):
	jitter = 0.0
	for _ in range(8):
		try:
			if jitter == 0.0:
				return cho_factor(mat, lower=True, check_finite=False)
			return cho_factor(mat + (jitter*np.eye(mat.shape[0])), lower=True, check_finite=False)
		except np.linalg.LinAlgError:
			if jitter == 0.0:
				jitter = 1e-8
			else:
				jitter = jitter*10.0
	print('Cholesky decomposition failed')
	pdb.set_trace()


def compute_sorted_predicted_prior_vars(all_anno, gamma):
	return np.sort(np.exp(np.dot(all_anno, gamma)))


def negative_log_marginal_likelihood_and_grad(params, opt_data, mu_prior_var=100.0, gamma_prior_var=.1):
	n_gamma = opt_data[0]['anno'].shape[1]
	n_genes = len(opt_data)
	mu, gamma, log_resid_vars = unpack_params(params, n_gamma, n_genes)

	objective = 0.0
	grad_mu = 0.0
	grad_gamma = np.zeros(n_gamma)
	grad_log_resid = np.zeros(n_genes)

	for gene_iter, gene_data in enumerate(opt_data):
		anno = gene_data['anno']
		delta = gene_data['borzoi']
		xtx = gene_data['xtx']
		xty = gene_data['xty']
		yty = gene_data['yty']
		n_samples = gene_data['n_samples']

		eta = np.dot(anno, gamma)
		prior_var = np.exp(eta)
		prior_prec = np.exp(-eta)
		resid_var = np.exp(log_resid_vars[gene_iter])
		prior_mean = mu*delta

		b_mat = np.diag(prior_prec) + (xtx/resid_var)
		chol_b = robust_cholesky(b_mat)

		r_vec = xty - np.dot(xtx, prior_mean)
		z_vec = r_vec/resid_var
		rtr = yty - (2.0*np.dot(prior_mean, xty)) + np.dot(prior_mean, np.dot(xtx, prior_mean))
		solve_z = cho_solve(chol_b, z_vec, check_finite=False)
		logdet_b = 2.0*np.sum(np.log(np.diag(chol_b[0])))
		objective = objective + 0.5*(
			n_samples*np.log(2.0*np.pi*resid_var) +
			np.sum(eta) +
			logdet_b +
			(rtr/resid_var) -
			np.dot(z_vec, solve_z)
		)

		rhs = (prior_prec*prior_mean) + (xty/resid_var)
		posterior_mean = cho_solve(chol_b, rhs, check_finite=False)
		posterior_cov = cho_solve(chol_b, np.eye(len(delta)), check_finite=False)
		posterior_var = np.diag(posterior_cov)
		expected_sq_resid = posterior_var + np.square(posterior_mean - prior_mean)

		grad_mu = grad_mu + np.sum(delta*(prior_mean - posterior_mean)*prior_prec)
		grad_gamma = grad_gamma + 0.5*np.dot(anno.T, (1.0 - (expected_sq_resid*prior_prec)))

		expected_resid_ss = yty - (2.0*np.dot(posterior_mean, xty))
		expected_resid_ss = expected_resid_ss + np.dot(posterior_mean, np.dot(xtx, posterior_mean))
		expected_resid_ss = expected_resid_ss + np.sum(posterior_cov*xtx)
		grad_log_resid[gene_iter] = 0.5*(n_samples - (expected_resid_ss/resid_var))

	if np.isfinite(mu_prior_var):
		objective = objective + (0.5*np.square(mu)/mu_prior_var)
		grad_mu = grad_mu + (mu/mu_prior_var)
	if np.isfinite(gamma_prior_var):
		objective = objective + (0.5*np.sum(np.square(gamma[1:]))/gamma_prior_var)
		grad_gamma[1:] = grad_gamma[1:] + (gamma[1:]/gamma_prior_var)

	grad = np.hstack((np.asarray([grad_mu]), grad_gamma, grad_log_resid))
	return objective, grad


class ObjectiveWrapper(object):
	def __init__(self, opt_data, mu_prior_var=100.0, gamma_prior_var=100.0):
		self.opt_data = opt_data
		self.mu_prior_var = mu_prior_var
		self.gamma_prior_var = gamma_prior_var
		self.n_eval = 0
		self.last_objective = np.nan

	def __call__(self, params):
		self.n_eval = self.n_eval + 1
		objective, grad = negative_log_marginal_likelihood_and_grad(
			params,
			self.opt_data,
			mu_prior_var=self.mu_prior_var,
			gamma_prior_var=self.gamma_prior_var
		)
		self.last_objective = objective
		if self.n_eval == 1 or np.mod(self.n_eval, 5) == 0:
			print('Objective eval ' + str(self.n_eval) + ': ' + str(objective), flush=True)
		return objective, grad


def fit_hyperparameters(opt_data, output_stem, mu_prior_var=100.0, gamma_prior_var=100.0, maxiter=200):
	n_gamma = opt_data[0]['anno'].shape[1]
	n_genes = len(opt_data)
	all_anno = np.vstack([gene_data['anno'] for gene_data in opt_data])
	init_mu = 1.0
	init_gamma = np.zeros(n_gamma)
	init_log_resid_vars = np.asarray([gene_data['init_log_resid_var'] for gene_data in opt_data])
	init_params = np.hstack((np.asarray([init_mu]), init_gamma, init_log_resid_vars))

	objective_wrapper = ObjectiveWrapper(opt_data, mu_prior_var=mu_prior_var, gamma_prior_var=gamma_prior_var)
	bounds = [(None, None)]*(1 + n_gamma) + [(-10.0, 10.0)]*n_genes

	def callback(params):
		mu, gamma, log_resid_vars = unpack_params(params, n_gamma, n_genes)
		print('Iter callback: mu=' + str(mu) + ', objective=' + str(objective_wrapper.last_objective), flush=True)
		print(compute_sorted_predicted_prior_vars(all_anno, gamma), flush=True)

	result = minimize(
		fun=objective_wrapper,
		x0=init_params,
		method='L-BFGS-B',
		jac=True,
		bounds=bounds,
		callback=callback,
		options={'disp': True, 'maxiter': maxiter}
	)
	return result


def save_results(output_stem, result, opt_data):
	n_gamma = opt_data[0]['anno'].shape[1]
	n_genes = len(opt_data)
	mu, gamma, log_resid_vars = unpack_params(result.x, n_gamma, n_genes)

	summary_output_file = output_stem + '_marginalized_fit_summary.txt'
	t = open(summary_output_file, 'w')
	t.write('field\tvalue\n')
	t.write('success\t' + str(result.success) + '\n')
	t.write('status\t' + str(result.status) + '\n')
	t.write('message\t' + str(result.message).replace('\n', ' ') + '\n')
	t.write('objective\t' + str(result.fun) + '\n')
	t.write('n_iterations\t' + str(result.nit) + '\n')
	t.write('n_function_evals\t' + str(result.nfev) + '\n')
	t.write('mu\t' + str(mu) + '\n')
	for gamma_iter, gamma_val in enumerate(gamma):
		t.write('gamma_' + str(gamma_iter) + '\t' + str(gamma_val) + '\n')
	t.close()

	resid_output_file = output_stem + '_marginalized_gene_resid_vars.txt'
	t = open(resid_output_file, 'w')
	t.write('gene_id\tlog_resid_var\tresid_var\n')
	for gene_iter, gene_data in enumerate(opt_data):
		log_resid_var = log_resid_vars[gene_iter]
		t.write(
			gene_data['gene_name'] + '\t' +
			str(log_resid_var) + '\t' +
			str(np.exp(log_resid_var)) + '\n'
		)
	t.close()


def main():
	borzoi_effect_file = sys.argv[1]
	genotype_stem = sys.argv[2]
	genotype_sample_mapping_file = sys.argv[3]
	expr_file = sys.argv[4]
	distribution = sys.argv[5]
	output_stem = sys.argv[6]
	standardize_geno = sys.argv[7]
	if len(sys.argv) > 8:
		maxiter = int(sys.argv[8])
	else:
		maxiter = 200

	gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)
	genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)
	gene_id_to_expression_vector = create_mapping_from_gene_id_to_expression_vector(expr_file)

	inf_input_data_obj = generate_input_data_object_for_bayesian_inf(
		gene_id_to_est_borzoi_effects,
		genotype_sample_indices,
		gene_id_to_expression_vector,
		genotype_stem,
		standardize_geno
	)
	del gene_id_to_est_borzoi_effects
	del genotype_sample_indices
	del gene_id_to_expression_vector

	pickle_output_file = output_stem + '_marginalized_inf_input_data_obj.pkl'
	with open(pickle_output_file, 'wb') as f:
		pickle.dump(inf_input_data_obj, f)

	opt_data = create_optimization_data(inf_input_data_obj)
	del inf_input_data_obj
	del distribution

	t0 = time.time()
	result = fit_hyperparameters(opt_data, output_stem, maxiter=maxiter)
	t1 = time.time()
	print('Optimization runtime: ' + str(t1 - t0), flush=True)
	save_results(output_stem, result, opt_data)


if __name__ == "__main__":
	main()
