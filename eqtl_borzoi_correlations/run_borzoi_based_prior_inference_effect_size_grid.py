import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time

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


def generate_input_data_object_for_bayesian_inf(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem):
	# Initialize output object
	gene_to_data = {}
	
	# Loop through chromsomes
	for chrom_num in range(8,23):
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
			if gene_id not in gene_id_to_variant_gene_anno:
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
						print('need to flix flip. should flip genotypes not snps')
						pdb.set_trace()
				if borzoi_variant_alleles[var_index,:][0] != borzoi_anno_variant_alleles[var_index,:][0]:
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
			geno_mat = (geno_mat - snp_means[:, None])/snp_sdevs[:, None]
			borzoi_effects_standardized = borzoi_effects_unstandardized*snp_sdevs

			gene_to_data[gene_id] = {}
			gene_to_data[gene_id]['expression'] = gene_id_to_expression_vector[gene_id]
			gene_to_data[gene_id]['anno'] = variant_anno[valid_snps, :]
			gene_to_data[gene_id]['borzoi'] = borzoi_effects_standardized[valid_snps]
			gene_to_data[gene_id]['genotype'] = np.ascontiguousarray(geno_mat[valid_snps, :])


	return gene_to_data



def initialize_data(inf_input_data_obj, grid_size):
	gene_names = np.asarray([*inf_input_data_obj])
	n_bins = inf_input_data_obj[gene_names[0]]['anno'].shape[1]
	n_genes = len(gene_names)

	# Initialize priors
	effect_size_grid_pis = np.ones((n_bins, grid_size))/grid_size
	effect_size_grid_pis[:,4:9] = effect_size_grid_pis[:,4:9] + 1.0
	for row_iter in range(effect_size_grid_pis.shape[0]):
		effect_size_grid_pis[row_iter, :] = effect_size_grid_pis[row_iter, :]/np.sum(effect_size_grid_pis[row_iter, :])


	#prior_var_fixed_grid = np.logspace(-7, 0, num=grid_size)
	if grid_size != 13:
		print('assumption error: not yet implemented')
		pdb.set_trace()
	'''
	effect_size_grid = np.array([
		-0.5, -0.2, -0.075, -0.03, -.01, -0.005, 0.0,
		 0.005,.01, 0.03, 0.075, 0.2, 0.5
	])
	'''
	effect_size_grid = np.array([
		-0.5, -0.2, -0.075, -0.025, -0.005, -0.0001, 0.0,
		 0.0001, 0.005,.025, 0.075, 0.2, 0.5
	])

	# residual variance for each gene
	resid_vars = np.ones(n_genes)

	# Initialize causal effects
	betas = []
	valid_genes = []
	gene_annos = []
	latent_states = []

	expression_resids = []

	for gene_name in gene_names:

		borzoi_preds = inf_input_data_obj[gene_name]['borzoi']

		n_snps = len(borzoi_preds)

		if n_snps < 10:
			print('asusmption error: too few snps')
			pdb.set_trace()

		init_gene_causal_effects = np.zeros(len(borzoi_preds))
		betas.append(init_gene_causal_effects)
		latent_states.append(init_gene_causal_effects.astype(int))

		resid_E = np.copy(inf_input_data_obj[gene_name]['expression']) - np.dot(init_gene_causal_effects, inf_input_data_obj[gene_name]['genotype'])
		expression_resids.append(resid_E)

		anno_mat = inf_input_data_obj[gene_name]['anno']
		tmp_anno = []
		for row_iter in range(anno_mat.shape[0]):

			indices = np.where(anno_mat[row_iter, :] == 1)[0]
			if len(indices) != 1:
				print('assumptioneornronro')
				pdb.set_trace()
			tmp_anno.append(indices[0])

		gene_annos.append(np.asarray(tmp_anno))


	return betas, latent_states, effect_size_grid, effect_size_grid_pis, resid_vars, expression_resids, gene_annos, gene_names,


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

@njit(cache=True)
def sample_categorical_from_log_probs(log_probs):
	max_log_prob = np.max(log_probs)
	total = 0.0
	probs = np.empty(len(log_probs))
	for ii in range(len(log_probs)):
		probs[ii] = np.exp(log_probs[ii] - max_log_prob)
		total = total + probs[ii]

	u = np.random.random()*total
	cdf = 0.0
	for ii in range(len(probs)):
		cdf = cdf + probs[ii]
		if u <= cdf:
			return ii
	return len(probs) - 1

@njit(cache=True)
def update_causal_effect_for_single_gene_ashr_shell(cur_gene_beta, cur_gene_latent_state, prior_means, prior_var_pis, prior_var_fixed_grid, anno_vec, gene_resid_var, resid_expr_vec, genotype_mat):
	n_snps = len(cur_gene_beta)
	n_samples = len(resid_expr_vec)
	n_grid = len(prior_var_fixed_grid)
	log_probs = np.empty(n_grid)

	for snp_iter in np.random.permutation(n_snps):
		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] + (geno_val*cur_gene_beta[snp_iter])

		anno_bin = anno_vec[snp_iter]
		prior_mean = prior_means[snp_iter]
		diff = cur_gene_beta[snp_iter] - prior_mean
		for grid_iter in range(n_grid):
			prior_var = prior_var_fixed_grid[grid_iter]
			log_probs[grid_iter] = np.log(prior_var_pis[anno_bin, grid_iter] + 1e-300) - (0.5*np.log(2.0*np.pi*prior_var)) - (0.5*(diff*diff)/prior_var)

		new_latent_state = sample_categorical_from_log_probs(log_probs)
		cur_gene_latent_state[snp_iter] = new_latent_state

		prior_var = prior_var_fixed_grid[new_latent_state]
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

		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] - (geno_val*new_beta)
		cur_gene_beta[snp_iter] = new_beta

	return cur_gene_beta, cur_gene_latent_state, resid_expr_vec

@njit(cache=True)
def update_causal_effect_for_single_gene_effect_size_grid_shell(cur_gene_beta, cur_gene_latent_state, fixed_effect_size_grid, effect_size_grid_pis, gene_resid_var, resid_expr_vec, anno_vec, genotype_mat):
	n_snps = len(cur_gene_beta)
	n_samples = len(resid_expr_vec)
	n_grid = len(fixed_effect_size_grid)
	log_probs = np.empty(n_grid)

	for snp_iter in np.random.permutation(n_snps):
		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] + (geno_val*cur_gene_beta[snp_iter])

		anno_bin = anno_vec[snp_iter]
		sum_sq_geno = 0.0
		geno_resid_dot = 0.0
		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			sum_sq_geno = sum_sq_geno + (geno_val*geno_val)
			geno_resid_dot = geno_resid_dot + (geno_val*resid_expr_vec[sample_iter])

		for grid_iter in range(n_grid):
			grid_effect = fixed_effect_size_grid[grid_iter]
			log_probs[grid_iter] = np.log(effect_size_grid_pis[anno_bin, grid_iter] + 1e-300) - (.5/gene_resid_var)*((grid_effect*grid_effect)*sum_sq_geno - (2.0*grid_effect*geno_resid_dot))

		new_latent_state = sample_categorical_from_log_probs(log_probs)
		new_beta = fixed_effect_size_grid[new_latent_state]

		cur_gene_latent_state[snp_iter] = new_latent_state
		cur_gene_beta[snp_iter] = new_beta

		for sample_iter in range(n_samples):
			geno_val = genotype_mat[snp_iter, sample_iter]
			resid_expr_vec[sample_iter] = resid_expr_vec[sample_iter] - (geno_val*new_beta)

	return cur_gene_beta, cur_gene_latent_state, resid_expr_vec

def update_causal_effect_for_single_gene(cur_gene_beta, cur_gene_latent_state, fixed_effect_size_grid, effect_size_grid_pis, gene_resid_var, resid_expr_vec, anno_vec, genotype_mat):
	return update_causal_effect_for_single_gene_effect_size_grid_shell(cur_gene_beta, cur_gene_latent_state, fixed_effect_size_grid, effect_size_grid_pis, gene_resid_var, resid_expr_vec, anno_vec, genotype_mat)

def update_resid_var_for_single_gene(resid_expr_vec, v0=1e-3, s_sq0=1e-3):
	n_samples = len(resid_expr_vec)
	shape_post = v0 + (n_samples/2.0)
	scale_post = s_sq0 + (np.sum(np.square(resid_expr_vec))/2.0)
	return 1.0/np.random.gamma(shape_post, 1.0/scale_post)

def update_causal_effect_estimates(betas, latent_states, fixed_effect_size_grid, effect_size_grid_pis, resid_vars, expression_resids, gene_snp_annos, gene_names, inf_input_data_obj):
	# Loop through genes (perform update seperately for each gene)
	for gene_iter, gene_name in enumerate(gene_names):

		#####################
		# Load in data for this gene
		gene_beta, gene_latent_state, gene_expression_resid = update_causal_effect_for_single_gene(betas[gene_iter], latent_states[gene_iter], fixed_effect_size_grid, effect_size_grid_pis, resid_vars[gene_iter], expression_resids[gene_iter], gene_snp_annos[gene_iter], inf_input_data_obj[gene_name]['genotype'])
		betas[gene_iter] = gene_beta
		latent_states[gene_iter] = gene_latent_state
		expression_resids[gene_iter] = gene_expression_resid
		resid_vars[gene_iter] = update_resid_var_for_single_gene(gene_expression_resid)

	return betas, latent_states, expression_resids, resid_vars


def update_priors(latent_states, effect_size_grid_pis, fixed_effect_size_grid, gene_snp_annos):
	dirichlet_alpha = np.ones(len(fixed_effect_size_grid))*.1

	all_anno = np.hstack(gene_snp_annos)
	all_latent_states = np.hstack(latent_states)

	n_bins = effect_size_grid_pis.shape[0]
	n_grid = len(fixed_effect_size_grid)
	for bin_iter in range(n_bins):
		bin_latent_states = all_latent_states[all_anno == bin_iter]
		if len(bin_latent_states) == 0:
			effect_size_grid_pis[bin_iter, :] = np.random.dirichlet(dirichlet_alpha)
			continue

		counts = np.zeros(n_grid)
		for grid_iter in range(n_grid):
			counts[grid_iter] = np.sum(bin_latent_states == grid_iter)
		effect_size_grid_pis[bin_iter, :] = np.random.dirichlet(dirichlet_alpha + counts)

	return effect_size_grid_pis

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

def save_prior_effect_size_grid(output_stem, prior_var_fixed_grid):
	grid_output_file = output_stem + '_fixed_effect_size_grid.txt'
	t = open(grid_output_file, 'w')
	t.write('grid_index\teffect_size\n')
	for grid_iter, grid_val in enumerate(prior_var_fixed_grid):
		t.write(str(grid_iter) + '\t' + str(grid_val) + '\n')
	t.close()

def create_posterior_param_names(n_bins, n_grid):
	param_names = []
	for param_iter in range(n_bins):
		for grid_iter in range(n_grid):
			param_names.append('effect_size_grid_pi_' + str(param_iter) + '_' + str(grid_iter))
	return param_names

def run_inference(inf_input_data_obj, output_stem, burn_in_iters=50, total_iters = 100, grid_size=13):
	# Initialize model parameters
	betas, latent_states, fixed_effect_size_grid, effect_size_grid_pis, resid_vars, expression_resids, gene_snp_annos, gene_names = initialize_data(inf_input_data_obj, grid_size)

	save_prior_effect_size_grid(output_stem, fixed_effect_size_grid)


	print(len(gene_names))
	posterior_samples = []

	# Begin sampling
	for itera in range(total_iters):

		t1 = time.time()
		# Update causal effects
		betas, latent_states, expression_resids, resid_vars = update_causal_effect_estimates(betas, latent_states, fixed_effect_size_grid, effect_size_grid_pis, resid_vars, expression_resids, gene_snp_annos, gene_names, inf_input_data_obj)

		# Update priors
		if itera > 10:
			effect_size_grid_pis = update_priors(latent_states, effect_size_grid_pis, fixed_effect_size_grid, gene_snp_annos)
		t2 = time.time()

		print('###################', flush=True)
		print('Iteration ' + str(itera), flush=True)
		print(effect_size_grid_pis, flush=True)
		print(t2-t1, flush=True)


		if itera >= burn_in_iters:
			posterior_samples.append(np.hstack(effect_size_grid_pis))

	if len(posterior_samples) > 0:
		param_names = create_posterior_param_names(effect_size_grid_pis.shape[0], effect_size_grid_pis.shape[1])
		save_posterior_samples(output_stem, np.vstack(posterior_samples), param_names)


########################
# Command line args
########################
borzoi_effect_file = sys.argv[1]
borzoi_annotation_file = sys.argv[2]
genotype_stem = sys.argv[3]
genotype_sample_mapping_file = sys.argv[4]
expr_file = sys.argv[5]
distribution = sys.argv[6]
output_stem = sys.argv[7]

###########################
# Load in data
###########################
# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)

# Create mapping from gene id to vector of variant-gene annotations
gene_id_to_variant_gene_anno, anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(borzoi_annotation_file)

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

# Create mapping from gene id to expression vector
gene_id_to_expression_vector = create_mapping_from_gene_id_to_expression_vector(expr_file)


# Generate input data Object
inf_input_data_obj = generate_input_data_object_for_bayesian_inf(gene_id_to_est_borzoi_effects, gene_id_to_variant_gene_anno, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem)
del gene_id_to_est_borzoi_effects
del gene_id_to_variant_gene_anno
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

run_inference(inf_input_data_obj, output_stem, burn_in_iters=5000, total_iters = 5200, grid_size=13)


# CONSIDER: DIRICHLET PRIOR, VARIANCE PRIOR, letting it run a bit before updating priors
