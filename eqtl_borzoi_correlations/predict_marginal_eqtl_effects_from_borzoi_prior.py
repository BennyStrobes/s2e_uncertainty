import gzip
import hashlib
import os
import sys

import numpy as np
from pandas_plink import read_plink


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
		effect = float(data[6])

		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			raise ValueError('Repeated variant-gene pair in Borzoi file: ' + gene_id + ', ' + var_id)
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
			raise ValueError('Repeated gene in expression file: ' + gene_id)
		dicti[gene_id] = expr
	f.close()
	return dicti


def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	mapping = {}
	for snp_iter, snp_name in enumerate(ordered_snps):
		if snp_name in mapping:
			raise ValueError('Repeated variant in genotype file: ' + snp_name)
		mapping[snp_name] = snp_iter
	return mapping


def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
	dicti = {}
	for ii, snp_id in enumerate(snp_array):
		if snp_id in dicti:
			raise ValueError('Repeated variant in genotype file: ' + snp_id)
		dicti[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
	return dicti


def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	return var_id_to_est_borzoi_effects[var_id][2]


def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_borzoi_effects):
	unique_vars = np.unique([*var_to_est_borzoi_effects])
	final_vars = []
	for var in unique_vars:
		if var not in rsid_to_genotype_index:
			continue
		geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])
		borzoi_alleles = var_to_est_borzoi_effects[var][4:6]
		if set(geno_alleles) != set(borzoi_alleles):
			continue
		final_vars.append(var)
	return np.asarray(final_vars)


def load_in_snp_gene_data(ordered_cis_variants, var_to_est_borzoi_effects):
	effects = []
	alleles = []
	for variant_id in ordered_cis_variants:
		var_info = var_to_est_borzoi_effects[variant_id]
		effects.append(var_info[6])
		alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles)


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
		raise ValueError('Could not find grid_knots in ' + summary_file)
	if model_effect_scale != expected_model_effect_scale:
		raise ValueError('Expected model_effect_scale=' + expected_model_effect_scale + ' in ' + summary_file)

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
		raise ValueError('Could not find a_prior in ' + param_file)

	ordered_coef_indices = np.sort(np.asarray([*grid_coef_mapping]))
	if np.array_equal(ordered_coef_indices, np.arange(len(grid_knots) + 2)) == False:
		raise ValueError('Grid coefficient indices do not match expected basis size in ' + param_file)

	return {
		'distribution_name': distribution_name,
		'model_effect_scale': model_effect_scale,
		'a_prior': a_prior,
		'grid_knots': grid_knots,
		'grid_coefs': np.asarray([grid_coef_mapping[index] for index in ordered_coef_indices])
	}


def load_in_ldscore_grid_distribution_data(borzoi_based_prior_output_stem):
	return load_in_ldscore_grid_distribution_data_from_files(borzoi_based_prior_output_stem, 'ldscore_grid', 'allelic_grid')


def load_in_ldscore_grid_squared_distribution_data(borzoi_based_prior_output_stem):
	return load_in_ldscore_grid_distribution_data_from_files(borzoi_based_prior_output_stem, 'ldscore_grid_squared', 'allelic_grid_squared')


def infer_ldscore_grid_distribution(borzoi_based_prior_output_stem):
	has_ldscore_grid = os.path.exists(borzoi_based_prior_output_stem + '_ldscore_grid_summary.txt')
	has_ldscore_grid_squared = os.path.exists(borzoi_based_prior_output_stem + '_ldscore_grid_squared_summary.txt')
	if has_ldscore_grid and has_ldscore_grid_squared:
		raise ValueError('Both ldscore_grid and ldscore_grid_squared prior files exist for this stem; pass distribution explicitly')
	if has_ldscore_grid:
		return 'ldscore_grid'
	if has_ldscore_grid_squared:
		return 'ldscore_grid_squared'
	raise ValueError('Could not find ldscore_grid or ldscore_grid_squared summary files for stem: ' + borzoi_based_prior_output_stem)


def compute_borzoi_prior_variance_allelic(borzoi_effects_unstandardized, distribution_obj):
	if distribution_obj['model_effect_scale'] == 'allelic_grid':
		borzoi_grid_basis = create_continuous_piecewise_linear_basis(np.abs(borzoi_effects_unstandardized), distribution_obj['grid_knots'])
	elif distribution_obj['model_effect_scale'] == 'allelic_grid_squared':
		borzoi_grid_basis = create_continuous_piecewise_linear_squared_basis(borzoi_effects_unstandardized, distribution_obj['grid_knots'])
	else:
		raise ValueError('Unsupported model_effect_scale: ' + str(distribution_obj['model_effect_scale']))
	return np.dot(borzoi_grid_basis, distribution_obj['grid_coefs'])


def compute_marginal_eqtl_effects(genotype_mat, expr_vec):
	n_snps = genotype_mat.shape[0]
	n_samples = genotype_mat.shape[1]
	y = expr_vec - np.mean(expr_vec)
	x_std = genotype_mat/np.std(genotype_mat, axis=1)[:, None]
	xtx = np.sum(np.square(x_std), axis=1)
	beta = np.dot(x_std, y)/xtx
	resid = y[None, :] - (beta[:, None]*x_std)
	resid_var = np.sum(np.square(resid), axis=1)/(n_samples - 2.0)
	beta_se_sq = resid_var/xtx
	return beta, beta_se_sq


def greedy_ld_prune(LD, gene_id, r_sq_threshold, random_seed):
	keep = np.zeros(LD.shape[0], dtype=bool)
	kept_indices = []
	gene_seed = (int(random_seed) + int(hashlib.md5(gene_id.encode('utf-8')).hexdigest()[:8], 16)) % (2**32)
	rng = np.random.default_rng(gene_seed)
	random_order = rng.permutation(LD.shape[0])
	for var_index in random_order:
		if len(kept_indices) == 0:
			keep[var_index] = True
			kept_indices.append(var_index)
			continue
		if np.all(np.square(LD[var_index, kept_indices]) <= r_sq_threshold):
			keep[var_index] = True
			kept_indices.append(var_index)
	return keep


def write_variant_row(t, gene_id, variant_id, chrom_num, snp_pos, geno_a0, geno_a1, borzoi_a0, borzoi_a1, snp_sdev, observed_beta_std, observed_beta_se_sq_std, borzoi_ld_predicted_mean_std, borzoi_ld_predicted_variance_std, borzoi_effects_unstandardized, borzoi_prior_mean_allelic, borzoi_prior_var_allelic, n_samples):
	t.write(
		gene_id + '\t' +
		variant_id + '\t' +
		str(chrom_num) + '\t' +
		snp_pos + '\t' +
		geno_a0 + '\t' +
		geno_a1 + '\t' +
		borzoi_a0 + '\t' +
		borzoi_a1 + '\t' +
		str(snp_sdev) + '\t' +
		str(observed_beta_std) + '\t' +
		str(observed_beta_se_sq_std) + '\t' +
		str(borzoi_ld_predicted_mean_std) + '\t' +
		str(borzoi_ld_predicted_variance_std) + '\t' +
		str(borzoi_effects_unstandardized) + '\t' +
		str(borzoi_prior_mean_allelic) + '\t' +
		str(borzoi_prior_var_allelic) + '\t' +
		str(n_samples) + '\n'
	)


def run_marginal_effect_prediction(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, distribution_obj, output_file, ld_pruned_output_file, ld_prune_r_sq_threshold, ld_prune_random_seed):
	output_dir = os.path.dirname(output_file)
	if output_dir != '':
		os.makedirs(output_dir, exist_ok=True)
	pruned_output_dir = os.path.dirname(ld_pruned_output_file)
	if pruned_output_dir != '':
		os.makedirs(pruned_output_dir, exist_ok=True)
	header = 'gene_id\tvariant_id\tchrom\tsnp_pos\tgeno_a0\tgeno_a1\tborzoi_a0\tborzoi_a1\tsnp_sdev\tobserved_marginal_effect_std\tobserved_marginal_effect_se_sq_std\tborzoi_ld_predicted_mean_std\tborzoi_ld_predicted_variance_std\tborzoi_effect_size_allelic\tborzoi_prior_mean_allelic\tborzoi_prior_variance_allelic\tn_samples\n'
	t = gzip.open(output_file, 'wt')
	t.write(header)
	t_pruned = gzip.open(ld_pruned_output_file, 'wt')
	t_pruned.write(header)

	for chrom_num in range(1,23):
		print(chrom_num, flush=True)
		(bim, fam, G) = read_plink(genotype_stem + str(chrom_num))
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))

		for gene_id in [*gene_id_to_est_borzoi_effects]:
			gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
			if str(gene_chrom_num) != str(chrom_num):
				continue
			if gene_id not in gene_id_to_expression_vector:
				continue

			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])
			if len(ordered_cis_variants) == 0:
				continue

			borzoi_effects_unstandardized, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])

			cis_genotype_indices = []
			geno_a0 = []
			geno_a1 = []
			snp_pos = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				geno_a0.append(str(snp_info[0]))
				geno_a1.append(str(snp_info[1]))
				snp_pos.append(str(snp_info[3]))

				if borzoi_variant_alleles[var_index, :][0] == geno_alleles[1]:
					borzoi_effects_unstandardized[var_index] = -1.0*borzoi_effects_unstandardized[var_index]

			cis_genotype_indices = np.asarray(cis_genotype_indices)
			geno_mat = (G[cis_genotype_indices, :].compute())[:, genotype_sample_indices]
			row_means = np.nanmean(geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(geno_mat))
			geno_mat[nan_rows, nan_cols] = row_means[nan_rows]

			snp_means = np.mean(geno_mat, axis=1)
			snp_sdevs = np.std(geno_mat, axis=1)
			valid_snps = np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
			if np.sum(valid_snps) == 0:
				continue

			geno_mat = geno_mat[valid_snps, :]
			geno_mat = geno_mat - snp_means[valid_snps, None]
			snp_sdevs = snp_sdevs[valid_snps]
			borzoi_effects_unstandardized = borzoi_effects_unstandardized[valid_snps]
			borzoi_variant_alleles = borzoi_variant_alleles[valid_snps, :]
			ordered_cis_variants_valid = ordered_cis_variants[valid_snps]
			geno_a0_valid = np.asarray(geno_a0)[valid_snps]
			geno_a1_valid = np.asarray(geno_a1)[valid_snps]
			snp_pos_valid = np.asarray(snp_pos)[valid_snps]

			expr_vec = gene_id_to_expression_vector[gene_id]
			if len(expr_vec) != geno_mat.shape[1]:
				raise ValueError('Expression and genotype sample sizes do not match for gene ' + gene_id)
			if len(expr_vec) <= 2:
				continue

			observed_beta_std, observed_beta_se_sq_std = compute_marginal_eqtl_effects(geno_mat, expr_vec)
			geno_mat_std = geno_mat/snp_sdevs[:, None]
			if geno_mat_std.shape[0] == 1:
				LD = np.ones((1,1))
			else:
				LD = np.corrcoef(geno_mat_std)

			borzoi_prior_mean_allelic = distribution_obj['a_prior']*borzoi_effects_unstandardized
			borzoi_prior_var_allelic = compute_borzoi_prior_variance_allelic(borzoi_effects_unstandardized, distribution_obj)
			if np.any(borzoi_prior_var_allelic < 0.0):
				print('warning: predicted negative per-allele prior variances; flooring at 0.0', flush=True)
			borzoi_prior_var_allelic = np.maximum(borzoi_prior_var_allelic, 0.0)

			borzoi_prior_mean_std = snp_sdevs*borzoi_prior_mean_allelic
			borzoi_prior_var_std = np.square(snp_sdevs)*borzoi_prior_var_allelic
			borzoi_ld_predicted_mean_std = np.dot(LD, borzoi_prior_mean_std)
			borzoi_ld_predicted_variance_std = np.dot(np.square(LD), borzoi_prior_var_std)
			ld_pruned_snps = greedy_ld_prune(LD, gene_id, ld_prune_r_sq_threshold, ld_prune_random_seed)

			for var_index, variant_id in enumerate(ordered_cis_variants_valid):
				row_args = (
					gene_id,
					variant_id,
					chrom_num,
					snp_pos_valid[var_index],
					geno_a0_valid[var_index],
					geno_a1_valid[var_index],
					borzoi_variant_alleles[var_index, 0],
					borzoi_variant_alleles[var_index, 1],
					snp_sdevs[var_index],
					observed_beta_std[var_index],
					observed_beta_se_sq_std[var_index],
					borzoi_ld_predicted_mean_std[var_index],
					borzoi_ld_predicted_variance_std[var_index],
					borzoi_effects_unstandardized[var_index],
					borzoi_prior_mean_allelic[var_index],
					borzoi_prior_var_allelic[var_index],
					len(expr_vec)
				)
				write_variant_row(t, *row_args)
				if ld_pruned_snps[var_index]:
					write_variant_row(t_pruned, *row_args)
			t.flush()
			t_pruned.flush()
	t.close()
	t_pruned.close()


borzoi_effect_file = sys.argv[1]
genotype_stem = sys.argv[2]
genotype_sample_mapping_file = sys.argv[3]
expr_file = sys.argv[4]
borzoi_based_prior_output_stem = sys.argv[5]
output_file = sys.argv[6]
ld_pruned_output_file = sys.argv[7]
ld_prune_r_sq_threshold = float(sys.argv[8])
ld_prune_random_seed = int(sys.argv[9])
distribution = 'auto'
if len(sys.argv) > 10:
	distribution = sys.argv[10]
if distribution == 'auto':
	distribution = infer_ldscore_grid_distribution(borzoi_based_prior_output_stem)

if distribution == 'ldscore_grid':
	distribution_obj = load_in_ldscore_grid_distribution_data(borzoi_based_prior_output_stem)
elif distribution == 'ldscore_grid_squared':
	distribution_obj = load_in_ldscore_grid_squared_distribution_data(borzoi_based_prior_output_stem)
else:
	raise ValueError('Unsupported distribution: ' + distribution)
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)
genotype_sample_indices = np.atleast_1d(np.loadtxt(genotype_sample_mapping_file)).astype(int)
gene_id_to_expression_vector = create_mapping_from_gene_id_to_expression_vector(expr_file)

run_marginal_effect_prediction(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, distribution_obj, output_file, ld_pruned_output_file, ld_prune_r_sq_threshold, ld_prune_random_seed)
