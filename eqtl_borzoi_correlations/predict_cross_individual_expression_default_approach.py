import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time









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
			expr_sample_names = np.asarray(data[4:])
			continue
		gene_id = data[3].split('.')[0]
		expr = np.asarray(data[4:]).astype(float)
		if gene_id in dicti:
			print('assumption erororr')
			pdb.set_trace()
		dicti[gene_id] = expr
	f.close()
	return dicti, expr_sample_names

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

def compute_mean_causal_effects_from_distribution(borzoi_vec, anno, distribution_obj):
	if len(borzoi_vec) != anno.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	avg_mu = distribution_obj['avg_mu']
	bin_indices = distribution_obj['bin_indices']

	if anno.shape[1] != len(avg_mu):
		print('assumption eroror')
		pdb.set_trace()
	if np.array_equal(bin_indices, np.arange(len(avg_mu))) == False:
		print('assumption eroror')
		pdb.set_trace()

	snp_bins = np.argmax(anno, axis=1)
	if np.all(np.sum(anno, axis=1) == 1.0) == False:
		print('assumption eroror')
		pdb.set_trace()

	return avg_mu[snp_bins]*borzoi_vec

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


def run_expression_prediction(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, expr_sample_names, output_file, gene_individual_output_file):
	# Initialize output file
	t = open(output_file,'w')
	t.write('gene_id\traw_correlation\tregression_coef_output_obs\tregression_coef_output_pred\tmax_abs_borzoi\n')

	t2 = open(gene_individual_output_file,'w')
	t2.write('gene_id\tind_id\tobs_expr\tborzoi_pred_expr\n')


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


			expr_vec = gene_id_to_expression_vector[gene_id]
			borzoi_vec = borzoi_effects_standardized[valid_snps]
			borzoi_vec_unstandardized = borzoi_effects_unstandardized[valid_snps]
			genotype_mat = np.transpose(geno_mat[valid_snps, :])

			default_pred_expr = np.dot(genotype_mat, borzoi_vec)
			default_corry = np.corrcoef(expr_vec, default_pred_expr)[0,1]

			design_mat_pred = np.column_stack((np.ones(len(default_pred_expr)), default_pred_expr))
			regression_coef_output_obs = np.linalg.lstsq(design_mat_pred, expr_vec, rcond=None)[0][1]

			design_mat_expr = np.column_stack((np.ones(len(expr_vec)), expr_vec))
			regression_coef_output_pred = np.linalg.lstsq(design_mat_expr, default_pred_expr, rcond=None)[0][1]

			t.write(gene_id + '\t' + str(default_corry) + '\t' + str(regression_coef_output_obs) + '\t' + str(regression_coef_output_pred) + '\t' + str(np.max(np.abs(borzoi_vec_unstandardized))) + '\n')
			t.flush()

			for sample_iter, ind_id in enumerate(expr_sample_names):
				t2.write(gene_id + '\t' + ind_id + '\t' + str(expr_vec[sample_iter]) + '\t' + str(default_pred_expr[sample_iter]) + '\n')
			t2.flush()


	t.close()
	t2.close()

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



########################
# Command line args
########################
borzoi_effect_file = sys.argv[1]
genotype_stem = sys.argv[2]
genotype_sample_mapping_file = sys.argv[3]
expr_file = sys.argv[4]
output_file = sys.argv[5]
gene_individual_output_file = sys.argv[6]



###########################
# Load in data
###########################


# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)


# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

# Create mapping from gene id to expression vector
gene_id_to_expression_vector, expr_sample_names = create_mapping_from_gene_id_to_expression_vector(expr_file)


# run expression prediction
run_expression_prediction(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, genotype_stem, expr_sample_names, output_file, gene_individual_output_file)
