import gzip
import sys
import pdb
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
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, effect)
	f.close()
	return mapping


def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	mapping = {}
	for snp_iter, snp_name in enumerate(ordered_snps):
		mapping[snp_name] = snp_iter
	return mapping


def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
	mapping = {}
	for ii, snp_id in enumerate(snp_array):
		mapping[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
	return mapping


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


def create_continuous_piecewise_linear_squared_basis(delta, knots):
	delta_sq = np.square(np.asarray(delta))
	basis = [np.ones(len(delta_sq)), delta_sq]
	for knot in knots:
		basis.append(np.maximum(delta_sq - knot, 0.0))
	return np.transpose(np.vstack(basis))


def compute_ldscore_grid_squared_prior_variance(borzoi_vec, distribution_obj):
	grid_basis = create_continuous_piecewise_linear_squared_basis(borzoi_vec, distribution_obj['grid_knots'])
	return np.dot(grid_basis, distribution_obj['grid_coefs'])


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

	ordered_coef_indices = np.sort(np.asarray([*grid_coef_mapping]))
	grid_coefs = np.asarray([grid_coef_mapping[index] for index in ordered_coef_indices])
	return {'a_prior': a_prior, 'grid_knots': grid_knots, 'grid_coefs': grid_coefs}


def extract_specific_example_data(borzoi_effect_file, genotype_stem, gwas_ss_file, borzoi_based_prior_output_stem, gene_id, output_file):
	gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_effect_file)
	if gene_id not in gene_id_to_est_borzoi_effects:
		print('Could not find gene ' + gene_id + ' in ' + borzoi_effect_file)
		sys.exit(1)

	gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
	bim, fam, G = read_plink(genotype_stem + str(gene_chrom_num))
	rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
	rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))
	ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])

	borzoi_effects, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])
	for var_index, cis_variant in enumerate(ordered_cis_variants):
		geno_alleles = rsid_to_snp_info[cis_variant][:2]
		if borzoi_variant_alleles[var_index, :][0] == geno_alleles[0]:
			borzoi_effects[var_index] = -1.0*borzoi_effects[var_index]

	variant_to_gwas_sumstats = create_mapping_from_variant_id_to_gwas_ss(gwas_ss_file)
	gwas_effect_sizes = []
	gwas_ses = []
	for variant_id in ordered_cis_variants:
		if variant_id in variant_to_gwas_sumstats:
			gwas_effect_sizes.append(variant_to_gwas_sumstats[variant_id][0])
			gwas_ses.append(variant_to_gwas_sumstats[variant_id][1])
		else:
			gwas_effect_sizes.append(np.nan)
			gwas_ses.append(np.nan)
	gwas_effect_sizes = np.asarray(gwas_effect_sizes)
	gwas_ses = np.asarray(gwas_ses)
	gwas_z_scores = gwas_effect_sizes/gwas_ses

	distribution_obj = load_in_ldscore_grid_squared_distribution_data(borzoi_based_prior_output_stem)
	borzoi_variances = compute_ldscore_grid_squared_prior_variance(borzoi_effects, distribution_obj)
	borzoi_variances_floored = np.maximum(borzoi_variances, 0.0)

	t = open(output_file,'w')
	t.write('gene_id\tvariant_id\tchrom\tposition\tgenotype_a0\tgenotype_a1\tgwas_beta\tgwas_beta_se\tgwas_z\tborzoi_effect_aligned_to_genotype\tborzoi_prior_variance\tborzoi_prior_variance_floored\n')
	for var_index, variant_id in enumerate(ordered_cis_variants):
		snp_info = rsid_to_snp_info[variant_id]
		t.write(gene_id + '\t' + variant_id + '\t' + str(snp_info[2]) + '\t' + str(snp_info[3]) + '\t' + str(snp_info[0]) + '\t' + str(snp_info[1]) + '\t' + str(gwas_effect_sizes[var_index]) + '\t' + str(gwas_ses[var_index]) + '\t' + str(gwas_z_scores[var_index]) + '\t' + str(borzoi_effects[var_index]) + '\t' + str(borzoi_variances[var_index]) + '\t' + str(borzoi_variances_floored[var_index]) + '\n')
	t.close()
	print(output_file)


########################
# Command line args
########################
borzoi_effect_file = sys.argv[1]
genotype_stem = sys.argv[2]
gwas_ss_file = sys.argv[3]
borzoi_based_prior_output_stem = sys.argv[4]
gene_id = sys.argv[5]
output_file = sys.argv[6]

extract_specific_example_data(borzoi_effect_file, genotype_stem, gwas_ss_file, borzoi_based_prior_output_stem, gene_id, output_file)
