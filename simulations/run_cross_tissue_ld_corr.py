import gzip
import numpy as np
import sys
from pandas_plink import read_plink


def load_est_eqtl_effects(est_eqtl_effect_size_file):
	mapping = {}
	with gzip.open(est_eqtl_effect_size_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			gene_id = data[0]
			var_id = data[1]
			if gene_id not in mapping:
				mapping[gene_id] = {}
			if var_id in mapping[gene_id]:
				raise ValueError('Repeat SNP in eQTL file for ' + gene_id + '\t' + var_id)
			mapping[gene_id][var_id] = (
				gene_id,
				var_id,
				data[2],
				data[3],
				data[4],
				data[5],
				float(data[6]),
				float(data[7])
			)
	return mapping


def load_borzoi_effects(est_borzoi_effect_size_file):
	mapping = {}
	with gzip.open(est_borzoi_effect_size_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			gene_id = data[0]
			var_id = data[1]
			if gene_id not in mapping:
				mapping[gene_id] = {}
			if var_id in mapping[gene_id]:
				raise ValueError('Repeat SNP in Borzoi file for ' + gene_id + '\t' + var_id)
			mapping[gene_id][var_id] = (
				gene_id,
				var_id,
				data[2],
				data[3],
				data[4],
				data[5],
				float(data[6])
			)
	return mapping


def create_delta_borzoi_effects(gene_id_to_t1_borzoi, gene_id_to_t2_borzoi):
	mapping = {}
	for gene_id in gene_id_to_t1_borzoi:
		if gene_id not in gene_id_to_t2_borzoi:
			continue
		for var_id in gene_id_to_t1_borzoi[gene_id]:
			if var_id not in gene_id_to_t2_borzoi[gene_id]:
				continue
			t1_info = gene_id_to_t1_borzoi[gene_id][var_id]
			t2_info = gene_id_to_t2_borzoi[gene_id][var_id]
			if t1_info[:6] != t2_info[:6]:
				raise ValueError('Tissue Borzoi files are not aligned for ' + gene_id + '\t' + var_id)
			if gene_id not in mapping:
				mapping[gene_id] = {}
			mapping[gene_id][var_id] = t1_info[:6] + (t1_info[6] - t2_info[6],)
	return mapping


def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	return {snp_name: snp_iter for snp_iter, snp_name in enumerate(ordered_snps)}


def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
	return {snp_id: (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii]) for ii, snp_id in enumerate(snp_array)}


def extract_gene_chrom_num(var_id_to_effects):
	var_id = list(var_id_to_effects.keys())[0]
	return var_id_to_effects[var_id][2]


def extract_ordered_variants(rsid_to_genotype_index, rsid_to_snp_info, var_to_delta_borzoi, var_to_t1_eqtl, var_to_t2_eqtl):
	unique_vars = np.unique(np.hstack((list(var_to_delta_borzoi.keys()), list(var_to_t1_eqtl.keys()), list(var_to_t2_eqtl.keys()))))
	final_vars = []
	for var_id in unique_vars:
		if var_id not in rsid_to_genotype_index:
			continue
		geno_alleles = set(rsid_to_snp_info[var_id][:2])
		passing = True
		if var_id in var_to_delta_borzoi and set(var_to_delta_borzoi[var_id][4:6]) != geno_alleles:
			passing = False
		if var_id in var_to_t1_eqtl and set(var_to_t1_eqtl[var_id][4:6]) != geno_alleles:
			passing = False
		if var_id in var_to_t2_eqtl and set(var_to_t2_eqtl[var_id][4:6]) != geno_alleles:
			passing = False
		if passing:
			final_vars.append(var_id)
	return np.asarray(final_vars)


def compute_ld(geno_mat):
	row_means = np.nanmean(geno_mat, axis=1)
	nan_rows, nan_cols = np.where(np.isnan(geno_mat))
	geno_mat[nan_rows, nan_cols] = row_means[nan_rows]
	centered = geno_mat - np.mean(geno_mat, axis=1, keepdims=True)
	sdevs = np.std(geno_mat, axis=1, ddof=1)
	std_geno = np.zeros(geno_mat.shape)
	valid = (sdevs > 0.0) & np.isfinite(sdevs)
	std_geno[valid, :] = centered[valid, :]/sdevs[valid, None]
	return (std_geno @ std_geno.T)/(geno_mat.shape[1] - 1.0)


def load_delta_borzoi_for_ordered_variants(ordered_variants, var_to_delta_borzoi, rsid_to_snp_info):
	effects = []
	alleles = []
	for var_id in ordered_variants:
		if var_id not in var_to_delta_borzoi:
			effects.append(np.nan)
			alleles.append(('nan', 'nan'))
			continue
		var_info = var_to_delta_borzoi[var_id]
		effect = var_info[6]
		geno_alleles = rsid_to_snp_info[var_id][:2]
		if var_info[4] == geno_alleles[1]:
			effect = -1.0*effect
		effects.append(effect)
		alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles)


def load_delta_eqtl_for_ordered_variants(ordered_variants, var_to_t1_eqtl, var_to_t2_eqtl, rsid_to_snp_info, genotype_sdev_t1, genotype_sdev_t2):
	effects = []
	effect_ses = []
	for var_iter, var_id in enumerate(ordered_variants):
		if var_id not in var_to_t1_eqtl or var_id not in var_to_t2_eqtl:
			effects.append(np.nan)
			effect_ses.append(np.nan)
			continue
		t1_info = var_to_t1_eqtl[var_id]
		t2_info = var_to_t2_eqtl[var_id]
		if t1_info[3] != t2_info[3] or t1_info[4] != t2_info[4] or t1_info[5] != t2_info[5]:
			raise ValueError('Tissue eQTL files are not aligned for ' + t1_info[0] + '\t' + var_id)

		geno_alleles = rsid_to_snp_info[var_id][:2]
		t1_effect = t1_info[6]
		t1_se = t1_info[7]
		t2_effect = t2_info[6]
		t2_se = t2_info[7]
		if t1_info[4] == geno_alleles[1]:
			t1_effect = -1.0*t1_effect
		if t2_info[4] == geno_alleles[1]:
			t2_effect = -1.0*t2_effect

		t1_std_effect = t1_effect*genotype_sdev_t1[var_iter]
		t1_std_se = t1_se*genotype_sdev_t1[var_iter]
		t2_std_effect = t2_effect*genotype_sdev_t2[var_iter]
		t2_std_se = t2_se*genotype_sdev_t2[var_iter]
		effects.append(t1_std_effect - t2_std_effect)
		effect_ses.append(np.sqrt(np.square(t1_std_se) + np.square(t2_std_se)))
	return np.asarray(effects), np.asarray(effect_ses)


def generate_gene_ld_means(gene_id_to_delta_borzoi, gene_id_to_t1_eqtl, gene_id_to_t2_eqtl, onek_genomes_plink_filestem, t1_sample_indices, t2_sample_indices):
	gene_to_ld_means = {}
	ld_sample_indices = np.arange(max(np.max(t1_sample_indices), np.max(t2_sample_indices)) + 1)
	for chrom_num in range(1, 23):
		print(chrom_num)
		bim, fam, G = read_plink(onek_genomes_plink_filestem + str(chrom_num))
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))

		for gene_id in gene_id_to_delta_borzoi:
			if gene_id not in gene_id_to_t1_eqtl or gene_id not in gene_id_to_t2_eqtl:
				continue
			if str(extract_gene_chrom_num(gene_id_to_delta_borzoi[gene_id])) != str(chrom_num):
				continue

			ordered_variants = extract_ordered_variants(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_delta_borzoi[gene_id], gene_id_to_t1_eqtl[gene_id], gene_id_to_t2_eqtl[gene_id])
			if len(ordered_variants) < 10:
				continue

			cis_genotype_indices = np.asarray([rsid_to_genotype_index[var_id] for var_id in ordered_variants])
			full_geno_mat = G[cis_genotype_indices, :].compute()
			geno_mat_for_ld = full_geno_mat[:, ld_sample_indices]
			geno_mat_t1 = full_geno_mat[:, t1_sample_indices]
			geno_mat_t2 = full_geno_mat[:, t2_sample_indices]
			genotype_sdev_t1 = np.std(geno_mat_t1, axis=1, ddof=1)
			genotype_sdev_t2 = np.std(geno_mat_t2, axis=1, ddof=1)

			eqtl_effects, eqtl_ses = load_delta_eqtl_for_ordered_variants(ordered_variants, gene_id_to_t1_eqtl[gene_id], gene_id_to_t2_eqtl[gene_id], rsid_to_snp_info, genotype_sdev_t1, genotype_sdev_t2)
			borzoi_effects, _ = load_delta_borzoi_for_ordered_variants(ordered_variants, gene_id_to_delta_borzoi[gene_id], rsid_to_snp_info)
			LD = compute_ld(geno_mat_for_ld)

			observed_eqtl = np.isnan(eqtl_effects) == False
			eqtl_effects = eqtl_effects[observed_eqtl]
			eqtl_ses = eqtl_ses[observed_eqtl]
			LD = LD[observed_eqtl, :]

			observed_borzoi = np.isnan(borzoi_effects) == False
			borzoi_effects = borzoi_effects[observed_borzoi]
			LD = LD[:, observed_borzoi]
			if len(eqtl_effects) == 0 or len(borzoi_effects) == 0:
				continue

			n_samples = geno_mat_for_ld.shape[1]
			adj_sq_ld = np.square(LD) - ((1.0 - np.square(LD))/(n_samples - 2.0))
			gene_to_ld_means[gene_id] = {
				'eQTL_effect_sizes': eqtl_effects,
				'eQTL_effect_ses': eqtl_ses,
				'borzoi_effects': borzoi_effects,
				'ld_means': LD @ borzoi_effects,
				'ld_scores': adj_sq_ld @ np.ones(len(borzoi_effects))
			}
	return gene_to_ld_means


def compute_calibration_slope(gene_id_to_ld_means, ordered_gene_names):
	y = []
	x = []
	for gene_name in ordered_gene_names:
		y.append(gene_id_to_ld_means[gene_name]['eQTL_effect_sizes'])
		x.append(gene_id_to_ld_means[gene_name]['ld_means'])
	y = np.hstack(y)
	x = np.hstack(x)
	valid = np.isfinite(y) & np.isfinite(x)
	if np.sum(valid) == 0 or np.sum(np.square(x[valid])) == 0.0:
		return np.nan
	return np.sum(x[valid]*y[valid])/np.sum(np.square(x[valid]))


def compute_borzoi_variance(gene_id_to_ld_means, ordered_gene_names):
	borzoi = []
	for gene_name in ordered_gene_names:
		borzoi.append(gene_id_to_ld_means[gene_name]['borzoi_effects'])
	borzoi = np.hstack(borzoi)
	valid = np.isfinite(borzoi)
	if np.sum(valid) == 0:
		return np.nan
	return np.mean(np.square(borzoi[valid]))


def compute_eqtl_variance(gene_id_to_ld_means, ordered_gene_names):
	y = []
	x = []
	for gene_name in ordered_gene_names:
		y.append(np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_sizes']) - np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_ses']))
		x.append(gene_id_to_ld_means[gene_name]['ld_scores'])
	y = np.hstack(y)
	x = np.hstack(x)
	valid = np.isfinite(y) & np.isfinite(x)
	if np.sum(valid) < 2:
		return np.nan
	xx = np.vstack((np.ones(np.sum(valid)), x[valid])).T
	coefs, _, _, _ = np.linalg.lstsq(xx, y[valid], rcond=None)
	return coefs[1]


def run_ld_corr(ordered_gene_names, gene_id_to_ld_means):
	calibration_slope = compute_calibration_slope(gene_id_to_ld_means, ordered_gene_names)
	per_snp_eqtl_h2 = compute_eqtl_variance(gene_id_to_ld_means, ordered_gene_names)
	borzoi_variance = compute_borzoi_variance(gene_id_to_ld_means, ordered_gene_names)
	return calibration_slope, per_snp_eqtl_h2, borzoi_variance


def compute_correlation(calibration_slope, per_snp_eqtl_h2, borzoi_variance):
	sqrt_term = borzoi_variance/per_snp_eqtl_h2
	if np.isfinite(sqrt_term) and sqrt_term >= 0.0:
		return calibration_slope*np.sqrt(sqrt_term)
	return np.nan


t1_est_borzoi_effect_size_file = sys.argv[1]
t2_est_borzoi_effect_size_file = sys.argv[2]
t1_est_eqtl_effect_size_file = sys.argv[3]
t2_est_eqtl_effect_size_file = sys.argv[4]
onek_genomes_plink_filestem = sys.argv[5]
t1_eqtl_sample_size = int(sys.argv[6])
t2_eqtl_sample_size = int(sys.argv[7])
ld_corr_output_stem = sys.argv[8]
n_bootstraps = int(sys.argv[9]) if len(sys.argv) > 9 else 100


gene_id_to_t1_eqtl = load_est_eqtl_effects(t1_est_eqtl_effect_size_file)
gene_id_to_t2_eqtl = load_est_eqtl_effects(t2_est_eqtl_effect_size_file)
gene_id_to_t1_borzoi = load_borzoi_effects(t1_est_borzoi_effect_size_file)
gene_id_to_t2_borzoi = load_borzoi_effects(t2_est_borzoi_effect_size_file)
gene_id_to_delta_borzoi = create_delta_borzoi_effects(gene_id_to_t1_borzoi, gene_id_to_t2_borzoi)

t1_sample_indices = np.arange(t1_eqtl_sample_size)
t2_sample_indices = np.arange(t2_eqtl_sample_size)
gene_id_to_ld_means = generate_gene_ld_means(gene_id_to_delta_borzoi, gene_id_to_t1_eqtl, gene_id_to_t2_eqtl, onek_genomes_plink_filestem, t1_sample_indices, t2_sample_indices)

ordered_gene_names = np.asarray(list(gene_id_to_ld_means.keys()))
if len(ordered_gene_names) == 0:
	raise ValueError('No genes were available for cross-tissue LD-corr inference.')

calibration_slope, per_snp_eqtl_h2, borzoi_variance = run_ld_corr(ordered_gene_names, gene_id_to_ld_means)
observed_corr = compute_correlation(calibration_slope, per_snp_eqtl_h2, borzoi_variance)

bs_calibration_slopes = []
bs_per_snp_eqtl_h2 = []
bs_borzoi_variances = []
bs_corrs = []
for bs_iter in range(n_bootstraps):
	print(bs_iter)
	bs_gene_names = np.random.choice(ordered_gene_names, size=len(ordered_gene_names), replace=True)
	bs_calibration_slope, tmp_bs_per_snp_eqtl_h2, bs_borzoi_variance = run_ld_corr(bs_gene_names, gene_id_to_ld_means)
	bs_calibration_slopes.append(bs_calibration_slope)
	bs_per_snp_eqtl_h2.append(tmp_bs_per_snp_eqtl_h2)
	bs_borzoi_variances.append(bs_borzoi_variance)
	bs_corrs.append(compute_correlation(bs_calibration_slope, tmp_bs_per_snp_eqtl_h2, bs_borzoi_variance))

bs_calibration_slopes = np.asarray(bs_calibration_slopes)
bs_per_snp_eqtl_h2 = np.asarray(bs_per_snp_eqtl_h2)
bs_borzoi_variances = np.asarray(bs_borzoi_variances)
bs_corrs = np.asarray(bs_corrs)

results_output_file = ld_corr_output_stem + '_bootstrap_summary.txt'
with open(results_output_file, 'w') as t:
	t.write('sample_label\tsample_iter\tannotation_name\tcalibration_slope\tper_snp_eqtl_h2\tborzoi_variance\tcorrelation\n')
	t.write('observed\t-1\tall_variant_gene_pairs\t' + str(calibration_slope) + '\t' + str(per_snp_eqtl_h2) + '\t' + str(borzoi_variance) + '\t' + str(observed_corr) + '\n')
	for bs_iter in range(n_bootstraps):
		t.write('bootstrap\t' + str(bs_iter) + '\tall_variant_gene_pairs\t' + str(bs_calibration_slopes[bs_iter]) + '\t' + str(bs_per_snp_eqtl_h2[bs_iter]) + '\t' + str(bs_borzoi_variances[bs_iter]) + '\t' + str(bs_corrs[bs_iter]) + '\n')

summary_stats_output_file = ld_corr_output_stem + '_bootstrap_stats.txt'
with open(summary_stats_output_file, 'w') as t:
	t.write('annotation_name\toutput_name\tmean\tbootstrapped_mean\tbootstrap_se\tgaussian_z_score\tempirical_ci_lower\tempirical_ci_upper\n')
	for output_name, observed_value, bootstrap_distribution in [
		('calibration_slope', calibration_slope, bs_calibration_slopes),
		('per_snp_eqtl_h2', per_snp_eqtl_h2, bs_per_snp_eqtl_h2),
		('borzoi_variance', borzoi_variance, bs_borzoi_variances),
		('correlation', observed_corr, bs_corrs)
	]:
		bootstrap_mean = np.nanmean(bootstrap_distribution)
		bootstrap_se = np.nanstd(bootstrap_distribution, ddof=1)
		if bootstrap_se == 0.0 or np.isnan(bootstrap_se):
			gaussian_z_score = np.nan
		else:
			gaussian_z_score = observed_value/bootstrap_se
		empirical_ci = np.nanquantile(bootstrap_distribution, [.025, .975])
		t.write('all_variant_gene_pairs\t' + output_name + '\t' + str(observed_value) + '\t' + str(bootstrap_mean) + '\t' + str(bootstrap_se) + '\t' + str(gaussian_z_score) + '\t' + str(empirical_ci[0]) + '\t' + str(empirical_ci[1]) + '\n')

print(results_output_file)
print(summary_stats_output_file)
