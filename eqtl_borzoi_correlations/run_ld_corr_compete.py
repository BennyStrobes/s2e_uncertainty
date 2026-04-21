import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time
from sparse_multivariate_sumstat_regression import SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES
from single_causal_tissue_sumstat_regression import SINGLE_CAUSAL_TISSUE_SUMSTAT_REGRESSION




def extract_borzoi_track_names_and_effect_files(borzoi_standardized_tracks_file):
	track_names = []
	effect_files = []
	f = open(borzoi_standardized_tracks_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		track_names.append(data[0])
		effect_files.append(data[1])
	f.close()
	return np.asarray(track_names), np.asarray(effect_files)






def create_mapping_from_gene_id_to_causal_effects(borzoi_standardized_tracks_file, vg_pair_to_maf, borzoi_thresh):
	track_names, effect_files = extract_borzoi_track_names_and_effect_files(borzoi_standardized_tracks_file)
	if len(effect_files) == 0:
		print('assumption eroroor')
		pdb.set_trace()

	file_handles = []
	for effect_file in effect_files:
		file_handles.append(gzip.open(effect_file, 'rt'))

	for file_handle in file_handles:
		head_line = file_handle.readline().rstrip().split('\t')
		if len(head_line) < 7:
			print('assumption eroroor')
			pdb.set_trace()

	mapping = {}
	for line in file_handles[0]:
		line = line.rstrip()
		data = line.split('\t')
		gene_id = data[0]
		var_id = data[1]
		chrom_num = data[2]
		snp_pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		effects = np.zeros(len(file_handles))
		effects[0] = float(data[6])
		for file_iter in range(1, len(file_handles)):
			line2 = file_handles[file_iter].readline().rstrip()
			data2 = line2.split('\t')
			if np.array_equal(np.asarray(data[:6]), np.asarray(data2[:6])) == False:
				print('assumption eroroor')
				pdb.set_trace()
			effects[file_iter] = float(data2[6])
		
		vg_pair = var_id + ':' + gene_id
		if vg_pair not in vg_pair_to_maf:
			continue

		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()

		maf = vg_pair_to_maf[vg_pair]

		std_effects = effects*np.sqrt(2.0*maf*(1.0-maf))

		std_effects[np.abs(effects) < borzoi_thresh] = 0.0

		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, std_effects)
	for file_handle in file_handles:
		file_handle.close()
	return mapping

def create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file):	
	variant_gene_to_maf = {}
	f = gzip.open(est_eqtl_effect_size_file,'rt')
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
		pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		effect = float(data[6])
		se = float(data[7])
		eqtl_sample_size = float(data[8])
		maf = float(data[9])
		variant_gene_to_maf[var_id + ':' + gene_id] = maf
		desired_se = np.sqrt(1.0/eqtl_sample_size)
		if gene_id not in mapping:
			mapping[gene_id] = {}
		zed = effect/se
		if np.square(zed) > 80:
			continue
		std_effect = zed*desired_se
		if var_id in mapping[gene_id]:
			print('repeat snp error')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, pos, a0, a1, std_effect, desired_se)
	f.close()
	return mapping, variant_gene_to_maf



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

def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_borzoi_effects, var_to_est_eqtl_effects):
	unique_vars = np.unique(np.hstack(([*var_to_est_borzoi_effects],[*var_to_est_eqtl_effects])))
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
		if var in var_to_est_eqtl_effects:
			eqtl_alleles = var_to_est_eqtl_effects[var][4:6]
			if set(geno_alleles) != set(eqtl_alleles):
				passing = False

		if passing == False:
			continue
		final_vars.append(var)
	return np.asarray(final_vars)


def load_in_snp_gene_eqtl_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []
	effect_ses = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.nan)
			effect_ses.append(np.nan)
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			effect_ses.append(var_info[7])
			alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles), np.asarray(effect_ses)

def load_in_snp_gene_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []
	n_tracks = len(var_to_est_eqtl_effects[[*var_to_est_eqtl_effects][0]][6])

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.full(n_tracks, np.nan))
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.vstack(effects), np.asarray(alleles)
def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num

def generate_gene_ld_means(gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, onek_genomes_plink_filestem):
	# Initialize output object
	gene_to_ld_means = {}
	
	# Loop through chromsomes
	for chrom_num in range(1,23):
		print(chrom_num)

		##################################
		# Load in per-chrom-genotype data
		##################################
		# string of chromosome name
		chrom_string = 'chr' + str(chrom_num)
		# Load in chromosome plink data
		(bim, fam, G) = read_plink(onek_genomes_plink_filestem + str(chrom_num))
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

			# Gene needs both borzoi effects AND eQTLs
			if gene_id not in gene_id_to_est_eqtl_effects:
				continue

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id], gene_id_to_est_eqtl_effects[gene_id])
			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 10:
				continue

			# Load in data for gene
			# eQTL
			eqtl_effects, eqtl_variant_alleles, eqtl_effect_ses = load_in_snp_gene_eqtl_data(ordered_cis_variants, gene_id_to_est_eqtl_effects[gene_id])
			# Borzoi
			borzoi_effects, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])

			# Load in LD
			cis_genotype_indices = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				
				# Also flip signs of eqtls to match LD
				if np.isnan(eqtl_effects[var_index]) == False:
					if eqtl_variant_alleles[var_index,:][0] == geno_alleles[1]:
						eqtl_effects[var_index] = -1.0*eqtl_effects[var_index]
				if np.isnan(borzoi_effects[var_index,0]) == False:
					if borzoi_variant_alleles[var_index,:][0] == geno_alleles[1]:
						borzoi_effects[var_index,:] = -1.0*borzoi_effects[var_index,:]

			
			# Extract genotype
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = G[cis_genotype_indices,:].compute()
			row_means = np.nanmean(geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(geno_mat))
			geno_mat[nan_rows, nan_cols] = row_means[nan_rows]
			LD = np.corrcoef(geno_mat)

			# Subset LD by missingness
			# A. on eQTL end
			observed_eqtl_indices = np.isnan(eqtl_effects) == False
			eqtl_effects = eqtl_effects[observed_eqtl_indices]
			eqtl_effect_ses = eqtl_effect_ses[observed_eqtl_indices]
			LD = LD[observed_eqtl_indices, :]
			# B. on borzoi end
			observed_borzoi_indices = np.isnan(borzoi_effects[:,0]) == False
			borzoi_effects = borzoi_effects[observed_borzoi_indices,:]
			LD = LD[:, observed_borzoi_indices]

			# Add to dictionary
			if gene_id in gene_to_ld_means:
				print('assumption erororo')
				pdb.set_trace()

			gene_to_ld_means[gene_id] = {}
			gene_to_ld_means[gene_id]['eQTL_effect_sizes'] = eqtl_effects
			gene_to_ld_means[gene_id]['eQTL_effect_ses'] = eqtl_effect_ses
			gene_to_ld_means[gene_id]['borzoi_effects'] = borzoi_effects
			gene_to_ld_means[gene_id]['ld_means'] = LD @ (borzoi_effects)
			gene_to_ld_means[gene_id]['variant_anno'] = np.ones(borzoi_effects.shape)
			gene_to_ld_means[gene_id]['ld_scores'] = (LD**2) @ np.ones((LD.shape[1],1))


	return gene_to_ld_means


def compute_calibration_coefs(gene_id_to_ld_means, ordered_gene_names, joint_regression=False):
	yy = []
	xx = []
	annos = []
	for gene_name in ordered_gene_names:
		yy.append(gene_id_to_ld_means[gene_name]['eQTL_effect_sizes'])
		xx.append(gene_id_to_ld_means[gene_name]['ld_means'])
		annos.append(gene_id_to_ld_means[gene_name]['variant_anno'])
	yy = np.hstack(yy)
	xx = np.vstack(xx)
	annos = np.vstack(annos)


	if joint_regression:
		calibration_coefs, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)

	else:
		calibration_coefs = []
		for col_iter in range(xx.shape[1]):
			tmp_coef, _, _, _ = np.linalg.lstsq(xx[:, col_iter:(col_iter+1)], yy, rcond=None)
			calibration_coefs.append(tmp_coef[0])
		calibration_coefs = np.asarray(calibration_coefs)


	return calibration_coefs


def compute_borzoi_variances(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	for gene_name in ordered_gene_names:
		yy.append(np.square(gene_id_to_ld_means[gene_name]['borzoi_effects']))
		xx.append(gene_id_to_ld_means[gene_name]['variant_anno'])
	yy = np.vstack(yy)
	xx = np.vstack(xx)

	'''
	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]
	'''

	avg_pred_vars = []
	for anno_iter in range(xx.shape[1]):
		borzoi_var_coefs, _, _, _ = np.linalg.lstsq(xx[:,anno_iter:(anno_iter+1)], yy[:,anno_iter], rcond=None)
	
		pred_vars = xx[:,anno_iter]*borzoi_var_coefs[0]
		
		averager = np.sum(xx[:, anno_iter]*pred_vars)/np.sum(xx[:, anno_iter])
		avg_pred_vars.append(averager)

	return np.asarray(avg_pred_vars)

def compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	for gene_name in ordered_gene_names:
		yy.append(np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_sizes']) - np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_ses']))
		xx.append(gene_id_to_ld_means[gene_name]['ld_scores'])

	yy = np.hstack(yy)
	xx = np.vstack(xx)

	'''
	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]
	'''

	taus, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	return taus[0]


def run_ld_corr(ordered_gene_names, gene_id_to_ld_means, model):

	if model == 'joint':
		avg_calibration_slopes = compute_calibration_coefs(gene_id_to_ld_means, ordered_gene_names, joint_regression=True)
	elif model == 'independent':
		avg_calibration_slopes = compute_calibration_coefs(gene_id_to_ld_means, ordered_gene_names, joint_regression=False)
	else:
		print('not yet implemented' + model)
		pdb.set_trace()

	avg_borzoi_variances = compute_borzoi_variances(gene_id_to_ld_means, ordered_gene_names)
	avg_per_snp_h2 = compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names)
	return avg_calibration_slopes, avg_per_snp_h2, avg_borzoi_variances




def reload_in_bootstrapped_results(results_output_file):
	arr = []
	f = open(results_output_file)
	ordered_tiss = []
	head_count = 0
	mean_arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[1] == '-1':
			ordered_tiss.append(data[2])
			mean_arr.append(float(data[3]))
			continue
		arr.append(float(data[3]))
	f.close()
	ordered_tiss = np.asarray(ordered_tiss)
	mat = np.asarray(arr).reshape((len(arr) // len(ordered_tiss), len(ordered_tiss)))
	
	return mat, ordered_tiss, np.asarray(mean_arr)


def save_sparse_mv_regression_results(output_stem, ordered_track_names, sparse_mv_model):
	per_track_output_file = output_stem + '_sparse_mv_regression_per_track.txt'
	with open(per_track_output_file, 'w') as t:
		t.write('track_name\tposterior_mean\tposterior_variance\tcomponent_variance\n')
		for track_iter, track_name in enumerate(ordered_track_names):
			t.write(
				str(track_name) + '\t' +
				str(sparse_mv_model.beta_mu[track_iter]) + '\t' +
				str(sparse_mv_model.beta_cov[track_iter, track_iter]) + '\t' +
				str(sparse_mv_model.component_variances[track_iter]) + '\n'
			)

	cov_output_file = output_stem + '_sparse_mv_regression_posterior_cov.txt'
	np.savetxt(cov_output_file, sparse_mv_model.beta_cov, fmt='%.18e', delimiter='\t')

	convergence_output_file = output_stem + '_sparse_mv_regression_convergence.txt'
	with open(convergence_output_file, 'w') as t:
		t.write('iteration\tmean_abs_diff\n')
		for iter_num, diff in enumerate(sparse_mv_model.convergence_tracker):
			t.write(str(iter_num) + '\t' + str(diff) + '\n')

	return per_track_output_file, cov_output_file, convergence_output_file


def save_single_causal_tissue_regression_results(output_stem, ordered_track_names, single_causal_model):
	output_file = output_stem + '_single_causal_tissue_regression.txt'
	with open(output_file, 'w') as t:
		t.write('track_name\tis_winner\tposterior_mean\tposterior_variance\tcomponent_variance\tobjective\n')
		for track_iter, track_name in enumerate(ordered_track_names):
			t.write(
				str(track_name) + '\t' +
				str(int(track_iter == single_causal_model.best_tissue_index)) + '\t' +
				str(single_causal_model.beta_mu[track_iter]) + '\t' +
				str(single_causal_model.beta_cov[track_iter, track_iter]) + '\t' +
				str(single_causal_model.component_variances[track_iter]) + '\t' +
				str(single_causal_model.objective_per_tissue[track_iter]) + '\n'
			)
	return output_file



##########################
# Command line args
##########################
borzoi_standardized_tracks_file = sys.argv[1]
est_eqtl_effect_size_file = sys.argv[2]
onek_genomes_plink_filestem = sys.argv[3]
ld_corr_output_stem = sys.argv[4]
model = sys.argv[5]
borzoi_thresh = float(sys.argv[6])

##############################
# Load in data
##############################
# Create mapping from gene id to vector of est (standardized) eqtl effect sizes
gene_id_to_est_eqtl_effects, vg_pair_to_maf = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file)

# Load in track names
track_names = np.loadtxt(borzoi_standardized_tracks_file, dtype=str, delimiter='\t')[1:,0]
# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_standardized_tracks_file, vg_pair_to_maf, borzoi_thresh)


# Generate per gene ld-means
gene_id_to_ld_means = generate_gene_ld_means(gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, onek_genomes_plink_filestem)
del gene_id_to_est_eqtl_effects
del gene_id_to_est_borzoi_effects



##############################
# Run regressions
##############################

# Run standard regression
ordered_gene_names = np.asarray([*gene_id_to_ld_means])
avg_calibration_slopes, avg_per_snp_eqtl_h2, avg_borzoi_variances = run_ld_corr(ordered_gene_names, gene_id_to_ld_means, model)

# Bootstrap estimates!
n_bs = 100
bs_calibration_slopes = []
bs_per_snp_eqtl_h2 = []
bs_borzoi_variances = []
for bs_iter in range(n_bs):
	bs_genes = np.random.choice(ordered_gene_names, size=len(ordered_gene_names), replace=True)
	tmp_bs_calibration_slopes, tmp_bs_per_snp_eqtl_h2, tmp_bs_borzoi_variances = run_ld_corr(bs_genes, gene_id_to_ld_means, model)

	bs_calibration_slopes.append(tmp_bs_calibration_slopes)
	bs_per_snp_eqtl_h2.append(tmp_bs_per_snp_eqtl_h2)
	bs_borzoi_variances.append(tmp_bs_borzoi_variances)

bs_calibration_slopes = np.asarray(bs_calibration_slopes)
bs_per_snp_eqtl_h2 = np.asarray(bs_per_snp_eqtl_h2)
bs_borzoi_variances = np.asarray(bs_borzoi_variances)




#############################
# Print results to output file
##############################

results_output_file = ld_corr_output_stem + '_bootstrap_summary.txt'
with open(results_output_file, 'w') as t:
	t.write('sample_label\tsample_iter\ttrack_name\tcalibration_slope\tborzoi_variance\n')
	for track_iter, track_name in enumerate(track_names):
		t.write('observed\t-1\t' + str(track_name) + '\t' + str(avg_calibration_slopes[track_iter]) + '\t' + str(avg_borzoi_variances[track_iter]) + '\n')

	for bs_iter in range(n_bs):
		for track_iter, track_name in enumerate(track_names):
			t.write('bootstrap\t' + str(bs_iter) + '\t' + str(track_name) + '\t' + str(bs_calibration_slopes[bs_iter, track_iter]) + '\t' + str(bs_borzoi_variances[bs_iter, track_iter]) + '\n')

eqtl_var_summary_output_file = ld_corr_output_stem + '_eqtl_variance_bootstrap_summary.txt'
with open(eqtl_var_summary_output_file, 'w') as t:
	t.write('sample_label\tsample_iter\tper_snp_eqtl_h2\n')
	t.write('observed\t-1\t' + str(avg_per_snp_eqtl_h2) + '\n')
	for bs_iter in range(n_bs):
		t.write('bootstrap\t' + str(bs_iter) + '\t' + str(bs_per_snp_eqtl_h2[bs_iter]) + '\n')

summary_stats_output_file = ld_corr_output_stem + '_bootstrap_stats.txt'
with open(summary_stats_output_file, 'w') as t:
	t.write('track_name\toutput_name\tmean\tbootstrapped_mean\tbootstrap_se\tgaussian_z_score\tempirical_ci_lower\tempirical_ci_upper\n')
	output_names = ['calibration_slope', 'borzoi_variance']
	observed_mat = np.vstack((avg_calibration_slopes, avg_borzoi_variances)).T
	bootstrap_mats = [bs_calibration_slopes, bs_borzoi_variances]
	for track_iter, track_name in enumerate(track_names):
		for output_iter, output_name in enumerate(output_names):
			bootstrap_distribution = bootstrap_mats[output_iter][:, track_iter]
			bootstrap_mean = np.mean(bootstrap_distribution)
			bootstrap_se = np.std(bootstrap_distribution, ddof=1)
			if bootstrap_se == 0.0:
				gaussian_z_score = np.nan
			else:
				gaussian_z_score = observed_mat[track_iter, output_iter]/bootstrap_se
			empirical_ci = np.quantile(bootstrap_distribution, [.025, .975])
			t.write(str(track_name) + '\t' + output_name + '\t' + str(observed_mat[track_iter, output_iter]) + '\t' + str(bootstrap_mean) + '\t' + str(bootstrap_se) + '\t' + str(gaussian_z_score) + '\t' + str(empirical_ci[0]) + '\t' + str(empirical_ci[1]) + '\n')

eqtl_var_stats_output_file = ld_corr_output_stem + '_eqtl_variance_bootstrap_stats.txt'
with open(eqtl_var_stats_output_file, 'w') as t:
	t.write('output_name\tmean\tbootstrapped_mean\tbootstrap_se\tgaussian_z_score\tempirical_ci_lower\tempirical_ci_upper\n')
	bootstrap_mean = np.mean(bs_per_snp_eqtl_h2)
	bootstrap_se = np.std(bs_per_snp_eqtl_h2, ddof=1)
	if bootstrap_se == 0.0:
		gaussian_z_score = np.nan
	else:
		gaussian_z_score = avg_per_snp_eqtl_h2/bootstrap_se
	empirical_ci = np.quantile(bs_per_snp_eqtl_h2, [.025, .975])
	t.write('per_snp_eqtl_h2\t' + str(avg_per_snp_eqtl_h2) + '\t' + str(bootstrap_mean) + '\t' + str(bootstrap_se) + '\t' + str(gaussian_z_score) + '\t' + str(empirical_ci[0]) + '\t' + str(empirical_ci[1]) + '\n')







'''
results_output_file = ld_corr_output_stem + '_bootstrap_summary.txt'

bs_calibration_slopes, ordered_track_names, bs_mean_arr = reload_in_bootstrapped_results(results_output_file)

calibration_slope_mean = np.mean(bs_calibration_slopes,axis=0)
calibration_slope_cov = np.cov(np.transpose(bs_calibration_slopes))

sparse_mv_model = SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES()
sparse_mv_model.fit(bs_mean_arr, calibration_slope_cov)

print(sparse_mv_model.beta_mu)
print(ordered_track_names)



single_causal_model = SINGLE_CAUSAL_TISSUE_SUMSTAT_REGRESSION()
single_causal_model.fit(bs_mean_arr, calibration_slope_cov)
pdb.set_trace()
single_causal_output_file = save_single_causal_tissue_regression_results(ld_corr_output_stem, ordered_track_names, single_causal_model)
print(single_causal_output_file)
'''
