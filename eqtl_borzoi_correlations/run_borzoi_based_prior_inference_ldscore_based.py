import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time
from scipy.interpolate import BSpline
from scipy.optimize import minimize





def create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev, standardize=True):
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

		if var_id not in variant_id_to_genotype_sdev:
			continue

		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()

		geno_sdev = variant_id_to_genotype_sdev[var_id]
		if standardize:
			effect = effect*geno_sdev

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


def create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file):	
	variant_id_to_geno_sdev = {}
	f = gzip.open(est_eqtl_effect_size_file,'rt')
	mapping = {}
	head_count = 0
	bad_genes = {}
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
		if len(data) < 12:
			genotype_sdev = np.sqrt(2.0*maf*(1.0-maf))
		else:
			genotype_sdev = float(data[11])

		std_effect = effect*genotype_sdev
		std_se = se*genotype_sdev
		variant_id_to_geno_sdev[var_id] = genotype_sdev
		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('repeat snp error')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, pos, a0, a1, std_effect, std_se)
	f.close()
	return mapping, variant_id_to_geno_sdev

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

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.nan)
			alleles.append(('nan', 'nan'))
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
		else:
			var_info = var_to_variant_gene_anno[variant_id]
			annos.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.vstack(annos), np.asarray(alleles)

def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num

def generate_gene_ld_means(gene_id_to_est_borzoi_effects,gene_id_to_est_borzoi_effects_unstandardized, gene_id_to_est_eqtl_effects, onek_genomes_plink_filestem, genotype_sample_indices, variant_id_to_genotype_sdev, anno_method):
	# Initialize output object
	gene_to_ld_means = {}

	# Loop through chromsomes
	for chrom_num in range(18,23):
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
			if len(ordered_cis_variants) < 100:
				continue

			# Load in data for gene
			# eQTL
			eqtl_effects, eqtl_variant_alleles, eqtl_effect_ses = load_in_snp_gene_eqtl_data(ordered_cis_variants, gene_id_to_est_eqtl_effects[gene_id])
			# Borzoi
			borzoi_effects, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])
			borzoi_effects_unstandardized, borzoi_variant_alleles_unstandardized = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects_unstandardized[gene_id])


			# Load in LD
			cis_genotype_indices = []
			variant_sdevs = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				variant_sdevs.append(variant_id_to_genotype_sdev[cis_variant])
				
				# Also flip signs of eqtls to match LD
				if np.isnan(eqtl_effects[var_index]) == False:
					if eqtl_variant_alleles[var_index,:][0] == geno_alleles[1]:
						eqtl_effects[var_index] = -1.0*eqtl_effects[var_index]
				if np.isnan(borzoi_effects[var_index]) == False:
					if borzoi_variant_alleles[var_index,:][0] == geno_alleles[1]:
						borzoi_effects[var_index] = -1.0*borzoi_effects[var_index]
						borzoi_effects_unstandardized[var_index] = -1.0*borzoi_effects_unstandardized[var_index]
			variant_sdevs = np.asarray(variant_sdevs)

			# Extract genotype
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = (G[cis_genotype_indices,:].compute())[:, genotype_sample_indices]
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
			observed_borzoi_indices = np.isnan(borzoi_effects) == False
			borzoi_effects = borzoi_effects[observed_borzoi_indices]
			borzoi_effects_unstandardized = borzoi_effects_unstandardized[observed_borzoi_indices]
			LD = LD[:, observed_borzoi_indices]

			eqtl_sdevs = variant_sdevs[observed_eqtl_indices]
			borzoi_sdevs = variant_sdevs[observed_borzoi_indices]

			# Add to dictionary
			if gene_id in gene_to_ld_means:
				print('assumption erororo')
				pdb.set_trace()

			sq_LD = LD**2

			gene_anno = np.transpose(np.vstack((np.ones(len(borzoi_effects_unstandardized)), borzoi_effects_unstandardized)))


			gene_to_ld_means[gene_id] = {}
			gene_to_ld_means[gene_id]['eQTL_effect_sizes'] = eqtl_effects
			gene_to_ld_means[gene_id]['eQTL_effect_ses'] = eqtl_effect_ses
			gene_to_ld_means[gene_id]['eQTL_variant_sdevs'] = eqtl_sdevs
			gene_to_ld_means[gene_id]['borzoi_effects'] = borzoi_effects
			gene_to_ld_means[gene_id]['borzoi_variant_sdevs'] = borzoi_sdevs
			gene_to_ld_means[gene_id]['borzoi_effects_unstandardized'] = borzoi_effects_unstandardized
			gene_to_ld_means[gene_id]['ld_means'] = np.dot(LD, borzoi_effects)
			gene_to_ld_means[gene_id]['ld_score_intercept'] = np.sum(sq_LD, axis=1)
			gene_to_ld_means[gene_id]['ld_score_delta_sq'] = np.dot(sq_LD, np.square(borzoi_effects))
			gene_to_ld_means[gene_id]['ld_score_allelic_intercept'] = np.dot(sq_LD, np.square(borzoi_sdevs))
			gene_to_ld_means[gene_id]['ld_score_allelic_delta_sq'] = np.dot(sq_LD, np.square(borzoi_sdevs)*np.square(borzoi_effects_unstandardized))
			gene_to_ld_means[gene_id]['anno'] = gene_anno


	return gene_to_ld_means


def compute_calibration_coefs(gene_id_to_ld_means, ordered_gene_names):
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

	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]
	#annos = annos[valid_rows, :]


	calibration_coefs, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	
	pred_coefs = np.dot(annos, calibration_coefs)
	avg_pred_coefs = []
	for anno_iter in range(annos.shape[1]):
		averager = np.sum(annos[:, anno_iter]*pred_coefs)/np.sum(annos[:, anno_iter])
		avg_pred_coefs.append(averager)

	return np.asarray(avg_pred_coefs)


def compute_borzoi_variances(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	annos = []
	for gene_name in ordered_gene_names:
		yy.append(np.square(gene_id_to_ld_means[gene_name]['borzoi_effects']))
		xx.append(gene_id_to_ld_means[gene_name]['variant_anno'])
	yy = np.hstack(yy)
	xx = np.vstack(xx)

	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]

	borzoi_var_coefs, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	
	pred_vars = np.dot(xx, borzoi_var_coefs)
	avg_pred_vars = []
	for anno_iter in range(xx.shape[1]):

		averager = np.sum(xx[:, anno_iter]*pred_vars)/np.sum(xx[:, anno_iter])
		avg_pred_vars.append(averager)

	return np.asarray(avg_pred_vars)

def compute_unstandardized_borzoi_variances(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	annos = []
	for gene_name in ordered_gene_names:
		yy.append(np.square(gene_id_to_ld_means[gene_name]['borzoi_effects_unstandardized']))
		xx.append(gene_id_to_ld_means[gene_name]['variant_anno'])
	yy = np.hstack(yy)
	xx = np.vstack(xx)

	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]

	borzoi_var_coefs, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	
	pred_vars = np.dot(xx, borzoi_var_coefs)
	avg_pred_vars = []
	for anno_iter in range(xx.shape[1]):

		averager = np.sum(xx[:, anno_iter]*pred_vars)/np.sum(xx[:, anno_iter])
		avg_pred_vars.append(averager)

	return np.asarray(avg_pred_vars)


def compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	annos = []
	for gene_name in ordered_gene_names:
		yy.append(np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_sizes']) - np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_ses']))
		xx.append(gene_id_to_ld_means[gene_name]['ld_scores'])
		annos.append(gene_id_to_ld_means[gene_name]['variant_anno'])

	yy = np.hstack(yy)
	xx = np.vstack(xx)
	annos = np.vstack(annos)

	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]
	#annos = annos[valid_rows, :]


	taus, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	
	pred_per_snp_h2s = np.dot(annos, taus)
	avg_pred_h2s = []
	for anno_iter in range(annos.shape[1]):
		averager = np.sum(annos[:, anno_iter]*pred_per_snp_h2s)/np.sum(annos[:, anno_iter])
		avg_pred_h2s.append(averager)

	return np.asarray(avg_pred_h2s)


def run_ld_corr(ordered_gene_names, gene_id_to_ld_means):
	avg_calibration_slopes = compute_calibration_coefs(gene_id_to_ld_means, ordered_gene_names)
	avg_borzoi_variances = compute_borzoi_variances(gene_id_to_ld_means, ordered_gene_names)
	avg_unstandardized_borzoi_variances = compute_unstandardized_borzoi_variances(gene_id_to_ld_means, ordered_gene_names)
	avg_per_snp_h2 = compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names)

	return avg_calibration_slopes, avg_per_snp_h2, avg_borzoi_variances, avg_unstandardized_borzoi_variances




def compute_calibration_coefficient(gene_to_ld_means):
	yy = []
	xx = []
	annos = []
	for gene_name in [*gene_to_ld_means]:
		yy.append(gene_to_ld_means[gene_name]['eQTL_effect_sizes'])
		xx.append(gene_to_ld_means[gene_name]['ld_means'])
	yy = np.hstack(yy)
	xx = np.hstack(xx)

	valid_rows = np.isfinite(yy) & np.isfinite(xx)
	yy = yy[valid_rows]
	xx = xx[valid_rows]
	#annos = annos[valid_rows, :]

	coef = np.sum(xx * yy) / np.sum(xx * xx)	

	return coef


def compute_calibration_coefficient_bootstrap_var(gene_to_ld_means, n_bootstraps=500, random_seed=1):
	gene_names = np.asarray([*gene_to_ld_means])
	gene_numerators = []
	gene_denominators = []
	for gene_name in gene_names:
		yy = gene_to_ld_means[gene_name]['eQTL_effect_sizes']
		xx = gene_to_ld_means[gene_name]['ld_means']
		valid_rows = np.isfinite(yy) & np.isfinite(xx)
		yy = yy[valid_rows]
		xx = xx[valid_rows]
		gene_numerators.append(np.sum(xx*yy))
		gene_denominators.append(np.sum(xx*xx))

	gene_numerators = np.asarray(gene_numerators)
	gene_denominators = np.asarray(gene_denominators)
	calibration_coef = np.sum(gene_numerators)/np.sum(gene_denominators)

	if n_bootstraps <= 1 or len(gene_names) <= 1:
		return calibration_coef, 0.0, 0.0

	rng = np.random.default_rng(random_seed)
	bootstrap_coefs = []
	for bootstrap_iter in range(n_bootstraps):
		bs_indices = rng.integers(0, len(gene_names), size=len(gene_names))
		bs_denominator = np.sum(gene_denominators[bs_indices])
		if bs_denominator == 0.0:
			continue
		bs_numerator = np.sum(gene_numerators[bs_indices])
		bootstrap_coefs.append(bs_numerator/bs_denominator)

	if len(bootstrap_coefs) <= 1:
		return calibration_coef, 0.0, 0.0

	bootstrap_coefs = np.asarray(bootstrap_coefs)
	calibration_coef_var = np.var(bootstrap_coefs, ddof=1)
	calibration_coef_se = np.sqrt(calibration_coef_var)
	return calibration_coef, calibration_coef_var, calibration_coef_se


def create_bspline_basis_info(x, n_basis=6, degree=3, fixed_interior_knots=None):
	x = np.asarray(x)
	finite_x = x[np.isfinite(x)]
	if len(finite_x) == 0:
		print('assumption eroror: no finite spline inputs')
		pdb.set_trace()

	x_min = np.min(finite_x)
	x_max = np.max(finite_x)
	unique_x = np.unique(finite_x)

	if len(unique_x) == 1:
		return {
			'type': 'constant',
			'n_basis': 1,
			'degree': 0,
			'knots': np.asarray([x_min, x_max]),
			'x_min': x_min,
			'x_max': x_max
		}

	degree = min(degree, n_basis - 1, len(unique_x) - 1)
	if degree < 1:
		return {
			'type': 'constant',
			'n_basis': 1,
			'degree': 0,
			'knots': np.asarray([x_min, x_max]),
			'x_min': x_min,
			'x_max': x_max
		}

	if fixed_interior_knots is not None:
		lower_boundary = 0.0
		upper_boundary = x_max
		interior_knots = np.asarray(fixed_interior_knots).astype(float)
		interior_knots = np.unique(interior_knots[(interior_knots > lower_boundary) & (interior_knots < upper_boundary)])
	else:
		lower_boundary = x_min
		upper_boundary = x_max
		n_interior_target = max(0, n_basis - degree - 1)
		if n_interior_target == 0:
			interior_knots = np.asarray([])
		else:
			probs = np.linspace(0.0, 1.0, n_interior_target + 2)[1:-1]
			interior_knots = np.quantile(finite_x, probs)
			interior_knots = np.unique(interior_knots[(interior_knots > lower_boundary) & (interior_knots < upper_boundary)])

	knots = np.concatenate((np.repeat(lower_boundary, degree + 1), interior_knots, np.repeat(upper_boundary, degree + 1)))
	n_basis_actual = len(knots) - degree - 1
	if n_basis_actual < 1:
		return {
			'type': 'constant',
			'n_basis': 1,
			'degree': 0,
			'knots': np.asarray([lower_boundary, upper_boundary]),
			'x_min': lower_boundary,
			'x_max': x_max
		}

	return {
		'type': 'bspline',
		'n_basis': n_basis_actual,
		'degree': degree,
		'knots': knots,
		'x_min': lower_boundary,
		'x_max': x_max
	}


def evaluate_bspline_basis(x, basis_info):
	x = np.asarray(x)
	if basis_info['type'] == 'constant':
		return np.ones((len(x), 1))

	design_mat = np.zeros((len(x), basis_info['n_basis']))
	for basis_iter in range(basis_info['n_basis']):
		coef = np.zeros(basis_info['n_basis'])
		coef[basis_iter] = 1.0
		spline_obj = BSpline(basis_info['knots'], coef, basis_info['degree'], extrapolate=False)
		spline_vals = spline_obj(x)
		spline_vals[np.isnan(spline_vals)] = 0.0
		design_mat[:, basis_iter] = spline_vals

	row_sums = np.sum(design_mat, axis=1)
	valid_rows = row_sums > 0.0
	design_mat[valid_rows, :] = design_mat[valid_rows, :]/row_sums[valid_rows, None]
	return design_mat


def create_latent_design_mat_from_spline_basis(spline_basis_mat):
	spline_basis_mat = np.asarray(spline_basis_mat)
	n_rows, n_cols = spline_basis_mat.shape
	if n_cols <= 1:
		return np.ones((n_rows, 1))
	return np.hstack((np.ones((n_rows, 1)), spline_basis_mat[:, :-1]))


def softplus(x):
	x = np.asarray(x)
	return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)


def sigmoid(x):
	x = np.asarray(x)
	out = np.zeros_like(x)
	pos = x >= 0.0
	neg = np.logical_not(pos)
	out[pos] = 1.0/(1.0 + np.exp(-x[pos]))
	exp_x = np.exp(x[neg])
	out[neg] = exp_x/(1.0 + exp_x)
	return out


def softplus_inverse(y):
	y = np.asarray(y)
	return y + np.log(-np.expm1(-y))


def tau_sq_softplus_objective_and_gradient(theta, gene_regression_data, ridge_lambda):
	objective = 0.5*ridge_lambda*np.sum(np.square(theta[1:]))
	gradient = np.zeros(len(theta))
	gradient[1:] = ridge_lambda*theta[1:]

	for gene_data in gene_regression_data:
		latent_design_mat = gene_data['latent_design_mat']
		squared_ld = gene_data['squared_ld']
		sd_sq = gene_data['sd_sq']
		valid_rows = gene_data['valid_rows']
		response = gene_data['response']

		linear_pred = np.dot(latent_design_mat, theta)
		tau_sq_pred = softplus(linear_pred)
		tau_sq_grad = sigmoid(linear_pred)

		per_borzoi_weight = sd_sq*tau_sq_pred
		predicted_response = np.dot(squared_ld, per_borzoi_weight)
		residual = predicted_response[valid_rows] - response[valid_rows]

		objective = objective + 0.5*np.sum(np.square(residual))
		if np.sum(valid_rows) > 0:
			ld_weighted_residual = np.dot(squared_ld[valid_rows, :].T, residual)
			gradient = gradient + np.dot(latent_design_mat.T, (sd_sq*tau_sq_grad*ld_weighted_residual))

	return objective, gradient


def fit_tau_sq_spline_model(gene_id_to_ld_means, calibration_coef, n_basis=6, degree=3, ridge_lambda=0.0):
	all_abs_borzoi = []
	for gene_name in [*gene_id_to_ld_means]:
		all_abs_borzoi.append(np.abs(gene_id_to_ld_means[gene_name]['borzoi_effects_unstandardized']))
	all_abs_borzoi = np.hstack(all_abs_borzoi)

	fixed_interior_knots = np.asarray([0.001, 0.01, 0.1, 0.3, 0.6, 1.0])
	basis_info = create_bspline_basis_info(all_abs_borzoi, n_basis=n_basis, degree=degree, fixed_interior_knots=fixed_interior_knots)

	gene_regression_data = []
	all_response = []
	for gene_name in [*gene_id_to_ld_means]:
		gene_info = gene_id_to_ld_means[gene_name]
		residualized_eqtl = gene_info['eQTL_effect_sizes'] - (calibration_coef*gene_info['ld_means'])
		gene_response = np.square(residualized_eqtl) - np.square(gene_info['eQTL_effect_ses'])

		borzoi_abs = np.abs(gene_info['borzoi_effects_unstandardized'])
		spline_basis_mat = evaluate_bspline_basis(borzoi_abs, basis_info)
		latent_design_mat = create_latent_design_mat_from_spline_basis(spline_basis_mat)

		squared_ld = gene_info['squared_LD']
		sd_sq = np.square(gene_info['borzoi_variant_sdevs'])
		valid_rows = np.isfinite(gene_response)
		if np.sum(valid_rows) == 0:
			continue
		gene_regression_data.append({
			'latent_design_mat': latent_design_mat,
			'squared_ld': squared_ld,
			'sd_sq': sd_sq,
			'valid_rows': valid_rows,
			'response': gene_response
		})
		all_response.append(gene_response[valid_rows])

	if len(gene_regression_data) == 0:
		print('assumption eroror: no valid rows for tau spline regression')
		pdb.set_trace()

	all_response = np.hstack(all_response)
	all_spline_basis = evaluate_bspline_basis(all_abs_borzoi, basis_info)
	all_latent_design_mat = create_latent_design_mat_from_spline_basis(all_spline_basis)
	if basis_info['type'] == 'constant':
		mean_response = np.mean(all_response)
		init_tau_sq = np.maximum(mean_response, 1e-8)
		init_theta = np.asarray([softplus_inverse(np.asarray([init_tau_sq]))[0]])
	else:
		init_tau_sq = np.maximum(np.mean(all_response), 1e-8)
		init_theta = np.zeros(all_latent_design_mat.shape[1])
		init_theta[0] = softplus_inverse(np.asarray([init_tau_sq]))[0]

	def objective_wrapper(theta):
		return tau_sq_softplus_objective_and_gradient(theta, gene_regression_data, ridge_lambda)

	callback_info = {'iter': 0}

	def optimizer_callback(theta):
		callback_info['iter'] = callback_info['iter'] + 1
		objective, _ = objective_wrapper(theta)
		all_tau_sq_pred = softplus(np.dot(all_latent_design_mat, theta))
		tau_quantiles = np.quantile(all_tau_sq_pred, [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0])
		print(
			'Tau spline iter ' + str(callback_info['iter']) +
			': objective=' + str(objective) +
			', min=' + str(tau_quantiles[0]) +
			', p1=' + str(tau_quantiles[1]) +
			', p5=' + str(tau_quantiles[2]) +
			', median=' + str(tau_quantiles[3]) +
			', p95=' + str(tau_quantiles[4]) +
			', p99=' + str(tau_quantiles[5]) +
			', max=' + str(tau_quantiles[6]),
			flush=True
		)
		print(theta, flush=True)

	fit_res = minimize(
		fun=objective_wrapper,
		x0=init_theta,
		method='L-BFGS-B',
		jac=True,
		callback=optimizer_callback,
		options={'maxiter': 500}
	)
	tau_spline_theta = fit_res.x

	fitted_y = []
	for gene_data in gene_regression_data:
		linear_pred = np.dot(gene_data['latent_design_mat'], tau_spline_theta)
		tau_sq_pred = softplus(linear_pred)
		predicted_response = np.dot(gene_data['squared_ld'], gene_data['sd_sq']*tau_sq_pred)
		fitted_y.append(predicted_response[gene_data['valid_rows']])
	fitted_y = np.hstack(fitted_y)

	for gene_name in [*gene_id_to_ld_means]:
		gene_info = gene_id_to_ld_means[gene_name]
		borzoi_abs = np.abs(gene_info['borzoi_effects_unstandardized'])
		spline_basis_mat = evaluate_bspline_basis(borzoi_abs, basis_info)
		latent_design_mat = create_latent_design_mat_from_spline_basis(spline_basis_mat)
		tau_sq_pred = softplus(np.dot(latent_design_mat, tau_spline_theta))
		gene_info['tau_sq_pred'] = tau_sq_pred
		gene_info['tau_sq_pred_from_ld'] = np.dot(gene_info['squared_LD'], np.square(gene_info['borzoi_variant_sdevs'])*tau_sq_pred)
		gene_info['eqtl_residual_sq_minus_se_sq'] = np.square(gene_info['eQTL_effect_sizes'] - (calibration_coef*gene_info['ld_means'])) - np.square(gene_info['eQTL_effect_ses'])

	fit_info = {}
	fit_info['calibration_coef'] = calibration_coef
	fit_info['basis_info'] = basis_info
	fit_info['tau_spline_theta'] = tau_spline_theta
	fit_info['latent_design_n_cols'] = all_latent_design_mat.shape[1]
	fit_info['ridge_lambda'] = ridge_lambda
	fit_info['optimization_success'] = fit_res.success
	fit_info['optimization_status'] = fit_res.status
	fit_info['optimization_message'] = fit_res.message
	fit_info['optimization_objective'] = fit_res.fun
	fit_info['n_iterations'] = fit_res.nit
	fit_info['n_function_evals'] = fit_res.nfev
	fit_info['n_regression_obs'] = len(all_response)
	fit_info['rmse'] = np.sqrt(np.mean(np.square(all_response - fitted_y)))
	fit_info['mean_response'] = np.mean(all_response)
	fit_info['mean_fitted_response'] = np.mean(fitted_y)
	return gene_id_to_ld_means, fit_info


def save_tau_sq_spline_fit(output_stem, fit_info):
	summary_output_file = output_stem + '_tau_sq_spline_summary.txt'
	with open(summary_output_file, 'w') as t:
		t.write('field\tvalue\n')
		t.write('calibration_coef\t' + str(fit_info['calibration_coef']) + '\n')
		t.write('basis_type\t' + str(fit_info['basis_info']['type']) + '\n')
		t.write('n_basis\t' + str(fit_info['basis_info']['n_basis']) + '\n')
		t.write('degree\t' + str(fit_info['basis_info']['degree']) + '\n')
		t.write('x_min\t' + str(fit_info['basis_info']['x_min']) + '\n')
		t.write('x_max\t' + str(fit_info['basis_info']['x_max']) + '\n')
		t.write('knots\t' + ','.join(fit_info['basis_info']['knots'].astype(str)) + '\n')
		t.write('latent_design_n_cols\t' + str(fit_info['latent_design_n_cols']) + '\n')
		t.write('ridge_lambda\t' + str(fit_info['ridge_lambda']) + '\n')
		t.write('optimization_success\t' + str(fit_info['optimization_success']) + '\n')
		t.write('optimization_status\t' + str(fit_info['optimization_status']) + '\n')
		t.write('optimization_message\t' + str(fit_info['optimization_message']).replace('\n', ' ') + '\n')
		t.write('optimization_objective\t' + str(fit_info['optimization_objective']) + '\n')
		t.write('n_iterations\t' + str(fit_info['n_iterations']) + '\n')
		t.write('n_function_evals\t' + str(fit_info['n_function_evals']) + '\n')
		t.write('n_regression_obs\t' + str(fit_info['n_regression_obs']) + '\n')
		t.write('rmse\t' + str(fit_info['rmse']) + '\n')
		t.write('mean_response\t' + str(fit_info['mean_response']) + '\n')
		t.write('mean_fitted_response\t' + str(fit_info['mean_fitted_response']) + '\n')

	coef_output_file = output_stem + '_tau_sq_spline_coefs.txt'
	with open(coef_output_file, 'w') as t:
		t.write('basis_index\ttheta\n')
		for coef_iter, coef_val in enumerate(fit_info['tau_spline_theta']):
			t.write(str(coef_iter) + '\t' + str(coef_val) + '\n')

	grid_output_file = output_stem + '_tau_sq_spline_curve.txt'
	grid_x = np.linspace(fit_info['basis_info']['x_min'], fit_info['basis_info']['x_max'], 200)
	grid_spline_basis = evaluate_bspline_basis(grid_x, fit_info['basis_info'])
	grid_latent_design = create_latent_design_mat_from_spline_basis(grid_spline_basis)
	grid_tau = softplus(np.dot(grid_latent_design, fit_info['tau_spline_theta']))
	with open(grid_output_file, 'w') as t:
		t.write('abs_unstandardized_borzoi\tpred_tau_sq\n')
		for x_val, tau_val in zip(grid_x, grid_tau):
			t.write(str(x_val) + '\t' + str(tau_val) + '\n')

	fit_pickle_file = output_stem + '_tau_sq_spline_fit.pkl'
	with open(fit_pickle_file, 'wb') as f:
		pickle.dump(fit_info, f)


def unpack_tau_sq_nn_params(params, hidden_units):
	w1 = params[:hidden_units]
	b1 = params[hidden_units:(2*hidden_units)]
	w2 = params[(2*hidden_units):(3*hidden_units)]
	b2 = params[-1]
	return w1, b1, w2, b2


def tau_sq_nn_forward(z_input, params, hidden_units):
	z_input = np.asarray(z_input)
	w1, b1, w2, b2 = unpack_tau_sq_nn_params(params, hidden_units)
	hidden_linear = z_input[:, None]*w1[None, :] + b1[None, :]
	hidden_activation = np.tanh(hidden_linear)
	output_linear = np.dot(hidden_activation, w2) + b2
	tau_sq_pred = softplus(output_linear)
	return tau_sq_pred, hidden_activation, hidden_linear, output_linear


def tau_sq_nn_objective_and_gradient(params, gene_regression_data, hidden_units, ridge_lambda):
	w1, b1, w2, b2 = unpack_tau_sq_nn_params(params, hidden_units)
	objective = 0.5*ridge_lambda*(np.sum(np.square(w1)) + np.sum(np.square(w2)))
	grad_w1 = ridge_lambda*w1
	grad_b1 = np.zeros(hidden_units)
	grad_w2 = ridge_lambda*w2
	grad_b2 = 0.0

	for gene_data in gene_regression_data:
		z_input = gene_data['z_input']
		squared_ld = gene_data['squared_ld']
		sd_sq = gene_data['sd_sq']
		valid_rows = gene_data['valid_rows']
		response = gene_data['response']

		tau_sq_pred, hidden_activation, hidden_linear, output_linear = tau_sq_nn_forward(z_input, params, hidden_units)
		per_borzoi_weight = sd_sq*tau_sq_pred
		predicted_response = np.dot(squared_ld, per_borzoi_weight)
		residual = predicted_response[valid_rows] - response[valid_rows]

		objective = objective + 0.5*np.sum(np.square(residual))
		if np.sum(valid_rows) == 0:
			continue

		d_loss_d_predicted_response = np.zeros(len(predicted_response))
		d_loss_d_predicted_response[valid_rows] = residual
		d_loss_d_tau_sq = sd_sq*np.dot(squared_ld.T, d_loss_d_predicted_response)
		d_loss_d_output_linear = d_loss_d_tau_sq*sigmoid(output_linear)
		grad_b2 = grad_b2 + np.sum(d_loss_d_output_linear)
		grad_w2 = grad_w2 + np.dot(hidden_activation.T, d_loss_d_output_linear)

		d_loss_d_hidden = d_loss_d_output_linear[:, None]*w2[None, :]
		d_loss_d_hidden_linear = d_loss_d_hidden*(1.0 - np.square(hidden_activation))
		grad_b1 = grad_b1 + np.sum(d_loss_d_hidden_linear, axis=0)
		grad_w1 = grad_w1 + np.sum(z_input[:, None]*d_loss_d_hidden_linear, axis=0)

	gradient = np.hstack((grad_w1, grad_b1, grad_w2, np.asarray([grad_b2])))
	return objective, gradient


def fit_tau_sq_neural_net_model(gene_id_to_ld_means, calibration_coef, hidden_units=20, ridge_lambda=1e-6):
	all_abs_borzoi = []
	for gene_name in [*gene_id_to_ld_means]:
		all_abs_borzoi.append(np.abs(gene_id_to_ld_means[gene_name]['borzoi_effects_unstandardized']))
	all_abs_borzoi = np.hstack(all_abs_borzoi)
	x_mean = np.mean(all_abs_borzoi)
	x_sd = np.std(all_abs_borzoi)
	if np.isfinite(x_sd) is False or x_sd == 0.0:
		x_sd = 1.0
	all_z_input = (all_abs_borzoi - x_mean)/x_sd

	gene_regression_data = []
	all_response = []
	aa = []
	bb = []
	cc = []
	for gene_name in [*gene_id_to_ld_means]:
		gene_info = gene_id_to_ld_means[gene_name]
		residualized_eqtl = gene_info['eQTL_effect_sizes'] - (calibration_coef*gene_info['ld_means'])
		gene_response = np.square(residualized_eqtl) - np.square(gene_info['eQTL_effect_ses'])
		borzoi_abs = np.abs(gene_info['borzoi_effects_unstandardized'])
		z_input = (borzoi_abs - x_mean)/x_sd
		squared_ld = gene_info['squared_LD']
		sd_sq = np.square(gene_info['borzoi_variant_sdevs'])
		valid_rows = np.isfinite(gene_response)
		aa.append(gene_response)
		bb.append(squared_ld @ (borzoi_abs <= .75))
		cc.append(squared_ld @ (borzoi_abs > .75))
		if np.sum(valid_rows) == 0:
			continue
		gene_regression_data.append({
			'z_input': z_input,
			'squared_ld': squared_ld,
			'sd_sq': sd_sq,
			'valid_rows': valid_rows,
			'response': gene_response
		})
		all_response.append(gene_response[valid_rows])

	aa = np.hstack(aa)
	bb = np.hstack(bb)
	cc = np.hstack(cc)
	meany = np.dot(aa,bb)/np.dot(bb,bb)
	X = np.column_stack((bb, cc))
	beta_bb, beta_cc = np.linalg.lstsq(X, aa, rcond=None)[0]

	print(beta_bb, flush=True)
	print(beta_cc, flush=True)

	if len(gene_regression_data) == 0:
		print('assumption eroror: no valid rows for tau neural net regression')
		pdb.set_trace()

	all_response = np.hstack(all_response)
	init_tau_sq = np.maximum(np.mean(all_response), 1e-8)
	init_b2 = softplus_inverse(np.asarray([init_tau_sq]))[0]
	init_w1 = np.linspace(-1.0, 1.0, hidden_units)
	init_b1 = np.linspace(-0.5, 0.5, hidden_units)
	init_w2 = np.ones(hidden_units)*0.01
	init_params = np.hstack((init_w1, init_b1, init_w2, np.asarray([init_b2])))

	def objective_wrapper(params):
		return tau_sq_nn_objective_and_gradient(params, gene_regression_data, hidden_units, ridge_lambda)

	callback_info = {'iter': 0}

	def optimizer_callback(params):
		callback_info['iter'] = callback_info['iter'] + 1
		objective, _ = objective_wrapper(params)
		all_tau_sq_pred, _, _, _ = tau_sq_nn_forward(all_z_input, params, hidden_units)
		tau_quantiles = np.quantile(all_tau_sq_pred, [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0])
		print(
			'Tau NN iter ' + str(callback_info['iter']) +
			': objective=' + str(objective) +
			', min=' + str(tau_quantiles[0]) +
			', p1=' + str(tau_quantiles[1]) +
			', p5=' + str(tau_quantiles[2]) +
			', median=' + str(tau_quantiles[3]) +
			', p95=' + str(tau_quantiles[4]) +
			', p99=' + str(tau_quantiles[5]) +
			', max=' + str(tau_quantiles[6]),
			flush=True
		)
		print(params, flush=True)

	fit_res = minimize(
		fun=objective_wrapper,
		x0=init_params,
		method='L-BFGS-B',
		jac=True,
		callback=optimizer_callback,
		options={'maxiter': 2000}
	)
	nn_params = fit_res.x

	fitted_y = []
	for gene_data in gene_regression_data:
		tau_sq_pred, _, _, _ = tau_sq_nn_forward(gene_data['z_input'], nn_params, hidden_units)
		predicted_response = np.dot(gene_data['squared_ld'], gene_data['sd_sq']*tau_sq_pred)
		fitted_y.append(predicted_response[gene_data['valid_rows']])
	fitted_y = np.hstack(fitted_y)

	for gene_name in [*gene_id_to_ld_means]:
		gene_info = gene_id_to_ld_means[gene_name]
		borzoi_abs = np.abs(gene_info['borzoi_effects_unstandardized'])
		z_input = (borzoi_abs - x_mean)/x_sd
		tau_sq_pred, _, _, _ = tau_sq_nn_forward(z_input, nn_params, hidden_units)
		gene_info['tau_sq_pred'] = tau_sq_pred
		gene_info['tau_sq_pred_from_ld'] = np.dot(gene_info['squared_LD'], np.square(gene_info['borzoi_variant_sdevs'])*tau_sq_pred)
		gene_info['eqtl_residual_sq_minus_se_sq'] = np.square(gene_info['eQTL_effect_sizes'] - (calibration_coef*gene_info['ld_means'])) - np.square(gene_info['eQTL_effect_ses'])

	fit_info = {}
	fit_info['calibration_coef'] = calibration_coef
	fit_info['hidden_units'] = hidden_units
	fit_info['x_mean'] = x_mean
	fit_info['x_sd'] = x_sd
	fit_info['x_max'] = np.max(all_abs_borzoi)
	fit_info['nn_params'] = nn_params
	fit_info['ridge_lambda'] = ridge_lambda
	fit_info['optimization_success'] = fit_res.success
	fit_info['optimization_status'] = fit_res.status
	fit_info['optimization_message'] = fit_res.message
	fit_info['optimization_objective'] = fit_res.fun
	fit_info['n_iterations'] = fit_res.nit
	fit_info['n_function_evals'] = fit_res.nfev
	fit_info['n_regression_obs'] = len(all_response)
	fit_info['rmse'] = np.sqrt(np.mean(np.square(all_response - fitted_y)))
	fit_info['mean_response'] = np.mean(all_response)
	fit_info['mean_fitted_response'] = np.mean(fitted_y)
	return gene_id_to_ld_means, fit_info


def save_tau_sq_neural_net_fit(output_stem, fit_info):
	summary_output_file = output_stem + '_tau_sq_neural_net_summary.txt'
	with open(summary_output_file, 'w') as t:
		t.write('field\tvalue\n')
		t.write('calibration_coef\t' + str(fit_info['calibration_coef']) + '\n')
		t.write('hidden_units\t' + str(fit_info['hidden_units']) + '\n')
		t.write('x_mean\t' + str(fit_info['x_mean']) + '\n')
		t.write('x_sd\t' + str(fit_info['x_sd']) + '\n')
		t.write('x_max\t' + str(fit_info['x_max']) + '\n')
		t.write('ridge_lambda\t' + str(fit_info['ridge_lambda']) + '\n')
		t.write('optimization_success\t' + str(fit_info['optimization_success']) + '\n')
		t.write('optimization_status\t' + str(fit_info['optimization_status']) + '\n')
		t.write('optimization_message\t' + str(fit_info['optimization_message']).replace('\n', ' ') + '\n')
		t.write('optimization_objective\t' + str(fit_info['optimization_objective']) + '\n')
		t.write('n_iterations\t' + str(fit_info['n_iterations']) + '\n')
		t.write('n_function_evals\t' + str(fit_info['n_function_evals']) + '\n')
		t.write('n_regression_obs\t' + str(fit_info['n_regression_obs']) + '\n')
		t.write('rmse\t' + str(fit_info['rmse']) + '\n')
		t.write('mean_response\t' + str(fit_info['mean_response']) + '\n')
		t.write('mean_fitted_response\t' + str(fit_info['mean_fitted_response']) + '\n')

	param_output_file = output_stem + '_tau_sq_neural_net_params.txt'
	with open(param_output_file, 'w') as t:
		w1, b1, w2, b2 = unpack_tau_sq_nn_params(fit_info['nn_params'], fit_info['hidden_units'])
		t.write('param_type\tindex\tvalue\n')
		for ii, val in enumerate(w1):
			t.write('w1\t' + str(ii) + '\t' + str(val) + '\n')
		for ii, val in enumerate(b1):
			t.write('b1\t' + str(ii) + '\t' + str(val) + '\n')
		for ii, val in enumerate(w2):
			t.write('w2\t' + str(ii) + '\t' + str(val) + '\n')
		t.write('b2\t0\t' + str(b2) + '\n')

	curve_output_file = output_stem + '_tau_sq_neural_net_curve.txt'
	grid_x = np.linspace(0.0, fit_info['x_max'], 200)
	grid_z = (grid_x - fit_info['x_mean'])/fit_info['x_sd']
	grid_tau_sq, _, _, _ = tau_sq_nn_forward(grid_z, fit_info['nn_params'], fit_info['hidden_units'])
	with open(curve_output_file, 'w') as t:
		t.write('abs_unstandardized_borzoi\tpred_tau_sq\n')
		for x_val, tau_val in zip(grid_x, grid_tau_sq):
			t.write(str(x_val) + '\t' + str(tau_val) + '\n')

	fit_pickle_file = output_stem + '_tau_sq_neural_net_fit.pkl'
	with open(fit_pickle_file, 'wb') as f:
		pickle.dump(fit_info, f)


def fit_scaled_variance_ldscore_model(gene_id_to_ld_means, calibration_coef, calibration_coef_var=0.0, calibration_n_bootstraps=0, model_effect_scale='standardized'):
	if model_effect_scale not in ('standardized', 'allelic'):
		print('assumption error: model_effect_scale must be standardized or allelic')
		pdb.set_trace()

	gene_regression_data = []
	all_response = []
	all_uncorrected_response = []
	all_calibration_correction = []
	all_delta_sq = []
	design_mats = []
	for gene_name in [*gene_id_to_ld_means]:
		gene_info = gene_id_to_ld_means[gene_name]
		residualized_eqtl = gene_info['eQTL_effect_sizes'] - (calibration_coef*gene_info['ld_means'])
		gene_uncorrected_response = np.square(residualized_eqtl) - np.square(gene_info['eQTL_effect_ses'])
		gene_calibration_correction = calibration_coef_var*np.square(gene_info['ld_means'])
		gene_response = gene_uncorrected_response - gene_calibration_correction
		if model_effect_scale == 'standardized':
			delta_sq = np.square(gene_info['borzoi_effects'])
			ld_score_intercept = gene_info['ld_score_intercept']
			ld_score_delta_sq = gene_info['ld_score_delta_sq']
		else:
			delta_sq = np.square(gene_info['borzoi_effects_unstandardized'])
			ld_score_intercept = gene_info['ld_score_allelic_intercept']
			ld_score_delta_sq = gene_info['ld_score_allelic_delta_sq']
		valid_rows = np.isfinite(gene_response)
		gene_design_mat = np.transpose(np.vstack((
			np.ones(len(gene_response)),
			ld_score_intercept,
			ld_score_delta_sq
		)))

		if np.sum(valid_rows) == 0:
			continue
		gene_regression_data.append({
			'delta_sq': delta_sq,
			'valid_rows': valid_rows,
			'response': gene_response,
			'design_mat': gene_design_mat
		})
		all_response.append(gene_response[valid_rows])
		all_uncorrected_response.append(gene_uncorrected_response[valid_rows])
		all_calibration_correction.append(gene_calibration_correction[valid_rows])
		all_delta_sq.append(delta_sq)
		design_mats.append(gene_design_mat[valid_rows, :])

	if len(gene_regression_data) == 0:
		print('assumption eroror: no valid rows for ldscore scaled variance regression')
		pdb.set_trace()

	all_response = np.hstack(all_response)
	all_uncorrected_response = np.hstack(all_uncorrected_response)
	all_calibration_correction = np.hstack(all_calibration_correction)
	all_delta_sq = np.hstack(all_delta_sq)
	design_mat = np.vstack(design_mats)

	coef, residuals, rank, singular_vals = np.linalg.lstsq(design_mat, all_response, rcond=None)
	ldscore_intercept = coef[0]
	b_prior = coef[1]
	c_prior = coef[2]

	fitted_y = []
	for gene_data in gene_regression_data:
		predicted_response = np.dot(gene_data['design_mat'], coef)
		fitted_y.append(predicted_response[gene_data['valid_rows']])
	fitted_y = np.hstack(fitted_y)

	for gene_name in [*gene_id_to_ld_means]:
		gene_info = gene_id_to_ld_means[gene_name]
		if model_effect_scale == 'standardized':
			per_snp_var = b_prior + (c_prior*np.square(gene_info['borzoi_effects']))
			predicted_from_ld = ldscore_intercept + (b_prior*gene_info['ld_score_intercept']) + (c_prior*gene_info['ld_score_delta_sq'])
		else:
			per_snp_var = np.square(gene_info['borzoi_variant_sdevs'])*(b_prior + (c_prior*np.square(gene_info['borzoi_effects_unstandardized'])))
			predicted_from_ld = ldscore_intercept + (b_prior*gene_info['ld_score_allelic_intercept']) + (c_prior*gene_info['ld_score_allelic_delta_sq'])
		gene_info['scaled_variance_pred'] = per_snp_var
		gene_info['scaled_variance_pred_from_ld'] = predicted_from_ld
		gene_info['eqtl_residual_sq_minus_se_sq_uncorrected'] = np.square(gene_info['eQTL_effect_sizes'] - (calibration_coef*gene_info['ld_means'])) - np.square(gene_info['eQTL_effect_ses'])
		gene_info['calibration_uncertainty_correction'] = calibration_coef_var*np.square(gene_info['ld_means'])
		gene_info['eqtl_residual_sq_minus_se_sq'] = gene_info['eqtl_residual_sq_minus_se_sq_uncorrected'] - gene_info['calibration_uncertainty_correction']

	fit_info = {}
	fit_info['calibration_coef'] = calibration_coef
	fit_info['model_effect_scale'] = model_effect_scale
	fit_info['calibration_coef_var'] = calibration_coef_var
	fit_info['calibration_coef_se'] = np.sqrt(calibration_coef_var)
	fit_info['calibration_n_bootstraps'] = calibration_n_bootstraps
	fit_info['ldscore_intercept'] = ldscore_intercept
	fit_info['b_prior'] = b_prior
	fit_info['c_prior'] = c_prior
	fit_info['optimization_success'] = True
	fit_info['optimization_status'] = 0
	fit_info['optimization_message'] = 'ordinary_least_squares'
	fit_info['optimization_objective'] = 0.5*np.sum(np.square(all_response - fitted_y))
	fit_info['n_iterations'] = 0
	fit_info['n_function_evals'] = np.nan
	fit_info['linear_model_rank'] = rank
	fit_info['linear_model_singular_values'] = ','.join(singular_vals.astype(str))
	fit_info['n_regression_obs'] = len(all_response)
	fit_info['rmse'] = np.sqrt(np.mean(np.square(all_response - fitted_y)))
	fit_info['mean_response'] = np.mean(all_response)
	fit_info['mean_uncorrected_response'] = np.mean(all_uncorrected_response)
	fit_info['mean_calibration_uncertainty_correction'] = np.mean(all_calibration_correction)
	fit_info['mean_fitted_response'] = np.mean(fitted_y)
	fit_info['mean_predicted_prior_var'] = np.mean(b_prior + (c_prior*all_delta_sq))
	fit_info['median_predicted_prior_var'] = np.median(b_prior + (c_prior*all_delta_sq))
	fit_info['min_predicted_prior_var'] = np.min(b_prior + (c_prior*all_delta_sq))
	fit_info['max_predicted_prior_var'] = np.max(b_prior + (c_prior*all_delta_sq))
	return gene_id_to_ld_means, fit_info


def save_scaled_variance_ldscore_fit(output_stem, fit_info):
	summary_output_file = output_stem + '_ldscore_scaled_variance_summary.txt'
	with open(summary_output_file, 'w') as t:
		t.write('field\tvalue\n')
		t.write('calibration_coef\t' + str(fit_info['calibration_coef']) + '\n')
		t.write('model_effect_scale\t' + str(fit_info['model_effect_scale']) + '\n')
		t.write('calibration_coef_var\t' + str(fit_info['calibration_coef_var']) + '\n')
		t.write('calibration_coef_se\t' + str(fit_info['calibration_coef_se']) + '\n')
		t.write('calibration_n_bootstraps\t' + str(fit_info['calibration_n_bootstraps']) + '\n')
		t.write('ldscore_intercept\t' + str(fit_info['ldscore_intercept']) + '\n')
		t.write('b_prior\t' + str(fit_info['b_prior']) + '\n')
		t.write('c_prior\t' + str(fit_info['c_prior']) + '\n')
		t.write('optimization_success\t' + str(fit_info['optimization_success']) + '\n')
		t.write('optimization_status\t' + str(fit_info['optimization_status']) + '\n')
		t.write('optimization_message\t' + str(fit_info['optimization_message']).replace('\n', ' ') + '\n')
		t.write('optimization_objective\t' + str(fit_info['optimization_objective']) + '\n')
		t.write('n_iterations\t' + str(fit_info['n_iterations']) + '\n')
		t.write('n_function_evals\t' + str(fit_info['n_function_evals']) + '\n')
		t.write('linear_model_rank\t' + str(fit_info['linear_model_rank']) + '\n')
		t.write('linear_model_singular_values\t' + str(fit_info['linear_model_singular_values']) + '\n')
		t.write('n_regression_obs\t' + str(fit_info['n_regression_obs']) + '\n')
		t.write('rmse\t' + str(fit_info['rmse']) + '\n')
		t.write('mean_response\t' + str(fit_info['mean_response']) + '\n')
		t.write('mean_uncorrected_response\t' + str(fit_info['mean_uncorrected_response']) + '\n')
		t.write('mean_calibration_uncertainty_correction\t' + str(fit_info['mean_calibration_uncertainty_correction']) + '\n')
		t.write('mean_fitted_response\t' + str(fit_info['mean_fitted_response']) + '\n')
		t.write('mean_predicted_prior_var\t' + str(fit_info['mean_predicted_prior_var']) + '\n')
		t.write('median_predicted_prior_var\t' + str(fit_info['median_predicted_prior_var']) + '\n')
		t.write('min_predicted_prior_var\t' + str(fit_info['min_predicted_prior_var']) + '\n')
		t.write('max_predicted_prior_var\t' + str(fit_info['max_predicted_prior_var']) + '\n')

	param_output_file = output_stem + '_ldscore_scaled_variance_params.txt'
	with open(param_output_file, 'w') as t:
		t.write('parameter\tvalue\n')
		t.write('a_prior\t' + str(fit_info['calibration_coef']) + '\n')
		t.write('model_effect_scale\t' + str(fit_info['model_effect_scale']) + '\n')
		t.write('a_prior_var\t' + str(fit_info['calibration_coef_var']) + '\n')
		t.write('a_prior_se\t' + str(fit_info['calibration_coef_se']) + '\n')
		t.write('ldscore_intercept\t' + str(fit_info['ldscore_intercept']) + '\n')
		t.write('b_prior\t' + str(fit_info['b_prior']) + '\n')
		t.write('c_prior\t' + str(fit_info['c_prior']) + '\n')

	fit_pickle_file = output_stem + '_ldscore_scaled_variance_fit.pkl'
	with open(fit_pickle_file, 'wb') as f:
		pickle.dump(fit_info, f)


##########################
# Command line args
##########################
est_borzoi_effect_size_file = sys.argv[1]
onek_genomes_plink_filestem = sys.argv[2]
genotype_sample_mapping_file = sys.argv[3]
est_eqtl_effect_size_file = sys.argv[4]
anno_method = sys.argv[5]
borzoi_based_prior_output_stem = sys.argv[6]
model_effect_scale = 'standardized'
if len(sys.argv) > 7 and sys.argv[7] in ('standardized', 'allelic'):
	model_effect_scale = sys.argv[7]
	if len(sys.argv) > 8:
		calibration_n_bootstraps = int(sys.argv[8])
	else:
		calibration_n_bootstraps = 200
elif len(sys.argv) > 7:
	calibration_n_bootstraps = int(sys.argv[7])
else:
	calibration_n_bootstraps = 200



##############################
# Load in data
##############################
# Create mapping from gene id to vector of est (standardized) eqtl effect sizes
gene_id_to_est_eqtl_effects, variant_id_to_genotype_sdev = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file)

# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev,standardize=True)
gene_id_to_est_borzoi_effects_unstandardized = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev,standardize=False)

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

# Generate per gene ld-means
gene_id_to_ld_means = generate_gene_ld_means(gene_id_to_est_borzoi_effects, gene_id_to_est_borzoi_effects_unstandardized, gene_id_to_est_eqtl_effects, onek_genomes_plink_filestem, genotype_sample_indices, variant_id_to_genotype_sdev, anno_method)
del gene_id_to_est_eqtl_effects
del gene_id_to_est_borzoi_effects
del gene_id_to_est_borzoi_effects_unstandardized
del variant_id_to_genotype_sdev


'''
pickle_output_file = borzoi_based_prior_output_stem + '_gene_id_to_ld_means.pkl'

with open(pickle_output_file, 'wb') as f:
	pickle.dump(gene_id_to_ld_means, f)

pickle_output_file = borzoi_based_prior_output_stem + '_gene_id_to_ld_means.pkl'
with open(pickle_output_file, 'rb') as f:
	gene_id_to_ld_means = pickle.load(f)
'''


calibration_coef, calibration_coef_var, calibration_coef_se = compute_calibration_coefficient_bootstrap_var(gene_id_to_ld_means, n_bootstraps=calibration_n_bootstraps)
print(calibration_coef, flush=True)
print('calibration_coef_se: ' + str(calibration_coef_se), flush=True)

gene_id_to_ld_means, scaled_variance_fit_info = fit_scaled_variance_ldscore_model(gene_id_to_ld_means, calibration_coef, calibration_coef_var=calibration_coef_var, calibration_n_bootstraps=calibration_n_bootstraps, model_effect_scale=model_effect_scale)
print('model_effect_scale: ' + str(scaled_variance_fit_info['model_effect_scale']), flush=True)
print('ldscore_intercept: ' + str(scaled_variance_fit_info['ldscore_intercept']), flush=True)
print('b_prior: ' + str(scaled_variance_fit_info['b_prior']), flush=True)
print('c_prior: ' + str(scaled_variance_fit_info['c_prior']), flush=True)
'''
with open(pickle_output_file, 'wb') as f:
	pickle.dump(gene_id_to_ld_means, f)
'''
save_scaled_variance_ldscore_fit(borzoi_based_prior_output_stem, scaled_variance_fit_info)
