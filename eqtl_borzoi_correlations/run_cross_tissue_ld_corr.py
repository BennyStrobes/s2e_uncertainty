import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time





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
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, pos, a0, a1, effect, std_se, genotype_sdev)
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

def generate_gene_ld_means(gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, gene_id_to_variant_gene_anno, onek_genomes_plink_filestem, anno_names, genotype_sample_indices):
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
			if gene_id not in gene_id_to_variant_gene_anno:
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

			# Anno
			variant_anno, borzoi_anno_variant_alleles = load_in_snp_gene_anno_data(ordered_cis_variants, gene_id_to_variant_gene_anno[gene_id], len(anno_names))

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
				if np.isnan(borzoi_effects[var_index]) == False:
					if borzoi_variant_alleles[var_index,:][0] == geno_alleles[1]:
						borzoi_effects[var_index] = -1.0*borzoi_effects[var_index]
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
			variant_anno = variant_anno[observed_borzoi_indices, :]
			LD = LD[:, observed_borzoi_indices]

			# Add to dictionary
			if gene_id in gene_to_ld_means:
				print('assumption erororo')
				pdb.set_trace()

			NNN = geno_mat.shape[1]
			sq_LD = LD**2
			adj_sq_ld = sq_LD - (1.0 - sq_LD)/(NNN - 2.0)

			gene_to_ld_means[gene_id] = {}
			gene_to_ld_means[gene_id]['eQTL_effect_sizes'] = eqtl_effects
			gene_to_ld_means[gene_id]['eQTL_effect_ses'] = eqtl_effect_ses
			gene_to_ld_means[gene_id]['borzoi_effects'] = borzoi_effects
			gene_to_ld_means[gene_id]['variant_anno'] = variant_anno
			gene_to_ld_means[gene_id]['ld_means'] = LD @ (variant_anno * borzoi_effects[:, None])
			gene_to_ld_means[gene_id]['ld_scores'] = (adj_sq_ld) @ variant_anno


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


def compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names, intercept=False):
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

	if intercept:
		xx_intercept = np.hstack((np.ones((xx.shape[0], 1)), xx))
		coefs, _, _, _ = np.linalg.lstsq(xx_intercept, yy, rcond=None)
		intercept = coefs[0]
		taus = coefs[1:]
	else:
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
	avg_per_snp_h2 = compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names, intercept=True)

	return avg_calibration_slopes, avg_per_snp_h2, avg_borzoi_variances


def create_mapping_from_gene_id_to_est_delta_eqtl_effect_sizes(gene_id_to_t1_est_eqtl_effects, gene_id_to_t2_est_eqtl_effects, variant_id_to_t1_genotype_sdev, variant_id_to_t2_genotype_sdev):
	mapping = {}
	variant_id_to_geno_sdev = {}
	# Loop through genes
	for gene_id in [*gene_id_to_t1_est_eqtl_effects]:
		if gene_id not in gene_id_to_t2_est_eqtl_effects:
			continue

		# Loop through variants
		for var_id in gene_id_to_t1_est_eqtl_effects[gene_id]:
			if var_id not in gene_id_to_t2_est_eqtl_effects[gene_id]:
				continue

			t1_info = gene_id_to_t1_est_eqtl_effects[gene_id][var_id]
			t2_info = gene_id_to_t2_est_eqtl_effects[gene_id][var_id]

			if t1_info[3] != t2_info[3] or t1_info[4] != t2_info[4] or t1_info[5] != t2_info[5]:
				print('assumption erroror')
				pdb.set_trace()


			t1_effect_size = t1_info[6]
			t1_effect_size_se = t1_info[7]
			t1_genotype_sdev = t1_info[8]
			t2_effect_size = t2_info[6]
			t2_effect_size_se = t2_info[7]
			t2_genotype_sdev = t2_info[8]

			# Should be more legit, but probs ok for now (more imporant that it is shared)
			joint_genotype_sdev = np.mean([t1_genotype_sdev, t2_genotype_sdev])

			t1_std_effect_size = t1_effect_size*joint_genotype_sdev
			t1_std_effect_size_se = t1_effect_size_se*joint_genotype_sdev
			t2_std_effect_size = t2_effect_size*joint_genotype_sdev
			t2_std_effect_size_se = t2_effect_size_se*joint_genotype_sdev

			delta_effect_size = t1_std_effect_size - t2_std_effect_size
			delta_effect_size_se = np.sqrt(np.square(t1_std_effect_size_se) + np.square(t2_std_effect_size_se))

			if gene_id not in mapping:
				mapping[gene_id] = {}
			if var_id in mapping[gene_id]:
				print('repeat snp error')
				pdb.set_trace()

			mapping[gene_id][var_id] = (gene_id, var_id, t1_info[2], t1_info[3], t1_info[4], t1_info[5], delta_effect_size, delta_effect_size_se)

			variant_id_to_geno_sdev[var_id] = joint_genotype_sdev

	return mapping, variant_id_to_geno_sdev


def create_mapping_from_gene_id_to_delta_causal_effects(gene_id_to_t1_est_borzoi_effects, gene_id_to_t2_est_borzoi_effects):
	mapping = {}

	for gene_id in [*gene_id_to_t1_est_borzoi_effects]:
		if gene_id not in gene_id_to_t2_est_borzoi_effects:
			continue
		for var_id in gene_id_to_t1_est_borzoi_effects[gene_id]:
			if var_id not in gene_id_to_t2_est_borzoi_effects[gene_id]:
				continue

			if gene_id not in mapping:
				mapping[gene_id] = {}
			if var_id in mapping[gene_id]:
				print('repeat snp error')
				pdb.set_trace()

			t1_info = gene_id_to_t1_est_borzoi_effects[gene_id][var_id]
			t2_info = gene_id_to_t2_est_borzoi_effects[gene_id][var_id]

			if t1_info[0] != t2_info[0] or t1_info[1] != t2_info[1] or t1_info[2] != t2_info[2] or t1_info[3] != t2_info[3] or t1_info[4] != t2_info[4] or t1_info[5] != t2_info[5]:
				print('asssumptionoenrron')
				pdb.set_trace()

			delta_score = t1_info[6] - t2_info[6]

			mapping[gene_id][var_id] = (t1_info[0], t1_info[1], t1_info[2], t1_info[3], t1_info[4], t1_info[5], delta_score)

	return mapping


def create_mapping_from_gene_id_to_delta_variant_gene_annotations(gene_id_to_t1_est_borzoi_effects, gene_id_to_t2_est_borzoi_effects):
	mapping = {}

	for gene_id in [*gene_id_to_t1_est_borzoi_effects]:
		if gene_id not in gene_id_to_t2_est_borzoi_effects:
			continue
		for var_id in gene_id_to_t1_est_borzoi_effects[gene_id]:
			if var_id not in gene_id_to_t2_est_borzoi_effects[gene_id]:
				continue

			if gene_id not in mapping:
				mapping[gene_id] = {}
			if var_id in mapping[gene_id]:
				print('repeat snp error')
				pdb.set_trace()

			t1_info = gene_id_to_t1_est_borzoi_effects[gene_id][var_id]
			t2_info = gene_id_to_t2_est_borzoi_effects[gene_id][var_id]

			if t1_info[0] != t2_info[0] or t1_info[1] != t2_info[1] or t1_info[2] != t2_info[2] or t1_info[3] != t2_info[3] or t1_info[4] != t2_info[4] or t1_info[5] != t2_info[5]:
				print('asssumptionoenrron')
				pdb.set_trace()

			top_bin = np.max([np.argmax(t1_info[6]), np.argmax(t2_info[6])])
			new_bins = np.zeros(len(t1_info[6]))
			new_bins[top_bin] = 1.0

			mapping[gene_id][var_id] = (t1_info[0], t1_info[1], t1_info[2], t1_info[3], t1_info[4], t1_info[5], new_bins)

	return mapping





##########################
# Command line args
##########################
t1_est_borzoi_effect_size_file = sys.argv[1]
t2_est_borzoi_effect_size_file = sys.argv[2]
t1_est_eqtl_effect_size_file = sys.argv[3]
t2_est_eqtl_effect_size_file = sys.argv[4]
t1_sim_variant_gene_annotation_file = sys.argv[5]
t2_sim_variant_gene_annotation_file = sys.argv[6]
onek_genomes_plink_filestem = sys.argv[7]
t1_genotype_sample_mapping_file = sys.argv[8]
t2_genotype_sample_mapping_file = sys.argv[9]
ld_corr_output_stem = sys.argv[10]




##############################
# Load in data
##############################
# Create mapping from gene id to vector of est (standardized) eqtl effect sizes in each tissue
gene_id_to_t1_est_eqtl_effects, variant_id_to_t1_genotype_sdev = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(t1_est_eqtl_effect_size_file)
gene_id_to_t2_est_eqtl_effects, variant_id_to_t2_genotype_sdev = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(t2_est_eqtl_effect_size_file)

# Create mapping from gene id to vector of est (standardized) delta eqtl effect sizes between two tissues
print('start')
gene_id_to_est_delta_eqtl_effects, variant_id_to_genotype_sdev = create_mapping_from_gene_id_to_est_delta_eqtl_effect_sizes(gene_id_to_t1_est_eqtl_effects, gene_id_to_t2_est_eqtl_effects, variant_id_to_t1_genotype_sdev, variant_id_to_t2_genotype_sdev)
del gene_id_to_t1_est_eqtl_effects
del gene_id_to_t2_est_eqtl_effects
del variant_id_to_t1_genotype_sdev
del variant_id_to_t2_genotype_sdev

'''
# TMP LOADING
delta_eqtl_effects_pickle_file = ld_corr_output_stem + '_gene_id_to_est_delta_eqtl_effects.pkl'
with open(delta_eqtl_effects_pickle_file, 'wb') as f:
	pickle.dump(gene_id_to_est_delta_eqtl_effects, f)

genotype_sdev_pickle_file = ld_corr_output_stem + '_variant_id_to_genotype_sdev.pkl'
with open(genotype_sdev_pickle_file, 'wb') as f:
	pickle.dump(variant_id_to_genotype_sdev, f)

delta_eqtl_effects_pickle_file = ld_corr_output_stem + '_gene_id_to_est_delta_eqtl_effects.pkl'
genotype_sdev_pickle_file = ld_corr_output_stem + '_variant_id_to_genotype_sdev.pkl'

with open(delta_eqtl_effects_pickle_file, 'rb') as f:
	gene_id_to_est_delta_eqtl_effects = pickle.load(f)

with open(genotype_sdev_pickle_file, 'rb') as f:
	variant_id_to_genotype_sdev = pickle.load(f)
'''



# Create mapping from gene id to vector of est borzoi effects
gene_id_to_t1_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(t1_est_borzoi_effect_size_file, variant_id_to_genotype_sdev,standardize=True)
gene_id_to_t2_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(t2_est_borzoi_effect_size_file, variant_id_to_genotype_sdev,standardize=True)
gene_id_to_delta_est_borzoi_effects = create_mapping_from_gene_id_to_delta_causal_effects(gene_id_to_t1_est_borzoi_effects, gene_id_to_t2_est_borzoi_effects)
del gene_id_to_t1_est_borzoi_effects
del gene_id_to_t2_est_borzoi_effects


# TMP Loading
'''
delta_borzoi_effects_pickle_file = ld_corr_output_stem + '_gene_id_to_est_delta_borzoi_effects.pkl'
with open(delta_borzoi_effects_pickle_file, 'wb') as f:
	pickle.dump(gene_id_to_delta_est_borzoi_effects, f)
delta_borzoi_effects_pickle_file = ld_corr_output_stem + '_gene_id_to_est_delta_borzoi_effects.pkl'
with open(delta_borzoi_effects_pickle_file, 'rb') as f:
	gene_id_to_delta_est_borzoi_effects = pickle.load(f)
'''


# Create mapping from gene id to vector of variant-gene annotations
gene_id_to_t1_variant_gene_anno, t1_anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(t1_sim_variant_gene_annotation_file)
gene_id_to_t2_variant_gene_anno, t2_anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(t2_sim_variant_gene_annotation_file)

if np.array_equal(t1_anno_names, t2_anno_names) == False:
	print('assumptione ororor')
	pdb.set_trace()
anno_names = np.copy(t1_anno_names)
gene_id_to_variant_gene_delta_anno = create_mapping_from_gene_id_to_delta_variant_gene_annotations(gene_id_to_t1_variant_gene_anno, gene_id_to_t2_variant_gene_anno)
del gene_id_to_t1_variant_gene_anno
del gene_id_to_t2_variant_gene_anno
del t1_anno_names
del t2_anno_names

'''
# TMP LOADING
delta_anno_pickle_file = ld_corr_output_stem + '_gene_id_to_delta_anno.pkl'
with open(delta_anno_pickle_file, 'wb') as f:
	pickle.dump(gene_id_to_variant_gene_delta_anno, f)
delta_anno_pickle_file = ld_corr_output_stem + '_gene_id_to_delta_anno.pkl'
with open(delta_anno_pickle_file, 'rb') as f:
	gene_id_to_variant_gene_delta_anno = pickle.load(f)
anno_names = np.asarray(['magnitude_bin0', 'magnitude_bin1', 'magnitude_bin2', 'magnitude_bin3', 'magnitude_bin4'])
'''


# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices1 = (np.loadtxt(t1_genotype_sample_mapping_file)).astype(int)
genotype_sample_indices2 = (np.loadtxt(t2_genotype_sample_mapping_file)).astype(int)
genotype_sample_indices = np.sort(np.unique(np.hstack((genotype_sample_indices1, genotype_sample_indices2))))
del genotype_sample_indices1
del genotype_sample_indices2


# Generate per gene ld-means
gene_id_to_ld_means = generate_gene_ld_means(gene_id_to_delta_est_borzoi_effects, gene_id_to_est_delta_eqtl_effects, gene_id_to_variant_gene_delta_anno, onek_genomes_plink_filestem, anno_names, genotype_sample_indices)
del gene_id_to_est_delta_eqtl_effects
del gene_id_to_variant_gene_delta_anno
del gene_id_to_delta_est_borzoi_effects


'''
# TMP Loading
pickle_file = ld_corr_output_stem + '.pkl'
with open(pickle_file, 'wb') as f:
	pickle.dump(gene_id_to_ld_means, f)
pickle_file = ld_corr_output_stem + '.pkl'
with open(pickle_file, 'rb') as f:
	gene_id_to_ld_means = pickle.load(f)
anno_names = np.asarray(['magnitude_bin0', 'magnitude_bin1', 'magnitude_bin2', 'magnitude_bin3', 'magnitude_bin4'])
'''


##############################
# Run regressions
##############################

# Run standard regression
ordered_gene_names = np.asarray([*gene_id_to_ld_means])
avg_calibration_slopes, avg_per_snp_eqtl_h2, avg_borzoi_variances = run_ld_corr(ordered_gene_names, gene_id_to_ld_means)
corrs = avg_calibration_slopes*np.sqrt(avg_borzoi_variances/avg_per_snp_eqtl_h2)



# Bootstrap estimates!
n_bs = 100
bs_calibration_slopes = []
bs_per_snp_eqtl_h2 = []
bs_borzoi_variances = []
for bs_iter in range(n_bs):
	print(bs_iter)
	bs_genes = np.random.choice(ordered_gene_names, size=len(ordered_gene_names), replace=True)
	tmp_bs_calibration_slopes, tmp_bs_per_snp_eqtl_h2, tmp_bs_borzoi_variances = run_ld_corr(bs_genes, gene_id_to_ld_means)

	bs_calibration_slopes.append(tmp_bs_calibration_slopes)
	bs_per_snp_eqtl_h2.append(tmp_bs_per_snp_eqtl_h2)
	bs_borzoi_variances.append(tmp_bs_borzoi_variances)

bs_calibration_slopes = np.asarray(bs_calibration_slopes)
bs_per_snp_eqtl_h2 = np.asarray(bs_per_snp_eqtl_h2)
bs_borzoi_variances = np.asarray(bs_borzoi_variances)

# Compute correlation information
observed_corr = np.full(len(anno_names), np.nan)
for anno_iter in range(len(anno_names)):
	sqrt_term = avg_borzoi_variances[anno_iter]/avg_per_snp_eqtl_h2[anno_iter]
	if np.isfinite(sqrt_term) and sqrt_term >= 0.0:
		observed_corr[anno_iter] = avg_calibration_slopes[anno_iter]*np.sqrt(sqrt_term)
bs_corr = np.full(bs_calibration_slopes.shape, np.nan)
for bs_iter in range(bs_calibration_slopes.shape[0]):
	for anno_iter in range(bs_calibration_slopes.shape[1]):
		sqrt_term = bs_borzoi_variances[bs_iter, anno_iter]/bs_per_snp_eqtl_h2[bs_iter, anno_iter]
		if np.isfinite(sqrt_term) and sqrt_term >= 0.0:
			bs_corr[bs_iter, anno_iter] = bs_calibration_slopes[bs_iter, anno_iter]*np.sqrt(sqrt_term)




#############################
# Print results to output file
##############################

results_output_file = ld_corr_output_stem + '_bootstrap_summary.txt'
with open(results_output_file, 'w') as t:
	t.write('sample_label\tsample_iter\tannotation_name\tcalibration_slope\tper_snp_eqtl_h2\tborzoi_variance\tcorrelation\n')
	for anno_iter, anno_name in enumerate(anno_names):
		t.write('observed\t-1\t' + str(anno_name) + '\t' + str(avg_calibration_slopes[anno_iter]) + '\t' + str(avg_per_snp_eqtl_h2[anno_iter]) + '\t' + str(avg_borzoi_variances[anno_iter]) + '\t' + str(observed_corr[anno_iter]) + '\n')

	for bs_iter in range(n_bs):
		for anno_iter, anno_name in enumerate(anno_names):
			t.write('bootstrap\t' + str(bs_iter) + '\t' + str(anno_name) + '\t' + str(bs_calibration_slopes[bs_iter, anno_iter]) + '\t' + str(bs_per_snp_eqtl_h2[bs_iter, anno_iter]) + '\t' + str(bs_borzoi_variances[bs_iter, anno_iter]) + '\t' + str(bs_corr[bs_iter, anno_iter]) + '\n')

summary_stats_output_file = ld_corr_output_stem + '_bootstrap_stats.txt'
with open(summary_stats_output_file, 'w') as t:
	t.write('annotation_name\toutput_name\tmean\tbootstrapped_mean\tbootstrap_se\tgaussian_z_score\tempirical_ci_lower\tempirical_ci_upper\n')
	output_names = ['calibration_slope', 'per_snp_eqtl_h2', 'borzoi_variance', 'correlation']
	observed_mat = np.vstack((avg_calibration_slopes, avg_per_snp_eqtl_h2, avg_borzoi_variances, observed_corr)).T
	bootstrap_mats = [bs_calibration_slopes, bs_per_snp_eqtl_h2, bs_borzoi_variances, bs_corr]
	for anno_iter, anno_name in enumerate(anno_names):
		for output_iter, output_name in enumerate(output_names):
			bootstrap_distribution = bootstrap_mats[output_iter][:, anno_iter]
			bootstrap_mean = np.mean(bootstrap_distribution)
			bootstrap_se = np.std(bootstrap_distribution, ddof=1)
			if bootstrap_se == 0.0:
				gaussian_z_score = np.nan
			else:
				gaussian_z_score = observed_mat[anno_iter, output_iter]/bootstrap_se
			empirical_ci = np.quantile(bootstrap_distribution, [.025, .975])
			t.write(str(anno_name) + '\t' + output_name + '\t' + str(observed_mat[anno_iter, output_iter]) + '\t' + str(bootstrap_mean) + '\t' + str(bootstrap_se) + '\t' + str(gaussian_z_score) + '\t' + str(empirical_ci[0]) + '\t' + str(empirical_ci[1]) + '\n')

print(results_output_file)
print(summary_stats_output_file)
