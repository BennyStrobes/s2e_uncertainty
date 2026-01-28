import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink1_bin
import pyarrow.parquet as pq
import pickle


def extract_dictionary_list_of_protein_coding_genes(protein_coding_genes_file, chrom_string):
	f = open(protein_coding_genes_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != chrom_string:
			continue
		dicti[data[3]] = 1

	f.close()
	return dicti


def extract_borzoi_effect_sizes(borzoi_effect_size_file_stem, chrom_num, pc_genes, borzoi_effect_size_file_parts):
	chrom_string = 'chr' + str(chrom_num)
	gene_obj = {}

	for borzoi_effect_size_file_part in borzoi_effect_size_file_parts:
		borzoi_effect_size_file = borzoi_effect_size_file_stem + borzoi_effect_size_file_part + '.txt'
		f = open(borzoi_effect_size_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count+1
				continue
			if chrom_string != data[0]:
				continue

			variant_name = data[2]
			gene_name = data[3]
			if gene_name not in pc_genes:
				continue

			mean_borzoi_effect = float(data[4])
			bs_borzoi_effects = np.asarray(data[5].split(';')).astype(float)

			bs_minor_count = np.min([np.sum(bs_borzoi_effects > 0.0), np.sum(bs_borzoi_effects < 0.0)])

			if gene_name not in gene_obj:
				gene_obj[gene_name] = []
			gene_obj[gene_name].append((variant_name, mean_borzoi_effect, np.std(bs_borzoi_effects), bs_minor_count, bs_borzoi_effects))


		f.close()


	return gene_obj



def extract_eqtl_ss_for_this_chromosome(eqtl_ss_file, borzoi_effect_size_gene_obj):
	eqtl_ss_obj = {}
	pf = pq.ParquetFile(eqtl_ss_file)
	for rg in range(pf.num_row_groups):
		table = pf.read_row_group(rg)   # this is a chunk

		# Process fields and print to output
		array = np.asarray(table)

		# Loop through snp-gene pairs 
		nrows = array.shape[0]
		for row_iter in range(nrows):
			data = array[row_iter, :]
			variant_id = data[1]
			gene_id = data[0]
			effect_size = float(data[-2])
			std_err = float(data[-1])
			if gene_id not in borzoi_effect_size_gene_obj:
				continue
			if gene_id not in eqtl_ss_obj:
				eqtl_ss_obj[gene_id] = []
			eqtl_ss_obj[gene_id].append((variant_id, effect_size, std_err))
	pf.close()
	return eqtl_ss_obj


def extract_eqtl_ss_for_gene(list_of_tuples):
	variants = []
	betas = []
	beta_ses = []

	for tupler in list_of_tuples:
		variants.append(tupler[0])
		betas.append(tupler[1])
		beta_ses.append(tupler[2])

	return np.asarray(variants), np.asarray(betas), np.asarray(beta_ses)


def extract_borzoi_effects_for_gene(list_of_tuples):
	variants = []
	betas = []
	beta_ses = []
	beta_mab = []
	beta_bs = []
	for tupler in list_of_tuples:
		variants.append(tupler[0])
		betas.append(tupler[1])
		beta_ses.append(tupler[2])
		beta_mab.append(tupler[3])
		beta_bs.append(tupler[4])

	return np.asarray(variants), np.asarray(betas), np.asarray(beta_ses), np.asarray(beta_mab), np.asarray(beta_bs)

def extract_dictionary_list_of_hm3_rsids(hm3_snp_list_file):
	f = open(hm3_snp_list_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		dicti[line] = 1
	f.close()
	return dicti

def extract_genotype_matrix_for_subset_of_snps(G_obj_geno, G_obj_gtex_ids_dicti, gene_eqtl_variants):
	valid_variants = []
	geno = []
	mafs = []
	for variant_id in gene_eqtl_variants:
		var_info = variant_id.split('_')
		alt_variant_id = var_info[0] + '_' + var_info[1] + '_' + var_info[3] + '_' + var_info[2] + '_b38'

		if variant_id in G_obj_gtex_ids_dicti:
			tmp_geno = G_obj_geno[:, G_obj_gtex_ids_dicti[variant_id]]
			geno.append(tmp_geno)
			valid_variants.append(True)
			mafs.append(np.sum(tmp_geno)/(2.0*len(tmp_geno)))
		elif alt_variant_id in G_obj_gtex_ids_dicti:
			tmp_geno = G_obj_geno[:, G_obj_gtex_ids_dicti[alt_variant_id]]
			geno.append(-tmp_geno)
			valid_variants.append(True)
			mafs.append(np.sum(tmp_geno)/(2.0*len(tmp_geno)))
		else:
			valid_variants.append(False)

	geno = np.asarray(geno)
	valid_variants = np.asarray(valid_variants)
	mafs = np.asarray(mafs)

	return geno, valid_variants, mafs


def rowwise_cross_correlation(A, B, eps=1e-12):
	"""
	Compute correlation matrix between rows of A and rows of B.

	Parameters
	----------
	A : array-like, shape (K1, N)
		First genotype matrix
	B : array-like, shape (K2, N)
		Second genotype matrix
	eps : float
		Small value to avoid division by zero

	Returns
	-------
	C : ndarray, shape (K1, K2)
		Correlation matrix
	"""
	# Convert to float for safety
	A = np.asarray(A, dtype=float)
	B = np.asarray(B, dtype=float)

	# Center rows
	A -= A.mean(axis=1, keepdims=True)
	B -= B.mean(axis=1, keepdims=True)

	# Compute row-wise std devs
	A_std = np.sqrt(np.sum(A**2, axis=1, keepdims=True))
	B_std = np.sqrt(np.sum(B**2, axis=1, keepdims=True))

	# Normalize
	A /= (A_std + eps)
	B /= (B_std + eps)

	# Correlation = dot product
	return A @ B.T

def mean_ci_sem(x, alpha=0.05):
	x = np.asarray(x, dtype=float)
	n = x.size
	
	mean = x.mean()
	sd = x.std(ddof=1)          # sample standard deviation
	sem = sd / np.sqrt(n)
	
	z = 1.96  # 95% CI normal approximation
	ci_lower = mean - z * sem
	ci_upper = mean + z * sem
	
	print(str(mean) + '\t' + str(ci_lower) + '\t' + str(ci_upper))
	return

def reorder_vector(new_variant_vec, old_effect_size_vec, old_variant_vec):
	dicti = {}
	for ii, variant_id in enumerate(old_variant_vec):
		dicti[variant_id] = old_effect_size_vec[ii]

	new_effect_size_vec = []
	for variant_id in new_variant_vec:
		new_effect_size_vec.append(dicti[variant_id])


	return np.asarray(new_effect_size_vec)

def reorder_matrix(new_variant_vec, old_effect_size_vec, old_variant_vec):
	dicti = {}
	for ii, variant_id in enumerate(old_variant_vec):
		dicti[variant_id] = old_effect_size_vec[ii,:]

	new_effect_size_vec = []
	for variant_id in new_variant_vec:
		new_effect_size_vec.append(dicti[variant_id])


	return np.asarray(new_effect_size_vec)

def rank_based_pvalue(bs_corry):
	np_p_plus = (1 + np.sum(bs_corry >= 0))/(1+len(bs_corry))
	np_p_minus = (1 + np.sum(bs_corry <=0))/(1+len(bs_corry))
	borzoi_pvalues = 2.0*np.min((np_p_plus, np_p_minus))
	return borzoi_pvalues

def rank_based_pvalue_array(bs_marginal_effects):
	np_p_plus = (1 + np.sum(bs_marginal_effects >= 0,axis=1))/(1+bs_marginal_effects.shape[1])
	np_p_minus = (1 + np.sum(bs_marginal_effects <=0,axis=1))/(1+bs_marginal_effects.shape[1])
	borzoi_pvalues = 2.0*np.min((np_p_plus, np_p_minus),axis=0)
	return borzoi_pvalues


def compute_correlation(filtered_std_gene_borzoi_effects, filtered_std_gene_eqtl_beta, R_mat):
	corry_num = np.dot(filtered_std_gene_borzoi_effects, filtered_std_gene_eqtl_beta) 
	corry_denom = np.sqrt(np.dot(np.dot(filtered_std_gene_borzoi_effects, R_mat), filtered_std_gene_borzoi_effects))
	corry = corry_num/corry_denom
	return corry


###########################
# Command line args
###########################
tissue_name = sys.argv[1]
borzoi_effect_size_file_stem = sys.argv[2]
eqtl_sumstats_dir = sys.argv[3]
protein_coding_genes_file = sys.argv[4]
genotype_1000G_plink_stem = sys.argv[5]
hm3_snp_list_file = sys.argv[6]
expr_correlation_summary_file = sys.argv[7]

borzoi_effect_size_file_parts = np.asarray(['0', '1'])

t = open(expr_correlation_summary_file,'w')
t.write('gene_id\texpression_correlation\tmax_abs_eqtl_z\tsignal_variance\tnoise_variance\tmax_abs_borzoi\tcorrelation_bs_pvalue\texpression_correlation_z_score_pred\texpression_correlation_ridge_predictor\tmin_borzoi_pvalue\texpression_correlation_uncertainty_denom\texpression_correlation_best_sig_marginal\texpression_correlation_sign_flip_weighted\texpression_correlation_low_ld_pred\n')


# Loop through chromosomes
vals = []
vals2 = []
for chrom_num in range(1,9):

	print(chrom_num)

	# Extract dictionary list of hapmap3 rsids
	hm3rsids = extract_dictionary_list_of_hm3_rsids(hm3_snp_list_file)


	# Extract protein coding genes
	pc_genes = extract_dictionary_list_of_protein_coding_genes(protein_coding_genes_file, 'chr' + str(chrom_num))

	# Extract borzoi effect sizes
	#borzoi_effect_size_gene_obj = {}
	borzoi_effect_size_gene_obj = extract_borzoi_effect_sizes(borzoi_effect_size_file_stem, chrom_num, pc_genes, borzoi_effect_size_file_parts)




	'''
	tmp_file = '/lab-share/CHIP-Strober-e2/Public/ben/tmp/tmper.pkl'
	f = open(tmp_file, "wb")
	pickle.dump(borzoi_effect_size_gene_obj, f)
	f.close()
	tmp_file = '/lab-share/CHIP-Strober-e2/Public/ben/tmp/tmper.pkl'
	f = open(tmp_file, "rb")
	borzoi_effect_size_gene_obj = pickle.load(f)
	f.close()
	'''


	# Extract eqtl sumstats file for this chromosome
	eqtl_ss_file = eqtl_sumstats_dir + tissue_name + '.v10.allpairs.chr' + str(chrom_num) + '.parquet'
	eqtl_ss_obj = extract_eqtl_ss_for_this_chromosome(eqtl_ss_file, borzoi_effect_size_gene_obj)

	pdb.set_trace()
	'''
	tmp_file = '/lab-share/CHIP-Strober-e2/Public/ben/tmp/tmper2.pkl'
	f = open(tmp_file, "wb")
	pickle.dump(eqtl_ss_obj, f)
	f.close()
	tmp_file = '/lab-share/CHIP-Strober-e2/Public/ben/tmp/tmper2.pkl'
	f = open(tmp_file, "rb")
	eqtl_ss_obj = pickle.load(f)
	f.close()
	'''

	# Load in reference panel genotype data
	G_obj = read_plink1_bin(genotype_1000G_plink_stem + str(chrom_num) + ".bed", genotype_1000G_plink_stem + str(chrom_num) +".bim", genotype_1000G_plink_stem + str(chrom_num) +".fam", verbose=False)
	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	G_obj_rsid = np.asarray(G_obj.snp)
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	G_obj_gtex_ids = []
	G_obj_gtex_ids_dicti = {}
	gtex_variant_to_hm3_bool = {}
	for ii, posser in enumerate(G_obj_pos):
		tmp_gtex_id = 'chr' + G_obj_chrom[ii] + '_' + str(posser) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii] + '_b38'
		G_obj_gtex_ids.append(tmp_gtex_id)
		G_obj_gtex_ids_dicti[tmp_gtex_id] = ii

		hm3_bool = False
		if G_obj_rsid[ii] in hm3rsids:
			hm3_bool = True
		gtex_variant_to_hm3_bool[tmp_gtex_id] = hm3_bool
		tmp_gtex_id_alt = 'chr' + G_obj_chrom[ii] + '_' + str(posser) + '_' + G_obj_a1[ii] + '_' + G_obj_a0[ii] + '_b38'
		gtex_variant_to_hm3_bool[tmp_gtex_id_alt] = hm3_bool

	G_obj_gtex_ids = np.asarray(G_obj_gtex_ids)



	# Loop through genes
	for geneid in [*borzoi_effect_size_gene_obj]:

		# Skip genes that we don't have summary stats for
		if geneid not in eqtl_ss_obj:
			continue

		# Extract eqtl sum stats for gene
		gene_eqtl_variants, gene_eqtl_beta, gene_eqtl_beta_se = extract_eqtl_ss_for_gene(eqtl_ss_obj[geneid])

		# Extract borzoi effects for gene
		gene_borzoi_variants, gene_borzoi_effects, gene_borzoi_effect_se, gene_borzoi_effect_mab, gene_borzoi_effect_bs = extract_borzoi_effects_for_gene(borzoi_effect_size_gene_obj[geneid])


		# Extract genotype matrix for eqtl variants
		eqtl_variant_geno, valid_eqtl_variants, eqtl_variant_afs = extract_genotype_matrix_for_subset_of_snps(G_obj_geno, G_obj_gtex_ids_dicti, gene_eqtl_variants)
		gene_eqtl_variants = gene_eqtl_variants[valid_eqtl_variants]
		gene_eqtl_beta = gene_eqtl_beta[valid_eqtl_variants]
		gene_eqtl_beta_se = gene_eqtl_beta_se[valid_eqtl_variants]
		hm3_valid = []
		for var in gene_eqtl_variants:
			#hm3_valid.append(gtex_variant_to_hm3_bool[var])
			hm3_valid.append(True)
		hm3_valid = np.asarray(hm3_valid)
		gene_eqtl_variants = gene_eqtl_variants[hm3_valid]
		gene_eqtl_beta = gene_eqtl_beta[hm3_valid]
		gene_eqtl_beta_se = gene_eqtl_beta_se[hm3_valid]
		eqtl_variant_geno = eqtl_variant_geno[hm3_valid,:]
		eqtl_variant_afs = eqtl_variant_afs[hm3_valid]
		gene_eqtl_zeds = gene_eqtl_beta/gene_eqtl_beta_se
		std_gene_eqtl_beta = gene_eqtl_beta*np.sqrt(2.0*eqtl_variant_afs*(1.0-eqtl_variant_afs))


		gene_eqtl_zed = gene_eqtl_beta/gene_eqtl_beta_se

		# Extract genotype matrix for borzoi variants
		borzoi_variant_geno, valid_borzoi_variants, borzoi_variant_afs = extract_genotype_matrix_for_subset_of_snps(G_obj_geno, G_obj_gtex_ids_dicti, gene_borzoi_variants)
		gene_borzoi_variants = gene_borzoi_variants[valid_borzoi_variants]
		gene_borzoi_effects = gene_borzoi_effects[valid_borzoi_variants]
		gene_borzoi_effect_se = gene_borzoi_effect_se[valid_borzoi_variants]
		gene_borzoi_effect_mab = gene_borzoi_effect_mab[valid_borzoi_variants]
		gene_borzoi_effect_bs = gene_borzoi_effect_bs[valid_borzoi_variants, :]
		std_gene_borzoi_effects = gene_borzoi_effects*np.sqrt(2.0*borzoi_variant_afs*(1.0-borzoi_variant_afs))
		std_gene_borzoi_effects_bs = np.transpose(np.transpose(gene_borzoi_effect_bs)*np.sqrt(2.0*borzoi_variant_afs*(1.0-borzoi_variant_afs)))

		# Get overlaps
		overlapping_variants = np.asarray(list(set(gene_eqtl_variants) & set(gene_borzoi_variants)))

		if len(overlapping_variants) < 5:
			continue


		filtered_std_gene_borzoi_effects = reorder_vector(overlapping_variants, std_gene_borzoi_effects, gene_borzoi_variants)
		filtered_gene_borzoi_effect_mab = reorder_vector(overlapping_variants, gene_borzoi_effect_mab, gene_borzoi_variants)
		filtered_gene_borzoi_effects = reorder_vector(overlapping_variants, gene_borzoi_effects, gene_borzoi_variants)
		filtered_std_gene_borzoi_effects_bs = reorder_matrix(overlapping_variants, std_gene_borzoi_effects_bs, gene_borzoi_variants)
		filtered_std_gene_eqtl_beta = reorder_vector(overlapping_variants, std_gene_eqtl_beta, gene_eqtl_variants)
		filtered_eqtl_zeds = reorder_vector(overlapping_variants, gene_eqtl_zeds, gene_eqtl_variants)

		

		overlapping_variant_geno, valid_overlapping_variants, overlapping_variant_afs = extract_genotype_matrix_for_subset_of_snps(G_obj_geno, G_obj_gtex_ids_dicti, overlapping_variants)

		R_mat = np.corrcoef(overlapping_variant_geno)


		'''
		#zero_indices = (filtered_gene_borzoi_effect_mab != 0.0) 
		zero_indices = (np.abs(filtered_gene_borzoi_effects) < .25)
		#zero_indices = (filtered_gene_borzoi_effect_mab != 0.0) | (np.abs(filtered_gene_borzoi_effects) < .05)

		thresh_std_gene_borzoi_effects = np.copy(filtered_std_gene_borzoi_effects)
		thresh_std_gene_borzoi_effects[zero_indices] = 0.0
		#if np.sum(zero_indices == False) == 0:
			#continue

		if np.std(thresh_std_gene_borzoi_effects) != 0.0:
			corry_num = np.dot(thresh_std_gene_borzoi_effects, filtered_std_gene_eqtl_beta) 
			corry_denom = np.sqrt(np.dot(np.dot(thresh_std_gene_borzoi_effects, R_mat), thresh_std_gene_borzoi_effects))
			corry_num2 = np.dot(filtered_std_gene_borzoi_effects, filtered_std_gene_eqtl_beta) 
			corry_denom2 = np.sqrt(np.dot(np.dot(filtered_std_gene_borzoi_effects, R_mat), filtered_std_gene_borzoi_effects))

			corry = corry_num/corry_denom
			corry2 = corry_num2/corry_denom2
			vals.append(corry)
			vals2.append(corry2)
			print('#####')
			print(corry)
			mean_ci_sem(vals)
			print(corry2)
			mean_ci_sem(vals2)
		'''


		bs_marginal_effects = np.dot(R_mat, filtered_std_gene_borzoi_effects_bs)
		same_sign = (np.all(bs_marginal_effects >= 0, axis=1) | np.all(bs_marginal_effects <= 0, axis=1))
		bs_marginal_effects_pvalue = rank_based_pvalue_array(bs_marginal_effects)

		np_p_plus = (1 + np.sum(filtered_std_gene_borzoi_effects_bs >= 0,axis=1))/(1+filtered_std_gene_borzoi_effects_bs.shape[1])
		np_p_minus = (1 + np.sum(filtered_std_gene_borzoi_effects_bs <=0,axis=1))/(1+filtered_std_gene_borzoi_effects_bs.shape[1])
		borzoi_pvalues = 2.0*np.min((np_p_plus, np_p_minus),axis=0)
		

		signal = np.dot(np.dot(filtered_std_gene_borzoi_effects, R_mat), filtered_std_gene_borzoi_effects)
		sampling_cov = np.cov(filtered_std_gene_borzoi_effects_bs)
		noise = np.trace(np.dot(R_mat, sampling_cov))
		signal_to_noise_ratio = signal/noise
		#print(signal_to_noise_ratio)
		
		corry = compute_correlation(filtered_std_gene_borzoi_effects, filtered_std_gene_eqtl_beta, R_mat)


		max_abs_borzoi = np.max(np.abs(filtered_gene_borzoi_effects))
		max_abs_eqtl_zed = np.max(np.abs(filtered_eqtl_zeds))

		# bs distribution on correlation
		bs_corry_num = filtered_std_gene_borzoi_effects_bs.T @ filtered_std_gene_eqtl_beta
		Rb = R_mat @ filtered_std_gene_borzoi_effects_bs
		bs_corry_denom = np.sqrt(np.sum(filtered_std_gene_borzoi_effects_bs * Rb, axis=0))
		bs_corry = bs_corry_num / bs_corry_denom
		bs_corry_pvalue = rank_based_pvalue(bs_corry)


		filtered_std_gene_borzoi_effects_z_scored = np.mean(filtered_std_gene_borzoi_effects_bs,axis=1)/np.std(filtered_std_gene_borzoi_effects_bs,axis=1)
		corry_z_score_pred = compute_correlation(filtered_std_gene_borzoi_effects_z_scored, filtered_std_gene_eqtl_beta, R_mat)



		tau_sq = 1e-10
		filtered_std_gene_borzoi_effects_ridged = np.mean(filtered_std_gene_borzoi_effects_bs,axis=1)*(tau_sq/(tau_sq + np.var(filtered_std_gene_borzoi_effects_bs,axis=1)))
		corry_ridge_pred = compute_correlation(filtered_std_gene_borzoi_effects_ridged, filtered_std_gene_eqtl_beta, R_mat)



		RW = R_mat @ filtered_std_gene_borzoi_effects_bs                  # (K, B)
		sampled_variances = np.sum(filtered_std_gene_borzoi_effects_bs * RW, axis=0)
		total_var = np.mean(sampled_variances)

		corry_num = np.dot(filtered_std_gene_borzoi_effects, filtered_std_gene_eqtl_beta) 
		corry_uncertainty_den = corry_num/np.sqrt(total_var)


		# Best sig marginal
		best_sig_marginal = np.zeros(bs_marginal_effects.shape[0])
		if np.min(bs_marginal_effects_pvalue) < .05:
			indices = bs_marginal_effects_pvalue < .05
			best_sig_marginal[indices] = filtered_std_gene_borzoi_effects[indices]
		corry_best_sig_marginal_pred = compute_correlation(best_sig_marginal, filtered_std_gene_eqtl_beta, R_mat)


		variances = np.var(filtered_std_gene_borzoi_effects_bs,axis=1)
		squared_means = np.square(filtered_std_gene_borzoi_effects)
		weights = 1.0/(1.0 + (variances/squared_means))
		corry_sign_flip_weighted = compute_correlation(filtered_std_gene_borzoi_effects*weights, filtered_std_gene_eqtl_beta, R_mat)



		indices = (borzoi_pvalues < .025) & (np.partition(np.abs(R_mat), -2, axis=0)[-2, :] < .8)
		low_ld_pred = np.zeros(bs_marginal_effects.shape[0])
		low_ld_pred[indices] = filtered_std_gene_borzoi_effects[indices]
		corry_low_ld_pred = compute_correlation(low_ld_pred, filtered_std_gene_eqtl_beta, R_mat)



		'''
		if geneid == 'ENSG00000163466.16' or geneid == 'ENSG00000069020.19' or geneid == 'ENSG00000152942.19' or geneid == 'ENSG00000165055.16':
			pdb.set_trace()
		'''

		t.write(geneid + '\t' + str(corry) + '\t' + str(max_abs_eqtl_zed) + '\t' + str(signal) + '\t' + str(noise) + '\t' + str(max_abs_borzoi) + '\t' + str(bs_corry_pvalue) + '\t' + str(corry_z_score_pred) + '\t' + str(corry_ridge_pred) + '\t' + str(np.min(borzoi_pvalues)) + '\t' + str(corry_uncertainty_den) + '\t' + str(corry_best_sig_marginal_pred) + '\t' + str(corry_sign_flip_weighted) + '\t' + str(corry_low_ld_pred) + '\n')
		t.flush()

		'''

		# Extract genotype matrix for borzoi variants
		borzoi_variant_geno, valid_borzoi_variants, borzoi_variant_afs = extract_genotype_matrix_for_subset_of_snps(G_obj_geno, G_obj_gtex_ids_dicti, gene_borzoi_variants)
		gene_borzoi_variants = gene_borzoi_variants[valid_borzoi_variants]
		gene_borzoi_effects = gene_borzoi_effects[valid_borzoi_variants]
		gene_borzoi_effect_se = gene_borzoi_effect_se[valid_borzoi_variants]
		gene_borzoi_effect_mab = gene_borzoi_effect_mab[valid_borzoi_variants]
		std_gene_borzoi_effects = gene_borzoi_effects*np.sqrt(2.0*borzoi_variant_afs*(1.0-borzoi_variant_afs))

		zero_indices = (gene_borzoi_effect_mab != 0.0)
		zero_indices = np.abs(gene_borzoi_effects) < .25

		thresh_std_gene_borzoi_effects = np.copy(std_gene_borzoi_effects)
		
		thresh_std_gene_borzoi_effects[zero_indices] = 0.0
		#if np.sum(zero_indices == False) == 0:
			#continue

		if eqtl_variant_geno.shape[0] < 1 or borzoi_variant_geno.shape[0] < 1:
			continue


		R_mat = rowwise_cross_correlation(eqtl_variant_geno, borzoi_variant_geno)


		if np.std(thresh_std_gene_borzoi_effects) != 0.0:
			#borzoi_full_pred = np.dot(R_mat, thresh_std_gene_borzoi_effects)
			#corry = np.corrcoef(gene_eqtl_zeds,borzoi_full_pred)[0,1]

			pdb.set_trace()
			corry_num = np.dot(thresh_std_gene_borzoi_effects, std_gene_eqtl_beta) 
			corry_denom = np.dot(np.dot(thresh_std_gene_borzoi_effects, R_mat), thresh_std_gene_borzoi_effects)

			#print(geneid + '\t' + 'borzoi_sig' + '\t' + str(corry))
			vals.append(corry)
			print('#####')
			print(corry)
			mean_ci_sem(vals)
		'''
t.close()





