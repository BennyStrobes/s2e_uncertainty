import numpy as np
import os
import sys
import pdb
import gzip
from pandas_plink import read_plink
from scipy.special import logsumexp
import pickle
import time
from numba import njit

def load_in_cis_variants(gene_cis_variant_file):
	f = open(gene_cis_variant_file)
	ordered_variant_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ordered_variant_names.append(data[1])
	f.close()

	return np.asarray(ordered_variant_names)


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
		effect = float(data[2])
		if gene_id not in mapping:
			mapping[gene_id] = []
		mapping[gene_id].append((var_id, effect))
	f.close()
	return mapping

def create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file, eqtl_sample_size):
	desired_se = np.sqrt(1.0/eqtl_sample_size)
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
		effect = float(data[2])
		se = float(data[3])
		if gene_id not in mapping:
			mapping[gene_id] = []
		zed = effect/se
		std_effect = zed*desired_se
		mapping[gene_id].append((var_id, std_effect))
	f.close()
	return mapping

def load_in_effects(tupler):
	var_arr = []
	effect_arr = []

	for tuply in tupler:
		var_arr.append(tuply[0])
		effect_arr.append(tuply[1])



	return np.asarray(effect_arr), np.asarray(var_arr)


def extract_ldsc_vectors(gene_ld_summary_file, onek_genomes_plink_filestem, gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, use_adjusted_LD=False):
	marginal_eqtl_vec = []
	marginal_borzoi_vec = []
	ld_score_vec = []
	gene_to_marginal_borzoi_vec = {}

	# Loop through chromosomes
	for chrom_num in range(21,22):
		print(chrom_num)
		# string of chromosome name
		chrom_string = 'chr' + str(chrom_num)

		# Load in chromosome plink data
		(bim, fam, G) = read_plink(onek_genomes_plink_filestem + str(chrom_num))

		# Create mapping from variant id to index
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))

		# Loop through genes / filter to only genes on this chromosome
		f = open(gene_ld_summary_file)
		head_count = 0
		counter = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			line_chrom_num = data[1]
			if line_chrom_num != chrom_string:
				continue
			gene_name = data[0]
			counter = counter + 1
			gene_cis_variant_file = data[4]

			cis_variant_names = load_in_cis_variants(gene_cis_variant_file)

			# Extract variant indices corresponding to cis_variant names
			cis_variant_indices = []
			for variant_name in cis_variant_names:
				cis_variant_indices.append(rsid_to_genotype_index[variant_name])
			cis_variant_indices = np.asarray(cis_variant_indices)

			# Extract genotype matrix
			geno_mat = G[cis_variant_indices,:].compute()
			n_samp = geno_mat.shape[1]

			# Create LD matrix
			LD = np.corrcoef(geno_mat)

			#LD2 = np.load(gene_id_to_ld_file[gene_name])
			#print(np.max(np.abs(LD-LD2)))

			# Get squared LD
			squared_LD = np.square(LD)
			squared_LD_adj = squared_LD - ((1.0-squared_LD)/(n_samp-2.0))

			# Compute Ld-scores
			if use_adjusted_LD:
				ld_scores = np.sum(squared_LD_adj,axis=1)
			else:
				ld_scores = np.sum(squared_LD,axis=1)
			
			# Load in est eqtl effect sizes
			est_eqtl_effect_sizes, tmp_variant_names1 = load_in_effects(gene_id_to_est_eqtl_effects[gene_name])

			# Load in est borzoi causal effect sizes
			est_borzoi_causal_effect_sizes, tmp_variant_names2 = load_in_effects(gene_id_to_est_borzoi_effects[gene_name])
			
			# Quick error checking
			if np.array_equal(tmp_variant_names1,cis_variant_names) == False or np.array_equal(tmp_variant_names2,cis_variant_names) == False:
				print('variant mismatch error')
				pdb.set_trace()

			# Get borzoi predicted marginal effect sizes
			est_borzoi_marginal_effects = np.dot(LD, est_borzoi_causal_effect_sizes)

			gene_to_marginal_borzoi_vec[gene_name] = est_borzoi_marginal_effects

			for var_iter, variant_name in enumerate(cis_variant_names):
				marginal_eqtl_vec.append(est_eqtl_effect_sizes[var_iter])
				marginal_borzoi_vec.append(est_borzoi_marginal_effects[var_iter])
				ld_score_vec.append(ld_scores[var_iter])

		f.close()
	return gene_to_marginal_borzoi_vec, np.asarray(marginal_eqtl_vec), np.asarray(marginal_borzoi_vec), np.asarray(ld_score_vec)


def create_mapping_from_gene_id_to_ld_file(gene_ld_summary_file):
	f = open(gene_ld_summary_file)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ens_id = data[0]
		ld_file = data[6]
		if ens_id in mapping:
			print('assumptioneornrornroo')
			pdb.set_trace()
		mapping[ens_id] = ld_file

	f.close()

	return mapping

def update_causal_effects_for_single_gene(cur_beta, resid_delta, omega_grid, log_pis, eqtl_sample_size, RR, LL):
	n_snps = len(cur_beta)
	gene_zz = np.zeros((n_snps,LL))

	# Component-specific posterior params (shared across all snps)
	precisions = eqtl_sample_size + (1.0 / omega_grid)
	post_vars = 1.0 / precisions


	# Looop through snps in gene (in random order)
	for snp_iter in np.random.permutation(np.arange(n_snps)):

		# c_j = r_j + R_jj * m_j (Un-residualize current causal effect)
		cj = resid_delta[snp_iter] + 1.0*cur_beta[snp_iter]

		# Component-specific posterior params
		post_means = (eqtl_sample_size * cj) / precisions

		# log responsibilities up to normalization
		log_weights = (
			log_pis
			- 0.5 * np.log(omega_grid)
			+ 0.5 * np.log(post_vars)
			+ 0.5 * (post_means ** 2) / post_vars
		)

		log_weights = log_weights - logsumexp(log_weights)
		resp = np.exp(log_weights)

		# posterior mean of eta_j
		old_beta = cur_beta[snp_iter]
		new_beta = np.sum(resp * post_means)
		cur_beta[snp_iter] = new_beta

		# MAP component assignment
		gene_zz[snp_iter, :] = resp

		# residual update: r <- r - R_.j * delta
		delta_beta = new_beta - old_beta
		if delta_beta != 0.0:
			resid_delta = resid_delta - RR[:, snp_iter] * delta_beta
	return cur_beta, gene_zz, resid_delta

@njit
def _update_causal_effects_for_single_gene_numba_core(
	cur_beta,
	resid_delta,
	RR,
	snp_order,
	eqtl_sample_size,
	precisions,
	post_vars,
	log_weight_const,
	resid_var
):
	n_snps = cur_beta.shape[0]
	LL = precisions.shape[0]

	gene_zz = np.zeros((n_snps, LL))

	post_means = np.empty(LL)
	log_weights = np.empty(LL)
	resp = np.empty(LL)

	for ii in range(n_snps):
		snp_iter = snp_order[ii]

		# assuming RR diagonal is 1
		cj = resid_delta[snp_iter] + cur_beta[snp_iter]

		max_logw = -1e300
		for ll in range(LL):
			pm = ((eqtl_sample_size/resid_var) * cj) / precisions[ll]
			post_means[ll] = pm
			lw = log_weight_const[ll] + 0.5 * (pm * pm) / post_vars[ll]
			log_weights[ll] = lw
			if lw > max_logw:
				max_logw = lw

		denom = 0.0
		for ll in range(LL):
			val = np.exp(log_weights[ll] - max_logw)
			resp[ll] = val
			denom += val

		new_beta = 0.0
		for ll in range(LL):
			resp[ll] = resp[ll] / denom
			gene_zz[snp_iter, ll] = resp[ll]
			new_beta += resp[ll] * post_means[ll]

		old_beta = cur_beta[snp_iter]
		cur_beta[snp_iter] = new_beta

		delta_beta = new_beta - old_beta
		if delta_beta != 0.0:
			for kk in range(n_snps):
				resid_delta[kk] -= RR[kk, snp_iter] * delta_beta

	return cur_beta, gene_zz, resid_delta


def update_causal_effects_for_single_gene_numba(
	cur_beta,
	resid_delta,
	omega_grid,
	log_pis,
	eqtl_sample_size,
	RR,
	LL,
	resid_var,
):
	cur_beta = np.asarray(cur_beta, dtype=np.float64).copy()
	resid_delta = np.asarray(resid_delta, dtype=np.float64).copy()
	omega_grid = np.asarray(omega_grid, dtype=np.float64)
	log_pis = np.asarray(log_pis, dtype=np.float64)
	RR = np.asarray(RR, dtype=np.float64)

	precisions = (eqtl_sample_size/resid_var) + (1.0 / omega_grid)
	post_vars = 1.0 / precisions
	log_weight_const = log_pis - 0.5 * np.log(omega_grid) + 0.5 * np.log(post_vars)

	snp_order = np.random.permutation(cur_beta.shape[0]).astype(np.int64)

	return _update_causal_effects_for_single_gene_numba_core(
		cur_beta,
		resid_delta,
		RR,
		snp_order,
		float(eqtl_sample_size),
		precisions,
		post_vars,
		log_weight_const,
		resid_var,
	)



def run_ashr_inference(gene_id_to_est_marginal_eqtl_effects, gene_id_to_est_marginal_borzoi_effects, gene_id_to_ld_file, eqtl_sample_size, ldsc_tau2, LL=20, max_iter=500, resid_var=1.0):
	###############################
	# Model initialization
	###############################
	gene_names = np.sort(np.asarray([*gene_id_to_est_marginal_borzoi_effects]))
	# Initialize causal effects
	gene_id_to_resid_delta = {}
	gene_id_to_beta = {}
	gene_id_to_zz = {}
	for gene_name in gene_names:
		resid_delta = np.copy(gene_id_to_est_marginal_borzoi_effects[gene_name]) - np.copy(gene_id_to_est_marginal_eqtl_effects[gene_name])
		gene_id_to_resid_delta[gene_name] = resid_delta
		gene_id_to_beta[gene_name] = np.zeros(len(resid_delta))
		gene_id_to_zz[gene_name] = np.zeros((len(resid_delta),LL))
	# Set up grid
	sigma_grid = np.concatenate([
		np.array([0.0]),
		np.exp(np.linspace(np.log(1e-3), np.log(1.0), LL-1))
	])
	omega_grid = sigma_grid**2 + 1e-10
	pis = np.zeros(LL) + .0005
	pis[np.argmin(np.abs(omega_grid - ldsc_tau2)):] = 1
	pis = pis/np.sum(pis)
	log_pis = np.log(pis)


	###############################
	# Iterative expectation-maximization algorithm
	###############################
	for itera in range(max_iter):

		print(itera)
		print(pis)
		print(np.sum(pis*omega_grid))
		t1 = time.time()
		###########################
		# E-step
		###########################
		# Loop through genes
		for gene_name in gene_names:
			#gene_beta, gene_zz, gene_resid_delta = update_causal_effects_for_single_gene(gene_id_to_beta[gene_name], gene_id_to_resid_delta[gene_name], omega_grid, log_pis, eqtl_sample_size, np.load(gene_id_to_ld_file[gene_name]), LL)
			gene_beta, gene_zz, gene_resid_delta = update_causal_effects_for_single_gene_numba(gene_id_to_beta[gene_name], gene_id_to_resid_delta[gene_name], omega_grid, log_pis, eqtl_sample_size, np.load(gene_id_to_ld_file[gene_name]), LL,resid_var)
			gene_id_to_beta[gene_name] = gene_beta
			gene_id_to_resid_delta[gene_name] = gene_resid_delta
			gene_id_to_zz[gene_name] = gene_zz

		###########################
		# M-step
		###########################
		if itera > 10:
			counts = np.zeros(LL) + 1e-10
			for gene_name in gene_names:
				counts = counts + np.sum(gene_id_to_zz[gene_name],axis=0)
			pis = counts/np.sum(counts)
			log_pis = np.log(pis)
		t2 = time.time()
		print(t2-t1)
	
	# Open output object
	results = {}
	results['omega2_grid'] = omega_grid
	results['pi'] = pis
	return results


#####################
# Command line args
#####################
est_borzoi_effect_size_file = sys.argv[1]
est_eqtl_effect_size_file = sys.argv[2]
gene_ld_summary_file = sys.argv[3]
onek_genomes_plink_filestem = sys.argv[4]
eqtl_sample_size = float(sys.argv[5])
results_file = sys.argv[6]


# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file)

# Create mapping from gene id to vector of est (standardized) eqtl effect sizes
gene_id_to_est_eqtl_effects = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file, eqtl_sample_size)

# Create mapping from gene id to ld npy
gene_id_to_ld_file = create_mapping_from_gene_id_to_ld_file(gene_ld_summary_file)



# Extract big vector (concatented across genes, variants) of standardized eqtl effects size, marginal pred borzoi effect sizes, and ldscores
gene_id_to_est_marginal_borzoi_effects, eqtl_marginal_effect_sizes, borzoi_marginal_effect_sizes, ldscores = extract_ldsc_vectors(gene_ld_summary_file, onek_genomes_plink_filestem, gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, use_adjusted_LD=False)


# Run LDSC-like regression
deltas = borzoi_marginal_effect_sizes - eqtl_marginal_effect_sizes
y = np.square(deltas) - (1.0/eqtl_sample_size)
x = ldscores
ldsc_tau2 = np.sum(x * y) / np.sum(x * x)

# Delete unecessary stuff
del y, x, eqtl_marginal_effect_sizes, borzoi_marginal_effect_sizes, ldscores



gene_id_to_est_marginal_eqtl_effects = {}
for gene_id in [*gene_id_to_est_eqtl_effects]:
	tupler = gene_id_to_est_eqtl_effects[gene_id]
	tmp_arr = []
	for ele in tupler:
		tmp_arr.append(ele[1])
	gene_id_to_est_marginal_eqtl_effects[gene_id] = np.asarray(tmp_arr)




# Run inference
resid_var = 1.0
results = run_ashr_inference(gene_id_to_est_marginal_eqtl_effects, gene_id_to_est_marginal_borzoi_effects, gene_id_to_ld_file, eqtl_sample_size, ldsc_tau2, max_iter=500, LL=50,resid_var=resid_var)



t = open(results_file,'w')
t.write('omega2\tpi\n')
for omega2, pi in zip(results["omega2_grid"], results["pi"]):
	t.write(str(omega2) + '\t' + str(pi) + '\n')
t.close()









'''
if True:
	with open('gene_id_to_est_eqtl_effects.pkl', 'wb') as f:
		pickle.dump(gene_id_to_est_eqtl_effects, f)
	with open('gene_id_to_est_marginal_borzoi_effects.pkl', 'wb') as f:
		pickle.dump(gene_id_to_est_marginal_borzoi_effects, f)


	with open('gene_id_to_est_eqtl_effects.pkl', 'rb') as f:
		gene_id_to_est_eqtl_effects = pickle.load(f)

	with open('gene_id_to_est_marginal_borzoi_effects.pkl', 'rb') as f:
		gene_id_to_est_marginal_borzoi_effects = pickle.load(f)

	ldsc_tau2 = 0.002054667819079888
	np.random.seed(1)
'''


'''
y = np.square(deltas) - (1.0/eqtl_sample_size)
x = ldscores

tau2 = np.sum(x * y) / np.sum(x * x)
t = open(results_file,'w')
t.write('borzoi_sampling_variance\n')
t.write(str(tau2) + '\n')
t.close()
'''


'''
pis = np.asarray([.4, .6])
sig_sqs = np.asarray([.1, .9])
scaler = sampling_variance/np.sum(pis*sig_sqs)
sig_sqs = sig_sqs*scaler
'''


