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
		mapping[gene_id].append(effect)
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


def create_mapping_from_gene_id_to_indiviudal_expression(est_individual_expression_file, eqtl_sample_size):
	f = gzip.open(est_individual_expression_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		ind_id = data[1]
		expr_val = float(data[2])
		if gene_id not in mapping:
			mapping[gene_id] = []
		mapping[gene_id].append(expr_val)
	f.close()

	for gene_id in [*mapping]:
		mapping[gene_id] = np.asarray(mapping[gene_id])

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

			# Load in est borzoi causal effect sizes
			est_borzoi_causal_effect_sizes, tmp_variant_names2 = load_in_effects(gene_id_to_est_borzoi_effects[gene_name])
			
			# Quick error checking
			if np.array_equal(tmp_variant_names1,cis_variant_names) == False or np.array_equal(tmp_variant_names2,cis_variant_names) == False:
				print('variant mismatch error')
				pdb.set_trace()

			# Get borzoi predicted marginal effect sizes
			est_borzoi_marginal_effects = np.dot(LD, est_borzoi_causal_effect_sizes)

			gene_to_marginal_borzoi_vec[gene_name] = est_borzoi_marginal_effects


		f.close()
	return gene_to_marginal_borzoi_vec


def create_mapping_from_gene_id_to_geno_file(gene_ld_summary_file):
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
		geno_file = data[7]
		if data[1] != 'chr10' and data[1] != 'chr5' and data[1] != 'chr2' and data[1] != 'chr1' and data[1] != 'chr4':
			continue
		if ens_id in mapping:
			print('assumptioneornrornroo')
			pdb.set_trace()
		mapping[ens_id] = geno_file

	f.close()

	return mapping

def update_causal_effects_for_single_gene(
	cur_beta, resid_expr, omega_grid, log_pis, X, resid_var
):
	n_snps = len(cur_beta)
	LL = len(omega_grid)
	gene_zz = np.zeros((n_snps, LL))

	for snp_iter in np.random.permutation(n_snps):
		xj = X[:, snp_iter]

		# add back old effect
		resid_expr = resid_expr + xj * cur_beta[snp_iter]

		xjTxj = np.dot(xj, xj) / resid_var
		xjTr = np.dot(xj, resid_expr) / resid_var

		precisions = xjTxj + 1.0 / omega_grid
		post_vars = 1.0 / precisions
		post_means = xjTr / precisions

		log_weights = (
			log_pis
			- 0.5 * np.log(omega_grid)
			+ 0.5 * np.log(post_vars)
			+ 0.5 * (post_means ** 2) / post_vars
		)
		log_weights -= logsumexp(log_weights)
		resp = np.exp(log_weights)

		new_beta = np.sum(resp * post_means)
		cur_beta[snp_iter] = new_beta
		gene_zz[snp_iter, :] = resp

		# remove new effect
		resid_expr = resid_expr - xj * new_beta

	return cur_beta, gene_zz, resid_expr


@njit
def _logsumexp_numba(arr):
	m = np.max(arr)
	total = 0.0
	for i in range(arr.shape[0]):
		total += np.exp(arr[i] - m)
	return m + np.log(total)


@njit
def update_causal_effects_for_single_gene_numba(
	cur_beta,
	resid_expr,
	inv_omega_grid,
	neg_half_log_omega_grid,
	log_pis,
	X,
	resid_var,
	snp_order
):
	n_snps = cur_beta.shape[0]
	LL = inv_omega_grid.shape[0]
	gene_zz = np.zeros((n_snps, LL))

	precisions = np.empty(LL)
	post_vars = np.empty(LL)
	post_means = np.empty(LL)
	log_weights = np.empty(LL)
	resp = np.empty(LL)

	for ii in range(n_snps):
		snp_iter = snp_order[ii]
		xj = X[:, snp_iter]

		old_beta = cur_beta[snp_iter]

		# add back old effect
		resid_expr = resid_expr + xj * old_beta

		xjTxj = np.dot(xj, xj) / resid_var
		xjTr = np.dot(xj, resid_expr) / resid_var

		for ll in range(LL):
			precisions[ll] = xjTxj + inv_omega_grid[ll]
			post_vars[ll] = 1.0 / precisions[ll]
			post_means[ll] = xjTr * post_vars[ll]
			log_weights[ll] = (
				log_pis[ll]
				+ neg_half_log_omega_grid[ll]
				+ 0.5 * np.log(post_vars[ll])
				+ 0.5 * (post_means[ll] * post_means[ll]) / post_vars[ll]
			)

		lse = _logsumexp_numba(log_weights)

		new_beta = 0.0
		for ll in range(LL):
			resp[ll] = np.exp(log_weights[ll] - lse)
			gene_zz[snp_iter, ll] = resp[ll]
			new_beta += resp[ll] * post_means[ll]

		cur_beta[snp_iter] = new_beta

		# remove new effect
		resid_expr = resid_expr - xj * new_beta

	return cur_beta, gene_zz, resid_expr



def run_ashr_inference(gene_id_to_ind_expr, gene_id_to_est_borzoi_effects, gene_id_to_geno_file, eqtl_sample_size, ldsc_tau2, LL=20, max_iter=500, resid_var_init=1.0):
	###############################
	# Model initialization
	###############################
	gene_names = np.sort(np.asarray([*gene_id_to_geno_file]))
	# Initialize causal effects
	gene_id_to_resid_expr = {}
	gene_id_to_geno_var = {}
	gene_id_to_beta = {}
	gene_id_to_zz = {}
	for gene_name in gene_names:
		gene_geno = np.load(gene_id_to_geno_file[gene_name])
		resid_expr = np.dot(gene_geno, np.copy(gene_id_to_est_borzoi_effects[gene_name])) - np.copy(gene_id_to_ind_expr[gene_name])
		gene_id_to_resid_expr[gene_name] = resid_expr
		gene_id_to_geno_var[gene_name] = np.diag(np.dot(np.transpose(gene_geno), gene_geno))
		gene_id_to_beta[gene_name] = np.zeros(len(gene_id_to_est_borzoi_effects[gene_name]))
		gene_id_to_zz[gene_name] = np.zeros((len(gene_id_to_est_borzoi_effects[gene_name]),LL))
	# Set up grid
	sigma_grid = np.concatenate([
		np.array([0.0]),
		np.exp(np.linspace(np.log(1e-3), np.log(1.0), LL-1))
	])
	omega_grid = sigma_grid**2 + 1e-10
	inv_omega_grid = 1.0 / omega_grid
	neg_half_log_omega_grid = -0.5 * np.log(omega_grid)
	pis = np.zeros(LL) + .0005
	pis[np.argmin(np.abs(omega_grid - ldsc_tau2)):] = 1
	pis = pis/np.sum(pis)
	log_pis = np.log(pis)
	resid_vars = np.zeros(len(gene_names)) + resid_var_init


	###############################
	# Iterative expectation-maximization algorithm
	###############################
	for itera in range(max_iter):
		print(itera)
		print(pis)
		print(np.sum(pis*omega_grid))
		#print(np.sort(resid_vars))
		t1 = time.time()
		###########################
		# E-step
		###########################
		# Loop through genes
		for gene_iter, gene_name in enumerate(gene_names):
			#gene_beta, gene_zz, gene_resid_expr = update_causal_effects_for_single_gene(gene_id_to_beta[gene_name], gene_id_to_resid_expr[gene_name], omega_grid, log_pis, np.load(gene_id_to_geno_file[gene_name]), resid_vars[gene_iter])
			order = np.random.permutation(gene_id_to_beta[gene_name].shape[0])

			gene_beta, gene_zz, gene_resid_expr = update_causal_effects_for_single_gene_numba(
				gene_id_to_beta[gene_name],
				gene_id_to_resid_expr[gene_name],
				inv_omega_grid,
				neg_half_log_omega_grid,
				log_pis,
				np.load(gene_id_to_geno_file[gene_name]).astype(np.float64),
				resid_vars[gene_iter],
				order
			)

			gene_id_to_beta[gene_name] = gene_beta
			gene_id_to_resid_expr[gene_name] = gene_resid_expr
			gene_id_to_zz[gene_name] = gene_zz

		###########################
		# M-step
		###########################
		if itera >= 5:
			# Pi Update
			counts = np.zeros(LL) + 1e-10
			for gene_name in gene_names:
				counts = counts + np.sum(gene_id_to_zz[gene_name],axis=0)
			pis = counts/np.sum(counts)
			log_pis = np.log(pis)

			# Resid var update
			for gene_iter, gene_name in enumerate(gene_names):
				resid = gene_id_to_resid_expr[gene_name]
				resid_vars[gene_iter] = np.mean(resid**2)

		t2 = time.time()
		print(t2-t1)
	
	# Open output object
	results = {}
	results['omega2_grid'] = omega_grid
	results['pi'] = pis
	return results



def update_causal_effects_for_single_gene_gibbs(
	cur_beta, resid_expr, omega_grid, log_pis, X, resid_var
):
	n_snps = len(cur_beta)
	LL = len(omega_grid)
	gene_zz = np.zeros(n_snps).astype(int)

	for snp_iter in np.random.permutation(n_snps):
		xj = X[:, snp_iter]

		# add back old effect
		resid_expr = resid_expr + xj * cur_beta[snp_iter]

		xjTxj = np.dot(xj, xj) / resid_var
		xjTr = np.dot(xj, resid_expr) / resid_var

		precisions = xjTxj + 1.0 / omega_grid
		post_vars = 1.0 / precisions
		post_means = xjTr / precisions

		log_weights = (
			log_pis
			- 0.5 * np.log(omega_grid)
			+ 0.5 * np.log(post_vars)
			+ 0.5 * (post_means ** 2) / post_vars
		)
		log_weights -= logsumexp(log_weights)
		resp = np.exp(log_weights)

		# Sample ZZ
		idx = np.random.choice(len(resp), p=resp)
		gene_zz[snp_iter] = idx

		# Sample beta
		new_beta = np.random.normal(loc=post_means[idx], scale=np.sqrt(post_vars[idx]))
		cur_beta[snp_iter] = new_beta

		# remove new effect
		resid_expr = resid_expr - xj * new_beta

	return cur_beta, gene_zz, resid_expr



@njit
def _sample_from_logweights(log_weights):
	# sample index proportional to exp(log_weights)
	max_logw = log_weights[0]
	for i in range(1, log_weights.shape[0]):
		if log_weights[i] > max_logw:
			max_logw = log_weights[i]

	total = 0.0
	for i in range(log_weights.shape[0]):
		total += np.exp(log_weights[i] - max_logw)

	u = np.random.random() * total
	csum = 0.0
	for i in range(log_weights.shape[0]):
		csum += np.exp(log_weights[i] - max_logw)
		if u <= csum:
			return i

	return log_weights.shape[0] - 1


@njit
def _fisher_yates_shuffle(arr):
	for i in range(arr.shape[0] - 1, 0, -1):
		j = np.random.randint(0, i + 1)
		tmp = arr[i]
		arr[i] = arr[j]
		arr[j] = tmp


@njit
def update_causal_effects_for_single_gene_gibbs_numba(
	cur_beta,
	resid_expr,
	omega_grid,
	log_pis,
	X,
	resid_var
):
	n_samples, n_snps = X.shape
	LL = omega_grid.shape[0]

	gene_zz = np.empty(n_snps, dtype=np.int64)

	inv_omega_grid = 1.0 / omega_grid
	neg_half_log_omega_grid = -0.5 * np.log(omega_grid)

	snp_order = np.arange(n_snps)
	_fisher_yates_shuffle(snp_order)

	log_weights = np.empty(LL)
	post_means = np.empty(LL)
	post_vars = np.empty(LL)

	for ii in range(n_snps):
		snp_iter = snp_order[ii]

		# add back old effect
		old_beta = cur_beta[snp_iter]
		for n in range(n_samples):
			resid_expr[n] += X[n, snp_iter] * old_beta

		# compute xjTxj and xjTr
		xjTxj = 0.0
		xjTr = 0.0
		for n in range(n_samples):
			xval = X[n, snp_iter]
			xjTxj += xval * xval
			xjTr += xval * resid_expr[n]
		xjTxj = xjTxj / resid_var
		xjTr = xjTr / resid_var

		# compute component-specific posterior quantities
		for ll in range(LL):
			prec = xjTxj + inv_omega_grid[ll]
			pv = 1.0 / prec
			pm = xjTr * pv

			post_vars[ll] = pv
			post_means[ll] = pm
			log_weights[ll] = (
				log_pis[ll]
				+ neg_half_log_omega_grid[ll]
				+ 0.5 * np.log(pv)
				+ 0.5 * (pm * pm) / pv
			)

		# sample mixture component
		idx = _sample_from_logweights(log_weights)
		gene_zz[snp_iter] = idx

		# sample beta
		new_beta = np.random.normal(post_means[idx], np.sqrt(post_vars[idx]))
		cur_beta[snp_iter] = new_beta

		# remove new effect
		for n in range(n_samples):
			resid_expr[n] -= X[n, snp_iter] * new_beta

	return cur_beta, gene_zz, resid_expr


def run_ashr_inference_gibbs(gene_id_to_ind_expr, gene_id_to_est_borzoi_effects, gene_id_to_geno_file, eqtl_sample_size, ldsc_tau2, LL=20, burn_in_iters=300, max_iter=500, resid_var_init=1.0, a0=1e-8, b0=1e-8):
	###############################
	# Model initialization
	###############################
	gene_names = np.sort(np.asarray([*gene_id_to_geno_file]))
	# Initialize causal effects
	gene_id_to_resid_expr = {}
	gene_id_to_beta = {}
	gene_id_to_zz = {}
	for gene_name in gene_names:
		gene_geno = np.load(gene_id_to_geno_file[gene_name])
		resid_expr = np.dot(gene_geno, np.copy(gene_id_to_est_borzoi_effects[gene_name])) - np.copy(gene_id_to_ind_expr[gene_name])
		gene_id_to_resid_expr[gene_name] = resid_expr
		gene_id_to_beta[gene_name] = np.zeros(len(gene_id_to_est_borzoi_effects[gene_name]))
		gene_id_to_zz[gene_name] = np.zeros(len(gene_id_to_est_borzoi_effects[gene_name]), dtype=int)
	# Set up grid
	sigma_grid = np.concatenate([
		np.array([0.0]),
		np.exp(np.linspace(np.log(1e-3), np.log(1.0), LL-1))
	])
	omega_grid = sigma_grid**2 + 1e-10
	inv_omega_grid = 1.0 / omega_grid
	neg_half_log_omega_grid = -0.5 * np.log(omega_grid)
	pis = np.zeros(LL) + .0005
	pis[np.argmin(np.abs(omega_grid - ldsc_tau2)):] = 1
	pis = pis/np.sum(pis)
	log_pis = np.log(pis)
	resid_vars = np.zeros(len(gene_names)) + resid_var_init

	# Sampled pis
	sampled_pis = []

	###############################
	# Iterative expectation-maximization algorithm
	###############################
	for itera in range(max_iter):
		print(itera)
		print(pis)
		print(np.sum(pis*omega_grid))
		#print(np.sort(resid_vars))
		t1 = time.time()
		###########################
		# E-step
		###########################
		# Loop through genes
		for gene_iter, gene_name in enumerate(gene_names):
			gene_beta, gene_zz, gene_resid_expr = update_causal_effects_for_single_gene_gibbs_numba(gene_id_to_beta[gene_name], gene_id_to_resid_expr[gene_name], omega_grid, log_pis, np.load(gene_id_to_geno_file[gene_name]), resid_vars[gene_iter])

			gene_id_to_beta[gene_name] = gene_beta
			gene_id_to_resid_expr[gene_name] = gene_resid_expr
			gene_id_to_zz[gene_name] = gene_zz

		###########################
		# M-step
		###########################
		if itera >= 5:
			# Pi Update
			counts = np.zeros(LL) + 1
			for gene_name in gene_names:
				zz = gene_id_to_zz[gene_name]  # shape: (n_snps,), values in [0, LL-1]
				counts += np.bincount(zz, minlength=LL)
			pis = np.random.dirichlet(counts)
			log_pis = np.log(pis)

			# Resid var update
			for gene_iter, gene_name in enumerate(gene_names):
				resid = gene_id_to_resid_expr[gene_name]
				n = len(resid)
				ss = np.sum(resid**2)

				a_post = a0 + n / 2.0
				b_post = b0 + 0.5 * ss

				# sample from inverse-gamma via 1 / Gamma
				sigma2 = 1.0 / np.random.gamma(shape=a_post, scale=1.0 / b_post)
				resid_vars[gene_iter] = sigma2

		if itera >= burn_in_iters:
			sampled_pis.append(np.copy(pis))


		t2 = time.time()
		print(t2-t1)
	
	# Open output object
	sampled_pis = np.asarray(sampled_pis)
	results = {}
	results['omega2_grid'] = omega_grid
	results['pi'] = np.mean(sampled_pis,axis=0)
	return results







#####################
# Command line args
#####################
est_borzoi_effect_size_file = sys.argv[1]
est_individual_expression_file = sys.argv[2]
gene_ld_summary_file = sys.argv[3]
onek_genomes_plink_filestem = sys.argv[4]
eqtl_sample_size = float(sys.argv[5])
results_file = sys.argv[6]


# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file)

# Create mapping from gene id to vector of est (standardized) eqtl effect sizes
gene_id_to_ind_expr = create_mapping_from_gene_id_to_indiviudal_expression(est_individual_expression_file, eqtl_sample_size)

# Create mapping from gene id to ld npy
gene_id_to_geno_file = create_mapping_from_gene_id_to_geno_file(gene_ld_summary_file)


# Run LDSC-like regression
'''
deltas = borzoi_marginal_effect_sizes - eqtl_marginal_effect_sizes
y = np.square(deltas) - (1.0/eqtl_sample_size)
x = ldscores
'''
ldsc_tau2 = 0.002

'''
# Delete unecessary stuff
del y, x, eqtl_marginal_effect_sizes, borzoi_marginal_effect_sizes, ldscores
'''



# Run inference
resid_var_init = 1.0
#results = run_ashr_inference(gene_id_to_ind_expr, gene_id_to_est_borzoi_effects, gene_id_to_geno_file, eqtl_sample_size, ldsc_tau2, max_iter=200, LL=20,resid_var_init=resid_var_init)



# Run inference
resid_var_init = 1.0
results = run_ashr_inference_gibbs(gene_id_to_ind_expr, gene_id_to_est_borzoi_effects, gene_id_to_geno_file, eqtl_sample_size, ldsc_tau2, max_iter=5000, burn_in_iters=4600, LL=20,resid_var_init=resid_var_init)





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


