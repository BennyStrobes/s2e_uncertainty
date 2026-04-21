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






def initialize_data(gene_data, gtex_tissue_index, magnitude_bins):
	n_bins = len(magnitude_bins) - 1
	n_genes = gene_data.shape[0]

	# Initialize priors
	mu_priors = np.ones(n_bins)
	sig_sq_priors = np.ones(n_bins)*.1/300


	# Initialize causal effects
	betas = []
	valid_genes = []
	snps_per_gene = []

	for gene_iter in range(n_genes):

		if gene_iter > 5000:
			break

		# Relevent files for gene
		eqtl_effect_file = gene_data[gene_iter,1]
		eqtl_effect_se_file = gene_data[gene_iter, 2]
		eqtl_maf_file = gene_data[gene_iter, 3]
		borzoi_effect_file = gene_data[gene_iter, 4]
		ld_file = gene_data[gene_iter, 5]

		# Load in data for gene
		gene_eqtl_data = np.load(eqtl_effect_file)[:, gtex_tissue_index]
		gene_eqtl_data_se = np.load(eqtl_effect_se_file)[:, gtex_tissue_index]
		zeds = gene_eqtl_data/gene_eqtl_data_se

		if np.sum(np.isnan(gene_eqtl_data)) > 0:
			betas.append(np.nan)

			if np.sum(np.isnan(gene_eqtl_data) == False) != 0:
				print('gene with only some missing values')
		elif np.max(np.square(zeds)) > 80:
			betas.append(np.nan)
			print('skip')
		else:
			valid_genes.append(gene_iter)
			if np.sum(np.isnan(gene_eqtl_data)) > 0:
				print('assumption error')
				pdb.set_trace()
			betas.append(np.zeros(len(gene_eqtl_data)))

			snps_per_gene.append(len(gene_eqtl_data))


			borzois = np.load(borzoi_effect_file)
			if np.sum(np.isnan(borzois)) > 0:
				print('assumption eroror')
				pdb.set_trace()

	valid_genes = np.asarray(valid_genes)
	snps_per_gene = np.asarray(snps_per_gene)

	return betas, valid_genes, mu_priors, sig_sq_priors


def load_in_data(valid_genes, gene_data, gtex_tissue_index, betas, magnitude_bins):
	valid_genes_dicti = {}
	for ele in valid_genes:
		valid_genes_dicti[ele] = 1


	n_genes = gene_data.shape[0]
	eqtl_effect_resids = []
	eqtl_effect_ses = []
	borzoi_effects = []
	snp_annos = []
	ld_files = []

	for gene_iter in range(n_genes):

		if gene_iter not in valid_genes_dicti:
			eqtl_effect_resids.append(np.nan)
			eqtl_effect_ses.append(np.nan)
			borzoi_effects.append(np.nan)
			ld_files.append('nan')
			snp_annos.append(np.nan)
		else:
			eqtl_effect_file = gene_data[gene_iter,1]
			eqtl_effect_se_file = gene_data[gene_iter, 2]
			eqtl_maf_file = gene_data[gene_iter, 3]
			borzoi_effect_file = gene_data[gene_iter, 4]
			ld_file = gene_data[gene_iter, 5]		

			gene_eqtl_betas = np.load(eqtl_effect_file)[:, gtex_tissue_index]
			gene_eqtl_beta_ses = np.load(eqtl_effect_se_file)[:, gtex_tissue_index]
			gene_eqtl_mafs = np.load(eqtl_maf_file)[:, gtex_tissue_index]
			gene_borzoi_effects = np.load(borzoi_effect_file)[:, gtex_tissue_index]
			LD_mat = np.load(ld_file)

			zeds = gene_eqtl_betas/gene_eqtl_beta_ses


			eqtl_effect_resids.append(gene_eqtl_betas - np.dot(LD_mat, betas[gene_iter]))
			eqtl_effect_ses.append(gene_eqtl_beta_ses)

			snp_sqrt_var = np.sqrt(2.0*gene_eqtl_mafs*(1.0-gene_eqtl_mafs))
			gene_std_borzoi_effects = gene_borzoi_effects * snp_sqrt_var
			borzoi_effects.append(gene_std_borzoi_effects)
			ld_files.append(ld_file)

			n_snps = len(gene_std_borzoi_effects)
			gene_snp_anno = np.zeros(n_snps)

			for snp_iter in range(n_snps):
				bin_index = np.digitize(np.abs(gene_borzoi_effects[snp_iter]), magnitude_bins) - 1
				if bin_index == -1 or bin_index == len(magnitude_bins):
					print('assumption eroror')
					pdb.set_trace()
				gene_snp_anno[snp_iter] = bin_index
			snp_annos.append(gene_snp_anno.astype(int))


	return eqtl_effect_resids, eqtl_effect_ses, borzoi_effects, np.asarray(ld_files), snp_annos

# @njit(cache=True)
def update_causal_effect_for_single_gene(cur_gene_beta, mu_priors, sig_sq_priors, cur_eqtl_effect_resid, eqtl_effect_se, borzoi_effects, LD_mat, snp_gene_anno):
	# Compute prior mean

	prior_means = (mu_priors[snp_gene_anno])*borzoi_effects
	per_snp_sig_sq_prior = sig_sq_priors[snp_gene_anno]


	Ns = 1.0/(np.square(eqtl_effect_se)*2.4)
	ld_diag = np.diag(LD_mat)

	posterior_vars = 1.0/((1.0/per_snp_sig_sq_prior) + (Ns*ld_diag))
	posterior_sds = np.sqrt(posterior_vars)

	# Number of snps per gene
	n_snps = len(prior_means)


	# Looop through snps in gene (in random order)
	for snp_iter in np.random.permutation(n_snps):
		N = Ns[snp_iter]
		ld_vec = LD_mat[snp_iter, :]

		# Add the current SNP contribution back into the residual before sampling it
		cur_eqtl_effect_resid = cur_eqtl_effect_resid + (ld_vec*cur_gene_beta[snp_iter])

		posterior_var = posterior_vars[snp_iter]
		posterior_mean = posterior_var*((prior_means[snp_iter]/per_snp_sig_sq_prior[snp_iter]) + (N*cur_eqtl_effect_resid[snp_iter]))

		new_beta = np.random.normal(loc=posterior_mean, scale=posterior_sds[snp_iter])

		# Remove the updated SNP contribution from the residual
		cur_eqtl_effect_resid = cur_eqtl_effect_resid - (ld_vec*new_beta)
		cur_gene_beta[snp_iter] = new_beta

	return cur_gene_beta, cur_eqtl_effect_resid



def update_causal_effect_estimates(betas, valid_genes, mu_priors, sig_sq_priors, eqtl_effect_resids, eqtl_effect_ses, borzoi_effects, ld_files, snp_gene_annos):
	# Loop through genes (perform update seperately for each gene)
	for gene_iter in valid_genes:

		#####################
		# Load in LD
		gene_ld_mat = np.load(ld_files[gene_iter])

		#####################
		# Update causal effects for a single gene
		gene_beta, gene_eqtl_effect_resid = update_causal_effect_for_single_gene(betas[gene_iter], mu_priors, sig_sq_priors, eqtl_effect_resids[gene_iter], eqtl_effect_ses[gene_iter], borzoi_effects[gene_iter], gene_ld_mat, snp_gene_annos[gene_iter])
		betas[gene_iter] = gene_beta
		eqtl_effect_resids[gene_iter] = gene_eqtl_effect_resid


	return betas, eqtl_effect_resids


def update_priors_non_informative(betas, borzoi_effects, valid_genes, sig_sq_priors, snp_gene_annos,itera):
	# Weak hyperparameters
	tau_sq = 1.0
	a0 = 1e-15
	b0 = 1e-15

	all_beta = []
	all_borzoi = []
	all_anno = []
	for gene_iter in valid_genes:
		all_beta.append(betas[gene_iter])
		all_borzoi.append(borzoi_effects[gene_iter])
		all_anno.append(snp_gene_annos[gene_iter])

	all_beta = np.hstack(all_beta)
	all_borzoi = np.hstack(all_borzoi)
	all_anno = np.hstack(all_anno)

	n_bins = len(sig_sq_priors)
	mu_priors = np.zeros(n_bins)
	r_squares = np.zeros(n_bins)
	for bin_iter in range(n_bins):
		# Subset to just this bin
		indices = all_anno == bin_iter

		bin_beta = all_beta[indices]
		bin_borzoi = all_borzoi[indices]
		if len(bin_beta) == 0:
			print('errror')
			pdb.set_trace()
			mu_priors[bin_iter] = np.random.normal(loc=0.0, scale=np.sqrt(tau_sq))
			sig_sq_priors[bin_iter] = 1.0/np.random.gamma(a0, 1.0/b0)
			continue

		posterior_var = 1.0/((np.sum(np.square(bin_borzoi))/sig_sq_priors[bin_iter]) + (1.0/tau_sq))
		posterior_mean = posterior_var*(np.sum(bin_borzoi*bin_beta)/sig_sq_priors[bin_iter])
		mu_priors[bin_iter] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))

		resid = bin_beta - (bin_borzoi*mu_priors[bin_iter])

		shape_post = a0 + (len(bin_beta)/2.0)
		scale_post = b0 + (0.5*np.sum(np.square(resid)))
		sig_sq_priors[bin_iter] = 1.0/np.random.gamma(shape_post, 1.0/scale_post)

		signal = np.var(mu_priors[bin_iter]*bin_borzoi)
		#signal = np.mean(np.square(mu_priors[bin_iter]*bin_borzoi))
		r_sq = signal/(signal + sig_sq_priors[bin_iter])
		r_squares[bin_iter] = r_sq

	return mu_priors, sig_sq_priors, r_squares



def update_priors(betas, borzoi_effects, valid_genes, sig_sq_priors, snp_gene_annos,itera):
	mu_priors, sig_sq_priors, r_squares = update_priors_non_informative(betas, borzoi_effects, valid_genes, sig_sq_priors, snp_gene_annos,itera)
	return mu_priors, sig_sq_priors, r_squares

def save_posterior_samples(output_stem, posterior_samples, param_names):
	sample_output_file = output_stem + '_posterior_samples.txt'

	t = open(sample_output_file,'w')
	t.write('sample_iter\t' + '\t'.join(param_names) + '\n')
	for sample_iter in range(posterior_samples.shape[0]):
		t.write(str(sample_iter) + '\t' + '\t'.join(posterior_samples[sample_iter, :].astype(str)) + '\n')
	t.close()

	summary_output_file = output_stem + '_posterior_summary.txt'
	t = open(summary_output_file,'w')
	t.write('parameter\tmean\tstandard_error\tci_2.5\tci_97.5\n')
	for param_iter, param_name in enumerate(param_names):
		param_samples = posterior_samples[:, param_iter]
		param_mean = np.mean(param_samples)
		if len(param_samples) > 1:
			param_se = np.std(param_samples, ddof=1)/np.sqrt(len(param_samples))
		else:
			param_se = 0.0
		param_lb = np.percentile(param_samples, 2.5)
		param_ub = np.percentile(param_samples, 97.5)
		t.write(param_name + '\t' + str(param_mean) + '\t' + str(param_se) + '\t' + str(param_lb) + '\t' + str(param_ub) + '\n')
	t.close()


def create_posterior_param_names(n_mu, prior_version, n_ard):
	param_names = []
	for param_iter in range(n_mu):
		param_names.append('mu_prior_' + str(param_iter))
	param_names.append('sig_sq_prior')
	return param_names


def run_ld_corr_compete_gibbs(gene_data, gtex_tissue_index, ordered_gtex_tissues, output_stem, magnitude_bins, burn_in_iters=400, total_iters = 500):
	# Initialize model parameters
	betas, valid_genes, mu_priors, sig_sq_priors = initialize_data(gene_data, gtex_tissue_index, magnitude_bins)


	# Load in relevent data
	eqtl_effect_resids, eqtl_effect_ses, borzoi_effects, ld_files, snp_gene_annos = load_in_data(valid_genes, gene_data, gtex_tissue_index, betas, magnitude_bins)

	posterior_samples = []

	# Begin sampling
	for itera in range(total_iters):

		t1 = time.time()
		# Update causal effects
		betas, eqtl_effect_resids = update_causal_effect_estimates(betas, valid_genes, mu_priors, sig_sq_priors, eqtl_effect_resids, eqtl_effect_ses, borzoi_effects, ld_files, snp_gene_annos)
		# Update priors
		mu_priors, sig_sq_priors, r_squares = update_priors(betas, borzoi_effects, valid_genes, sig_sq_priors, snp_gene_annos, itera)
		t2 = time.time()

		print('###################', flush=True)
		print('Iteration ' + str(itera), flush=True)
		print(mu_priors, flush=True)
		print(sig_sq_priors, flush=True)
		print(r_squares, flush=True)
		print(t2-t1, flush=True)

		if itera >= burn_in_iters:
			posterior_samples.append(np.hstack((mu_priors, np.asarray([sig_sq_prior]))))

	if len(posterior_samples) > 0:
		param_names = create_posterior_param_names(len(mu_priors), prior_version, 0 if ard_prior_vars is None else len(ard_prior_vars))
		save_posterior_samples(output_stem, np.vstack(posterior_samples), param_names)

##########################
# Command line args
##########################
bayes_input_data_summary_file = sys.argv[1]
tissue_name = sys.argv[2]
borzoi_gtex_independent_target_names_file = sys.argv[3]
bayes_ld_corr_output_stem = sys.argv[4]



magnitude_bins=np.asarray([0.0, 0.01, 0.075, 0.2, 10.0])


# Load in target file
target_df = np.loadtxt(borzoi_gtex_independent_target_names_file,dtype=str,delimiter='\t')
ordered_gtex_tissues = target_df[1:,-1]
borzoi_target_indices = target_df[1:,1].astype(int)

print(ordered_gtex_tissues,flush=True)


gtex_tissue_index = np.where(ordered_gtex_tissues == tissue_name)[0][0]

# Load in gene level summary file
gene_data = np.loadtxt(bayes_input_data_summary_file, dtype=str,delimiter='\t')[1:, :]




run_ld_corr_compete_gibbs(gene_data, gtex_tissue_index, ordered_gtex_tissues, bayes_ld_corr_output_stem, magnitude_bins)



