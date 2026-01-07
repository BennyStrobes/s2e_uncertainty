import numpy as np
import os
import sys
import pdb
import gzip
import time

def extract_gene_level_sumstat_data_structure(sardinia_raw_sumstat_file, chrom_string):
	chrom_word = chrom_string + '\t'
	dicti = {}
	head_count = 0
	f = gzip.open(sardinia_raw_sumstat_file, 'rt')
	for line in f:
		if head_count == 0:
			head_count = head_count + 1
			header_line = line.rstrip()
			continue
		if line.startswith(chrom_word) == False:
			continue
		line = line.rstrip()
		data = line.split('\t')
		'''
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_chrom = data[0]
		if line_chrom != chrom_string:
			continue
		'''

		gene_id = data[5]
		if gene_id not in dicti:
			dicti[gene_id] = []
		dicti[gene_id].append(data)
	f.close()
	return dicti, header_line


def wakefield_logabf(beta_hat, se, W):
	beta_hat = np.asarray(beta_hat, float)
	se = np.asarray(se, float)
	V = se**2
	return 0.5*(np.log(V) - np.log(V+W)) + 0.5*beta_hat**2 * (1.0/V - 1.0/(V+W))

def single_causal_pips(beta_hat, se, W=0.01, prior=None):
	log_abf = wakefield_logabf(beta_hat, se, W)

	m = len(log_abf)
	if prior is None:
		prior = np.ones(m) / m
	else:
		prior = np.asarray(prior, float)
		prior = prior / prior.sum()

	log_post = np.log(prior) + log_abf
	# stable normalization
	a = np.max(log_post)
	w = np.exp(log_post - a)
	pip = w / w.sum()
	return pip, log_abf

def credible_set(pip, alpha=0.95):
	pip = np.asarray(pip, float)
	order = np.argsort(-pip)
	cs = []
	s = 0.0
	for j in order:
		cs.append(j)
		s += pip[j]
		if s >= alpha:
			break
	return np.array(cs), pip[cs], s


def extract_effect_sizes_and_standard_errors_for_this_gene(per_gene_info):
	effect_sizes = []
	ses = []

	for ele in per_gene_info:
		effect_sizes.append(ele[-3])
		ses.append(ele[-2])



	return np.asarray(effect_sizes).astype(float), np.asarray(ses).astype(float)




#######################
# Command line args
#######################
sardinia_raw_sumstat_file = sys.argv[1]
fm_eqtl_sumstats_file = sys.argv[2] 
pip_threshold = float(sys.argv[3])

t = open(fm_eqtl_sumstats_file,'w')

# Loop through chromosomes
for chrom_num in range(1,23):
	print(chrom_num)
	# Extract gene level data structure where each key is a gene and values are all snps
	# But only on this chromosome
	chrom_gene_info, header_line = extract_gene_level_sumstat_data_structure(sardinia_raw_sumstat_file, str(chrom_num))
	if chrom_num == 1:
		t.write(header_line + '\n')
	# Extract names of genes
	gene_ids = np.asarray([*chrom_gene_info])

	# Loop through genes
	for gene_id in gene_ids:
		# Extract relevent info
		per_gene_info = chrom_gene_info[gene_id]

		# Extract effect sizes and standard errors for this gene
		effect_sizes, ses = extract_effect_sizes_and_standard_errors_for_this_gene(per_gene_info)
		if len(effect_sizes) < 10:
			continue
		# Run single causal variant fine-mapping
		pips, log_abfs = single_causal_pips(effect_sizes,ses)

		if np.max(pips) < pip_threshold:
			continue

		best_index = np.argmax(pips)
		best_data = np.asarray(per_gene_info[best_index])
		t.write('\t'.join(best_data) + '\n')
t.close()