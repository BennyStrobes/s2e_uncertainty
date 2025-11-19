import numpy as np
import os
import sys
import pdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import normaltest, spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler







def remove_low_expression_genes(tpm_mat, tpm_lower_thresh=0.1, frac_low_samples_thresh=0.2):

	n_low_samples = np.sum(tpm_mat <= tpm_lower_thresh,axis=1)

	frac_low_samples = n_low_samples/tpm_mat.shape[1]

	valid_gene_indices_boolean = frac_low_samples <= frac_low_samples_thresh

	return valid_gene_indices_boolean

def regress_out_n_pcs(expr_mat, n_pcs=10) -> pd.DataFrame:
	X = np.transpose(np.copy(expr_mat)) # (n_samples, n_genes)
	scaler = StandardScaler()
	X_scaled = scaler.fit_transform(X)

	pca = PCA(n_components=10)
	pcs = pca.fit_transform(X_scaled)
	X_corrected = X_scaled - pcs @ pca.components_
	X_corrected = scaler.inverse_transform(X_corrected)

	return np.transpose(X_corrected)


def print_expression_to_output(output_file, log_tpm_mat, gene_names_mat, header):
	t = open(output_file,'w')

	t.write('\t'.join(header) + '\n')

	for row_iter in range(log_tpm_mat.shape[0]):
		t.write('\t'.join(gene_names_mat[row_iter, :]) + '\t')
		t.write('\t'.join(log_tpm_mat[row_iter,:].astype(str)) + '\n')
	t.close()
	return



########################
# Command line args
########################
input_expression_file = sys.argv[1] # TPM expression file
output_expression_file_stem = sys.argv[2]
tpm_lower_thresh = float(sys.argv[3])
max_prop_sample_missing = float(sys.argv[4])


pseudocount = 1e-2
n_pcs = 10

# Most of processing taken from: https://github.com/ni-lab/finetuning-enformer/blob/main/process_geuvadis_data/log_tpm/create_corrected_log_tpm.ipynb


###########################
# Load in expression data
############################
# Load the data directly from the gzipped file
raw_tpm_mat = np.loadtxt(input_expression_file, dtype=str,delimiter='\t',skiprows=2)[:,1:]
# Parse expression fields
gene_names_mat = raw_tpm_mat[1:, :2]
header = raw_tpm_mat[0,:]
tpm_mat = raw_tpm_mat[1:, 2:].astype(float)



###########################
# Remove genes with low expression across samples
############################
valid_genes_boolean = remove_low_expression_genes(tpm_mat, tpm_lower_thresh=tpm_lower_thresh, frac_low_samples_thresh=max_prop_sample_missing)
print(np.sum(valid_genes_boolean))
# filter tpm mat and gene_names mat
gene_names_mat = gene_names_mat[valid_genes_boolean, :]
tpm_mat = tpm_mat[valid_genes_boolean, :]

###########################
# log transform
############################
log_tpm_mat = np.log(tpm_mat + pseudocount)


###########################
# Regress out n_pcs
############################
corrected_log_tpm_mat = regress_out_n_pcs(log_tpm_mat, n_pcs=n_pcs)


#############################
# Print to output file
#############################
uncorrected_output_file = output_expression_file_stem + '.txt'
print_expression_to_output(uncorrected_output_file, log_tpm_mat, gene_names_mat, header)
print(uncorrected_output_file)

corrected_output_file = output_expression_file_stem + '_pc_corrected.txt'
print_expression_to_output(corrected_output_file, corrected_log_tpm_mat, gene_names_mat, header)



