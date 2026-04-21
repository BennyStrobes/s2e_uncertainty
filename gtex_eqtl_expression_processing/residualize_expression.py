import numpy as np
import os
import sys
import pdb
import gzip




def load_in_covariate_data(input_covariate_file):
	f = open(input_covariate_file)
	sample_names = []
	cov_mat = []
	cov_names = []

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			sample_names = np.asarray(data[1:])
			continue
		cov_names.append(data[0])
		cov_mat.append(np.asarray(data[1:]).astype(float))
	f.close()

	sample_names = np.asarray(sample_names)
	cov_names = np.asarray(cov_names)
	cov_mat = np.asarray(cov_mat)

	return cov_mat, sample_names, cov_names



def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf):
	valid_chroms = {}
	for chrom_num in range(1,23):
		valid_chroms['chr' + str(chrom_num)] = 1


	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in valid_chroms:
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = 1

	f.close()

	return dicti



#########################
# Command line args
#########################
input_expression_file = sys.argv[1]
input_covariate_file = sys.argv[2]
residual_expression_file = sys.argv[3]
gtex_v10_pc_genes_gtf = sys.argv[4]


protein_coding_genes_dictionary = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)


###############################
# First load in covariate file
cov_mat, ordered_cov_sample_names, ordered_cov_names = load_in_covariate_data(input_covariate_file)
cov_design = np.transpose(cov_mat)
intercept = np.ones((cov_design.shape[0], 1))
cov_design = np.hstack((intercept, cov_design))




###############################
# Second, residualize out expression
f = gzip.open(input_expression_file, 'rt')
t = open(residual_expression_file,'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		ordered_expr_sample_names = np.asarray(data[4:])
		if np.array_equal(ordered_expr_sample_names, ordered_cov_sample_names) == False:
			print('assumption eroror')
			pdb.set_trace()
		t.write('\t'.join(data) + '\n')
		continue

	gene_id = data[3]

	if gene_id.startswith('ENSG') == False:
		print('assumption erororo')
		pdb.set_trace()
	if len(gene_id.split('.')) != 2:
		print('assumptioneorronror')
		pdb.set_trace()

	if gene_id.split('.')[0] not in protein_coding_genes_dictionary:
		continue

	expr_vec = np.asarray(data[4:]).astype(float)
	beta, _, _, _ = np.linalg.lstsq(cov_design, expr_vec, rcond=None)
	pred_expr = np.dot(cov_design, beta)
	resid_expr = expr_vec - pred_expr

	resid_std = np.std(resid_expr)
	if resid_std == 0.0:
		print('errorororoororororo')
		standardized_resid_expr = np.zeros(len(resid_expr))
	else:
		standardized_resid_expr = (resid_expr - np.mean(resid_expr))/resid_std
	t.write('\t'.join(data[:4]) + '\t' + '\t'.join(standardized_resid_expr.astype(str)) + '\n')
f.close()
t.close()
