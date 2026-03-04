import numpy as np
import os
import sys
import pdb
from scipy import sparse
import statsmodels.api as sm


def simulate_data(plof_variant_to_gene_effect, regulatory_variant_effects, regulatory_variants_per_class_in_gene, sample_size, n_total_genes, n_causal_genes, rvat_h2, n_plof_variants_per_gene, n_regulatory_variant_classes, n_detected):
	total_variants_per_gene = n_plof_variants_per_gene + (n_regulatory_variant_classes*regulatory_variants_per_class_in_gene)
	gene_info = []
	phenotype = np.zeros(sample_size)

	causal_genes = np.random.choice(n_total_genes, n_causal_genes, replace=False)
	causal_genes_dicti = {}
	for indexer in causal_genes:
		causal_genes_dicti[indexer] = 1

	rare_indis = np.arange(1,int(sample_size*.0001) + 1)

	annotation_to_gene_effects = np.hstack(([plof_variant_to_gene_effect],regulatory_variant_effects ))


	true_alphas = []
	observed_anno_mats = []

	for gene_iter in range(n_total_genes):

		variant_iter = 0

		
		# pLOF variants
		true_anno_weighted_dosage = np.zeros((sample_size, n_regulatory_variant_classes + 1))
		observed_anno_weighted_dosage = np.zeros((sample_size, n_regulatory_variant_classes + 1))

		# PLOF
		anno_counter = 0
		for ii in range(n_plof_variants_per_gene):

			n_indis = np.random.choice(rare_indis)
			dosage = np.zeros(sample_size)
			dosage[np.random.choice(np.arange(sample_size), n_indis)] = 1.0

			true_anno_weighted_dosage[:,0] = true_anno_weighted_dosage[:,0] + (-dosage)
			observed_anno_weighted_dosage[:,0] = observed_anno_weighted_dosage[:,0] + (-dosage)

		anno_counter = anno_counter + 1

		for varant_class in range(n_regulatory_variant_classes):

			detected_variants = np.random.choice(np.arange(regulatory_variants_per_class_in_gene), size=n_detected, replace=False)
			detected_variant_dicti = {}
			for ele in detected_variants:
				detected_variant_dicti[ele] = 1

			for var_iter in range(regulatory_variants_per_class_in_gene):

				n_indis = np.random.choice(rare_indis)
				dosage = np.zeros(sample_size)
				dosage[np.random.choice(np.arange(sample_size), n_indis)] = 1.0

				anno_effect = np.random.normal(loc=0,scale=1)
				true_anno_weighted_dosage[:, anno_counter] = true_anno_weighted_dosage[:, anno_counter] + (dosage*anno_effect)

				observed_anno_effect = np.random.normal(loc=anno_effect,scale=np.sqrt(.5))
				if var_iter in detected_variant_dicti:
					observed_anno_weighted_dosage[:,anno_counter] = observed_anno_weighted_dosage[:,anno_counter] + (dosage*observed_anno_effect)

			# Add extra variant
			n_indis = np.random.choice(rare_indis)
			dosage = np.zeros(sample_size)
			dosage[np.random.choice(np.arange(sample_size), n_indis)] = 1.0
			anno_effect = np.random.normal(loc=0,scale=1)
			observed_anno_weighted_dosage[:,anno_counter] = observed_anno_weighted_dosage[:,anno_counter] + (dosage*anno_effect)

			anno_counter = anno_counter + 1
		

		genetic_gene = np.dot(true_anno_weighted_dosage, annotation_to_gene_effects)
		std_genetic_gene = (genetic_gene - np.mean(genetic_gene))/np.std(genetic_gene)

		alpha = 0.0
		if gene_iter in causal_genes:
			alpha = np.random.normal(loc=0,scale=np.sqrt(rvat_h2/n_causal_genes))

		phenotype = phenotype + alpha*std_genetic_gene
		true_alphas.append(alpha)

		sparse_observed_anno_weighted_dosage = sparse.csr_matrix(observed_anno_weighted_dosage)

		observed_anno_mats.append(sparse_observed_anno_weighted_dosage)

	full_phenotype = np.random.normal(loc=phenotype, scale=np.sqrt(1.0-rvat_h2))


	true_alphas = np.asarray(true_alphas)

	return full_phenotype, observed_anno_mats, true_alphas


def learn_alpha(gene_mat, weights, sim_pheno):
	genetic_gene = np.dot(gene_mat.toarray(), weights)
	X = sm.add_constant(genetic_gene)  # adds intercept

	model = sm.OLS(sim_pheno, X).fit()

	beta = model.params[1]
	se = model.bse[1]
	pval = model.pvalues[1]	

	return beta, se, pval

def iterative_algorithm(train_gene_mats, sim_pheno, max_iter=10):
	# Initialize weights to 1.0 for plof, 0 else
	weights = np.zeros(train_gene_mats[0].shape[1])
	weights[0] = 1.0
	n_genes = len(train_gene_mats)

	components2 = np.zeros(train_gene_mats[0].shape)

	for itera in range(max_iter):
		print(itera)

		# Step 1
		alphas = []
		for gene_iter in range(n_genes):
			alpha, alpha_se, alpha_p = learn_alpha(train_gene_mats[gene_iter], weights, sim_pheno)
			alphas.append(alpha)
		alphas = np.asarray(alphas)


		# Step 2
		# zero out
		components2 = components2*0.0
		for gene_iter in range(n_genes):
			components2 = components2 + (train_gene_mats[gene_iter].toarray())*alphas[gene_iter]
		model = sm.OLS(sim_pheno, sm.add_constant(components2[:,1:])).fit()
		weights[1:] = model.params[1:]


	return -weights






######################
# Command line args
######################
sim_number = int(sys.argv[1])
sample_size = int(sys.argv[2])
n_detected = int(sys.argv[3])
output_stem = sys.argv[4]

n_total_genes = 500
n_causal_genes = 100


rvat_h2 = .03
n_plof_variants_per_gene = 5
n_regulatory_variant_classes = 10
regulatory_variants_per_class_in_gene=20
n_regulatory_variants_per_gene = np.repeat(regulatory_variants_per_class_in_gene,n_regulatory_variant_classes)
n_folds = 10


output_file = output_stem + '_sim_summary.txt'
t = open(output_file,'w')
t.write('gene_index\tmethod\tsim_alpha\test_alpha\test_alpha_se\test_alpha_p\n')

plof_variant_to_gene_effect = -np.abs(np.random.normal(0, scale=np.sqrt(.2/n_plof_variants_per_gene)))
regulatory_variant_effects = np.random.normal(0, scale=np.sqrt(.8/regulatory_variants_per_class_in_gene), size=n_regulatory_variant_classes)

sim_pheno, gene_mats, sim_alphas = simulate_data(plof_variant_to_gene_effect, regulatory_variant_effects, regulatory_variants_per_class_in_gene, sample_size, n_total_genes, n_causal_genes, rvat_h2, n_plof_variants_per_gene, n_regulatory_variant_classes, n_detected)


folds = np.array_split(np.arange(n_total_genes), n_folds)

burden_weights = np.zeros(n_regulatory_variant_classes + 1)
burden_weights[0] = -1

# Loop through folds
for fold_iter, held_out_genes in enumerate(folds):

	# Held out causal effects
	held_out_alphas = sim_alphas[held_out_genes]
	held_out_gene_mats = [gene_mats[i] for i in held_out_genes]

	held_out_set = set(held_out_genes)
	train_gene_mats = [
		gene_mats[i] for i in range(len(gene_mats))
		if i not in held_out_set
	]

	learned_anno_weights = iterative_algorithm(train_gene_mats, sim_pheno)


	for gene_iter, full_gene_number in enumerate(held_out_genes):


		alpha, alpha_se, alpha_p = learn_alpha(held_out_gene_mats[gene_iter], learned_anno_weights, sim_pheno)

		alpha_burden, alpha_burden_se, alpha_burden_p = learn_alpha(held_out_gene_mats[gene_iter], burden_weights, sim_pheno)

		t.write(str(full_gene_number) + '\t' + 'dRVAT' + '\t' + str(held_out_alphas[gene_iter]) + '\t' + str(alpha) + '\t' + str(alpha_se) + '\t' + str(alpha_p) + '\n')
		t.write(str(full_gene_number) + '\t' + 'burden' + '\t' + str(held_out_alphas[gene_iter]) + '\t' + str(alpha_burden) + '\t' + str(alpha_burden_se) + '\t' + str(alpha_burden_p) + '\n')


t.close()


print(output_file)
