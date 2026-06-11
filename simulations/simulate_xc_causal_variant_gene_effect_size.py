import numpy as np
import sys
import gzip





def load_in_variants(variant_file):
	f2 = open(variant_file)
	arr = []
	a0_arr = []
	a1_arr = []
	chr_arr = []
	pos_arr = []
	for line in f2:
		line = line.rstrip()
		data = line.split('\t')
		arr.append(data[1])
		a0_arr.append(data[4])
		a1_arr.append(data[5])
		chr_arr.append(data[0])
		pos_arr.append(data[3])
	f2.close()
	return np.asarray(arr), np.asarray(a0_arr), np.asarray(a1_arr), np.asarray(chr_arr), np.asarray(pos_arr)


def simulate_cross_context_causal_effect_sizes(tot_variants_in_gene, n_shared_causal_variants, n_tissue_specific_causal_variants, cis_snp_h2s):

	n_total_causal_variants = n_shared_causal_variants + (2*n_tissue_specific_causal_variants)
	if tot_variants_in_gene < n_total_causal_variants:
		raise ValueError('Need at least ' + str(n_total_causal_variants) + ' variants per gene to sample ' + str(n_shared_causal_variants) + ' shared and ' + str(n_tissue_specific_causal_variants) + ' tissue-specific variants per tissue.')

	t1_causal_effects = np.zeros(tot_variants_in_gene)
	t2_causal_effects = np.zeros(tot_variants_in_gene)

	causal_variant_indices = np.random.choice(np.arange(tot_variants_in_gene), size=n_total_causal_variants, replace=False)
	shared_causal_variant_indices = causal_variant_indices[:n_shared_causal_variants]
	t1_specific_causal_variant_indices = causal_variant_indices[n_shared_causal_variants:(n_shared_causal_variants+n_tissue_specific_causal_variants)]
	t2_specific_causal_variant_indices = causal_variant_indices[(n_shared_causal_variants+n_tissue_specific_causal_variants):]

	for index in shared_causal_variant_indices:
		shared_effect = np.random.normal(loc=0, scale=np.sqrt(np.random.choice(cis_snp_h2s)))
		t1_causal_effects[index] = shared_effect
		t2_causal_effects[index] = shared_effect

	for index in t1_specific_causal_variant_indices:
		t1_causal_effects[index] = np.random.normal(loc=0, scale=np.sqrt(np.random.choice(cis_snp_h2s)))

	for index in t2_specific_causal_variant_indices:
		t2_causal_effects[index] = np.random.normal(loc=0, scale=np.sqrt(np.random.choice(cis_snp_h2s)))

	return t1_causal_effects, t2_causal_effects




########################
# Command line args
########################
simulation_iter = int(sys.argv[1])
gene_ld_summary_file = sys.argv[2]
t1_causal_variant_gene_effect_size_file = sys.argv[3]
t2_causal_variant_gene_effect_size_file = sys.argv[4]


# Potential n causal variants per gene
n_shared_causal_variants=9
n_tissue_specific_causal_variants=4
n_total_causal_variants = n_shared_causal_variants + (2*n_tissue_specific_causal_variants)
cis_snp_h2s = np.asarray([.0025, .005, .015])

# Set random seed
np.random.seed(simulation_iter)

# Open output file handle 
t1 = gzip.open(t1_causal_variant_gene_effect_size_file,'wt')
t1.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\teffect_size\n')

# Open output file handle 
t2 = gzip.open(t2_causal_variant_gene_effect_size_file,'wt')
t2.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\teffect_size\n')


# Loop through genes
f = open(gene_ld_summary_file)
head_count = 0
n_genes = 0
n_skipped_genes = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	gene_id = data[0]
	variant_file = data[4]
	cis_variants, cis_variant_a0s, cis_variant_a1s, cis_chroms, cis_poss = load_in_variants(variant_file)

	tot_variants_in_gene = len(cis_variants)
	if tot_variants_in_gene < n_total_causal_variants:
		n_skipped_genes = n_skipped_genes + 1
		continue

	n_genes = n_genes + 1
	t1_causal_variant_to_gene_effect_sizes, t2_causal_variant_to_gene_effect_sizes = simulate_cross_context_causal_effect_sizes(tot_variants_in_gene, n_shared_causal_variants, n_tissue_specific_causal_variants, cis_snp_h2s)
	
	for var_iter, variant_name in enumerate(cis_variants):
		t1.write(gene_id + '\t' + variant_name + '\t' + cis_chroms[var_iter] + '\t' + cis_poss[var_iter] + '\t' + cis_variant_a0s[var_iter] + '\t' + cis_variant_a1s[var_iter] + '\t' + str(t1_causal_variant_to_gene_effect_sizes[var_iter]) + '\n')
		t2.write(gene_id + '\t' + variant_name + '\t' + cis_chroms[var_iter] + '\t' + cis_poss[var_iter] + '\t' + cis_variant_a0s[var_iter] + '\t' + cis_variant_a1s[var_iter] + '\t' + str(t2_causal_variant_to_gene_effect_sizes[var_iter]) + '\n')

f.close()
t1.close()
t2.close()

sys.stderr.write('Simulated cross-context causal effects for ' + str(n_genes) + ' genes; skipped ' + str(n_skipped_genes) + ' genes with fewer than ' + str(n_total_causal_variants) + ' cis variants.\n')
