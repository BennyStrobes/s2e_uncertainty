import math
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_specific_example_data(input_file):
	f = open(input_file)
	head_count = 0
	positions = []
	gwas_z_scores = []
	borzoi_effects = []
	gene_id = None
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = {}
			for index, name in enumerate(data):
				header[name] = index
			continue
		if gene_id is None:
			gene_id = data[header['gene_id']]
		positions.append(float(data[header['position']]))
		gwas_z_scores.append(float(data[header['gwas_z']]))
		borzoi_effects.append(float(data[header['borzoi_effect_aligned_to_genotype']]))
	f.close()
	return gene_id, np.asarray(positions), np.asarray(gwas_z_scores), np.asarray(borzoi_effects)


def save_specific_example_locus_plot(input_file, output_file):
	gene_id, positions, gwas_z_scores, borzoi_effects = load_specific_example_data(input_file)
	valid_indices = np.isfinite(positions) & np.isfinite(gwas_z_scores) & np.isfinite(borzoi_effects)
	positions = positions[valid_indices]
	gwas_z_scores = gwas_z_scores[valid_indices]
	borzoi_effects = borzoi_effects[valid_indices]

	pvalues = np.asarray([math.erfc(abs(z_score)/math.sqrt(2.0)) for z_score in gwas_z_scores])
	pvalues = np.maximum(pvalues, np.nextafter(0.0, 1.0))
	neg_log10_pvalues = -1.0*np.log10(pvalues)

	fig, ax = plt.subplots(figsize=(5.2, 3.6))
	scatter = ax.scatter(positions, neg_log10_pvalues, c=np.abs(borzoi_effects), cmap='viridis', s=18, alpha=0.85, linewidths=0.0)
	ax.set_xlabel('SNP position')
	ax.set_ylabel('-log10(GWAS p-value)')
	ax.set_title(gene_id)
	cbar = fig.colorbar(scatter, ax=ax)
	cbar.set_label('|Borzoi effect size|')
	fig.tight_layout()
	fig.savefig(output_file)
	plt.close(fig)
	print(output_file)


########################
# Command line args
########################
input_file = sys.argv[1]
output_file = sys.argv[2]

save_specific_example_locus_plot(input_file, output_file)
