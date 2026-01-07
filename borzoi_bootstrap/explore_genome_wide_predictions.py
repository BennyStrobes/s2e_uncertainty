import numpy as np
import os
import sys
import pdb
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.integrate import quad
import scipy.stats
from scipy.optimize import brentq
from scipy.stats import norm
import matplotlib.cm as cm
from matplotlib.colors import Normalize



def plot_pvalue_histogram(borzoi_gaus_pvalues, output_file):
	plt.figure(figsize=(6,4))
	plt.hist(borzoi_gaus_pvalues, bins=30, edgecolor="black", alpha=0.7)

	plt.xlabel("P-value")
	plt.ylabel("Count")

	plt.tight_layout()
	plt.savefig(output_file, dpi=300)
	plt.close()
	return

def extract_borzoi_data(genome_wide_pred_summary_file):
	f = open(genome_wide_pred_summary_file)
	head_count = 0
	counter = 0
	effects = []
	pvalues = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		bs_vals = np.asarray(data[5].split(';')).astype(float)

		borzoi_effect = np.mean(bs_vals)

		np_p_plus = (1 + np.sum(bs_vals >= 0))/(1+len(bs_vals))
		np_p_minus = (1 + np.sum(bs_vals <=0))/(1+len(bs_vals))
		borzoi_pvalue = 2.0*np.min((np_p_plus, np_p_minus))


		effects.append(borzoi_effect)
		pvalues.append(borzoi_pvalue)
		counter = counter + 1
		if counter > 200000:
			break
	f.close()

	return np.asarray(effects), np.asarray(pvalues)





####################
# Command line args
####################
genome_wide_pred_summary_file = sys.argv[1]
output_root = sys.argv[2]


# Load in data
borzoi_effects, borzoi_pvalues = extract_borzoi_data(genome_wide_pred_summary_file)


# Plot histogram of p-values
output_file = output_root + 'pvalue_histogram.png'
plot_pvalue_histogram(borzoi_pvalues, output_file)