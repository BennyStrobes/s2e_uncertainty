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


def get_bootstrapped_estimates(string_arr):
	bs_ests = []

	for ele in string_arr:
		bs_ests.append(np.asarray(ele.split(';')).astype(float))

	return np.asarray(bs_ests)



def save_bootstrap_kde_plot(bootstraps, output_file):
	"""
	Generate and save a histogram + KDE plot from bootstrap samples.

	Parameters
	----------
	bootstraps : array-like
		The bootstrap estimates (1D array).
	output_file : str
		The output path, e.g. 'plot.png' or 'plot.pdf'.
	"""
	boot = np.asarray(bootstraps)
	
	# KDE
	kde = gaussian_kde(boot)
	
	# Evaluation grid
	x = np.linspace(
		boot.min() - 3 * boot.std(ddof=1),
		boot.max() + 3 * boot.std(ddof=1),
		400
	)
	pdf_kde = kde(x)
	
	# Plot
	plt.figure(figsize=(6, 4))
	plt.hist(boot, bins=15, density=True, alpha=0.5, edgecolor="k", label="Bootstrap histogram")
	plt.plot(x, pdf_kde, linewidth=2, label="KDE density")
	plt.xlabel("θ")
	plt.ylabel("Density")
	plt.title("Bootstrap Sampling Distribution (KDE)")
	plt.legend()
	
	# Save
	plt.savefig(output_file, dpi=300, bbox_inches="tight")
	plt.close()

	print(output_file)

	return output_file


def bootstrap_kde_pvalue(boot, theta_hat_obs, significance=0.05):
	"""
	Compute two-sided p-value using KDE estimated sampling distribution.
	
	Parameters
	----------
	boot : array-like
		Bootstrap estimates.
	theta_obs : float
		Observed estimate.
	theta0 : float
		Null value. Default is 0.
		
	Returns
	-------
	p_value : float
		Two-sided bootstrap KDE p-value.
	"""
	boot = np.asarray(boot)
	kde = gaussian_kde(boot)


	lo = boot.min() - 5*boot.std()
	hi = boot.max() + 5*boot.std()

	def cdf(x):
		return kde.integrate_box_1d(lo, x)

	# solve cdf(x) = 0.025 and 0.975
	ci_lower = brentq(lambda x: cdf(x) - significance/2.0, lo, hi)
	ci_upper = brentq(lambda x: cdf(x) - (1.0 - significance/2.0), lo, hi)

	if np.sign(ci_lower) == np.sign(ci_upper):
		return True
	else:
		return False

def get_same_sign_indices(borzoi_preds, thresh=0):
	nrows = borzoi_preds.shape[0]
	big_K = borzoi_preds.shape[1]
	arr = []
	for row_iter in range(nrows):
		kk = np.sum(np.sign(borzoi_preds[row_iter, :]) == 1.0)
		kk = np.min([kk, big_K-kk])
		if kk <= thresh:
			arr.append(True)
		else:
			arr.append(False)
	return np.asarray(arr)

def plot_correlation_dependence_on_number_of_bootstraps(bs_range, borzoi_bs, borzoi_mean, susie, output_file, n_subsamples=20):
	ncols = borzoi_bs.shape[1]
	num_bootstraps_arr = []
	avg_correlation = []

	for n_bs in bs_range:
		bs_corrs = []
		for subsample_iter in range(n_subsamples):
			subset_cols = np.random.choice(np.arange(ncols), size=n_bs, replace=False)
			same_sign_indices = get_same_sign_indices(borzoi_bs[:, subset_cols],thresh=0)

			corry = np.corrcoef(borzoi_mean[same_sign_indices], susie[same_sign_indices])[0,1]
			bs_corrs.append(corry)
		bs_corrs = np.asarray(bs_corrs)

		num_bootstraps_arr.append(n_bs)
		avg_correlation.append(np.mean(bs_corrs))
	num_bootstraps_arr = np.asarray(num_bootstraps_arr)
	avg_correlation = np.asarray(avg_correlation)

	plt.figure()                     # start fresh
	plt.scatter(num_bootstraps_arr, avg_correlation)
	plt.xlabel("Number of bootraps")
	plt.ylabel("Corr(borzoi effects, eqtl effects)")
	plt.tight_layout()               # avoids cut-off labels
	plt.savefig(output_file, dpi=300)
	plt.close() 

	return

def plot_num_hits_dependence_on_number_of_bootstraps(bs_range, borzoi_bs, borzoi_mean, susie, output_file, n_subsamples=20):
	ncols = borzoi_bs.shape[1]
	num_bootstraps_arr = []
	avg_correlation = []

	for n_bs in bs_range:
		bs_corrs = []
		for subsample_iter in range(n_subsamples):
			subset_cols = np.random.choice(np.arange(ncols), size=n_bs, replace=False)
			same_sign_indices = get_same_sign_indices(borzoi_bs[:, subset_cols],thresh=0)

			bs_corrs.append(np.sum(same_sign_indices))
		bs_corrs = np.asarray(bs_corrs)

		num_bootstraps_arr.append(n_bs)
		avg_correlation.append(np.mean(bs_corrs))
	num_bootstraps_arr = np.asarray(num_bootstraps_arr)
	avg_correlation = np.asarray(avg_correlation)

	plt.figure()                     # start fresh
	plt.scatter(num_bootstraps_arr, avg_correlation)
	plt.xlabel("Number of bootraps")
	plt.ylabel("Number of hits")
	plt.tight_layout()               # avoids cut-off labels
	plt.savefig(output_file, dpi=300)
	plt.close() 

	return 

def make_correlation_plot_colored_by_abs_borzoi(susie, borzoi_mean, z_score,output_file):
	susie = np.asarray(susie)
	borzoi_mean = np.asarray(borzoi_mean)
	abs_z = np.abs(z_score)

	# Convert correlations to 4 significant figures
	#corry_fmt = float(f"{corry:.4g}")
	#new_corry_fmt = float(f"{new_corry:.4g}")

	plt.figure(figsize=(6, 4))

	# -----------------------
	# Color map: white → blue
	# -----------------------
	cmap = cm.Blues  # already goes light → dark
	norm = Normalize(vmin=np.min(abs_z), vmax=np.max(abs_z))

	# -----------------------
	# Scatter plot
	# -----------------------
	plt.scatter(
	    susie,
	    borzoi_mean,
	    c=abs_z,
	    cmap=cmap,
	    norm=norm,
	    s=13,
	    alpha=0.8,
	    edgecolor="black",
	    linewidth=0.3
	)

	# -----------------------
	# Axes + title
	# -----------------------
	plt.xlabel("SuSiE Effect Size")
	plt.ylabel("Borzoi Mean Effect Size")
	#plt.title(f"Total corr: {corry_fmt} / Sig subset corr: {new_corry_fmt} / Total sig: {n_total_sig}")

	# -----------------------
	# Colorbar
	# -----------------------
	cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
	cbar.set_label('|z-score|')

	plt.tight_layout()
	plt.savefig(output_file, dpi=300)
	plt.close()

	return


def make_correlation_plot_colored_by_sig(susie, borzoi_mean, same_sign_indices,output_file):
	susie = np.asarray(susie)
	borzoi_mean = np.asarray(borzoi_mean)
	same_sign_indices = np.asarray(same_sign_indices)

	# Correlations
	corry = np.corrcoef(susie, borzoi_mean)[0,1]
	new_corry = np.corrcoef(susie[same_sign_indices], borzoi_mean[same_sign_indices])[0,1]

	n_total_sig = np.sum(same_sign_indices)

	print(np.sum(same_sign_indices))
	print(scipy.stats.pearsonr(susie[same_sign_indices],borzoi_mean[same_sign_indices]))


	# Convert correlations to 4 significant figures
	corry_fmt = float(f"{corry:.4g}")
	new_corry_fmt = float(f"{new_corry:.4g}")

	plt.figure(figsize=(6, 4))

	# --- Plot FALSE first (gray) ---
	plt.scatter(
		susie[~same_sign_indices],
		borzoi_mean[~same_sign_indices],
		c="gray",
		s=13,
		alpha=0.7,
		edgecolor="black",
		linewidth=0.3
	)

	# --- Plot TRUE second (blue) so they appear ON TOP ---
	plt.scatter(
		susie[same_sign_indices],
		borzoi_mean[same_sign_indices],
		c="blue",
		s=13,
		alpha=0.7,
		edgecolor="black",
		linewidth=0.3
	)

	plt.xlabel("SuSiE Effect Size")
	plt.ylabel("Borzoi Mean Effect Size")
	plt.title(f"Total corr: {corry_fmt} / Sig subset corr: {new_corry_fmt} / Total sig: {n_total_sig}")

	# Legend
	from matplotlib.lines import Line2D
	legend_elements = [
		Line2D([0], [0], marker='o', color='w', label='Borzoi sig.',
			   markerfacecolor='blue', markersize=8),
		Line2D([0], [0], marker='o', color='w', label='Background',
			   markerfacecolor='gray', markersize=8)
	]
	plt.legend(handles=legend_elements, loc="upper left")

	plt.tight_layout()
	plt.savefig(output_file, dpi=300)
	plt.close()

	return

def plot_pvalue_histogram(borzoi_gaus_pvalues, output_file):
	plt.figure(figsize=(6,4))
	plt.hist(borzoi_gaus_pvalues, bins=30, edgecolor="black", alpha=0.7)

	plt.xlabel("P-value")
	plt.ylabel("Count")

	plt.tight_layout()
	plt.savefig(output_file, dpi=300)
	plt.close()
	return

def get_borzoi_bs_pvalues(borzoi_bs, theta0=0.0):
	n_rows = borzoi_bs.shape[0]

	ps = []

	for row_iter in range(n_rows):
		theta_boot = borzoi_bs[row_iter, :]

		theta_boot_null = theta_boot - theta_boot.mean() + theta0


		B = len(theta_boot_null)

		p_boot = (np.sum(np.abs(theta_boot_null - theta0) >= np.abs(theta_boot.mean() - theta0)) + 1) / (B + 1)
		ps.append(p_boot)
	ps = np.asarray(ps)

	return ps

###################
# Command line args
####################
input_summary_file = sys.argv[1]
output_dir = sys.argv[2]


# Load in data
raw_data = np.loadtxt(input_summary_file, dtype=str, delimiter='\t')
susie = raw_data[1:,-5].astype(float)
borzoi_mean = raw_data[1:,-3].astype(float)
borzoi_se = raw_data[1:, -2].astype(float)
borzoi_bs = get_bootstrapped_estimates(raw_data[1:,-1])
borzoi_zeds = borzoi_mean/borzoi_se


# assume z is your numpy array of z-scores
borzoi_gaus_pvalues = 2 * (1 - norm.cdf(np.abs(borzoi_zeds)))

borzoi_bs_pvalues = get_borzoi_bs_pvalues(borzoi_bs)


thresh = [1.5, 2.0, 2.5]
for thresh in thresh:
	sig_indices = np.abs(borzoi_zeds) > thresh
	output_file = output_dir + 'borzoi_vs_qtl_scatter_colored_by_borzoi_sig_z_' + str(thresh) + '.png'
	make_correlation_plot_colored_by_sig(susie, borzoi_mean, sig_indices,output_file)

output_file = output_dir + 'borzoi_vs_qtl_scatter_colored_by_borzoi_abs_z.png'
make_correlation_plot_colored_by_abs_borzoi(susie, borzoi_mean, borzoi_zeds,output_file)	


same_sign_indices = (get_same_sign_indices(borzoi_bs,thresh=0)) #& (np.min(np.abs(borzoi_bs),axis=1) > .01)
output_file = output_dir + 'borzoi_vs_qtl_scatter_colored_by_borzoi_sig.png'
make_correlation_plot_colored_by_sig(susie, borzoi_mean, same_sign_indices,output_file)


# Plot correlation dependence on number of bootstraps
ncols = borzoi_bs.shape[1]
bs_range = np.arange(5, ncols+1)
output_file = output_dir + 'num_bootstraps_vs_eqtl_corr_scatter.png'
#plot_correlation_dependence_on_number_of_bootstraps(bs_range, borzoi_bs, borzoi_mean, susie, output_file)

# Plot num_hits dependence on number of bootstraps
ncols = borzoi_bs.shape[1]
bs_range = np.arange(5, ncols+1)
output_file = output_dir + 'num_bootstraps_vs_num_hits_scatter.png'
#plot_num_hits_dependence_on_number_of_bootstraps(bs_range, borzoi_bs, borzoi_mean, susie, output_file)


# Plot histogram of p-values
output_file = output_dir + 'pvalue_histogram.png'
#plot_pvalue_histogram(borzoi_gaus_pvalues, output_file)

# Plot histogram of p-values
output_file = output_dir + 'pvalue_bs_histogram.png'
#plot_pvalue_histogram(borzoi_bs_pvalues, output_file)



#save_bootstrap_kde_plot(borzoi_bs[570,:], output_dir + 'column570.png')

'''

for siggy in [.1, .05, .025, .01]:
	print(siggy)
	sigs = []
	for itera in range(borzoi_bs.shape[0]):
		pp = bootstrap_kde_pvalue(borzoi_bs[itera,:], np.mean(borzoi_bs[itera,:]), significance=siggy)
		sigs.append(pp)
	sigs = np.asarray(sigs)

	sig_indices = np.copy(sigs)
	print(siggy)
	print(np.sum(sig_indices))
	print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))
'''



'''
subset_cols = np.random.choice(np.arange(ncols), size=10, replace=False)
same_sign_indices = get_same_sign_indices(borzoi_bs[:, subset_cols],thresh=0)

print(np.corrcoef(borzoi_mean[same_sign_indices],susie[same_sign_indices]))

pdb.set_trace()
'''
'''
save_bootstrap_kde_plot(borzoi_bs[0,:], output_dir + 'column0.png')


for siggy in [.1, .05, .025, .01]:
	print(siggy)
	sigs = []
	for itera in range(borzoi_bs.shape[0]):
		pp = bootstrap_kde_pvalue(borzoi_bs[itera,:], np.mean(borzoi_bs[itera,:]), significance=siggy)
		sigs.append(pp)
	sigs = np.asarray(sigs)

	sig_indices = np.copy(sigs)
	print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))
pdb.set_trace()


save_bootstrap_kde_plot(borzoi_bs[1567,:], output_dir + 'column1567.png')
save_bootstrap_kde_plot(borzoi_bs[7906,:], output_dir + 'column7906.png')
save_bootstrap_kde_plot(borzoi_bs[5827,:], output_dir + 'column5827.png')
'''


