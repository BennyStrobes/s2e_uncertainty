import numpy as np
import os
import sys
import pdb
import gzip
import matplotlib.pyplot as plt
from scipy.stats import norm


def load_in_sim_errors(sim_borzoi_error_file):
	f = gzip.open(sim_borzoi_error_file, 'rt')
	errors = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		errors.append(float(data[2]))
	f.close()
	return np.asarray(errors)


########################
# Command line args
########################
ashr_results_file = sys.argv[1]
sim_borzoi_error_file = sys.argv[2]
output_stem = sys.argv[3]

# Load in sim errors
sim_errors = load_in_sim_errors(sim_borzoi_error_file)

# Load in ashr-style results
tmp_data = np.loadtxt(ashr_results_file, dtype=str, delimiter='\t')
ashr_est_omega2 = tmp_data[1:, 0].astype(float)
ashr_est_pi = tmp_data[1:, 1].astype(float)

# Inputs
errors = sim_errors
omega2 = ashr_est_omega2
pis = ashr_est_pi


# Optional: renormalize in case of tiny numerical issues
pis = pis / np.sum(pis)

# Average variance implied by ASHR mixture
avg_var = np.sum(pis * omega2)
print(avg_var)

# Grid of x values
x_grid = np.linspace(np.percentile(errors, 0.1), np.percentile(errors, 99.9), 500)

############################
# Empirical CDF
############################
sorted_errors = np.sort(errors)
emp_cdf = np.searchsorted(sorted_errors, x_grid, side="right") / len(sorted_errors)

############################
# ASHR mixture CDF
############################
ashr_cdf = np.zeros_like(x_grid)

for pi, w2 in zip(pis, omega2):
	if w2 == 0:
		ashr_cdf += pi * (x_grid >= 0)
	else:
		ashr_cdf += pi * norm.cdf(x_grid, loc=0.0, scale=np.sqrt(w2))

############################
# Single Gaussian CDF with matched average variance
############################
if avg_var == 0:
	gaussian_cdf = (x_grid >= 0).astype(float)
else:
	gaussian_cdf = norm.cdf(x_grid, loc=0.0, scale=np.sqrt(avg_var))

############################
# Plot
############################
plt.figure(figsize=(6, 5))
plt.plot(x_grid, emp_cdf, label="Empirical CDF", linewidth=2)
plt.plot(x_grid, ashr_cdf, label="ASHR mixture CDF", linewidth=2)
plt.plot(x_grid, gaussian_cdf, label="Gaussian CDF (matched avg var)", linewidth=2)
plt.xlabel("Error")
plt.ylabel("CDF")
plt.legend()
plt.tight_layout()

plt.savefig(output_stem + '_cdf.pdf', dpi=300)
plt.close()

# Optional: print avg variance
print('ashr_average_variance\t' + str(avg_var))



