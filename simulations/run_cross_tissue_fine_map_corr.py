import gzip
import numpy as np
import sys


def load_effect_file(effect_file, expected_effect_column):
	mapping = {}
	with gzip.open(effect_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				if len(data) < 7 or data[6] != expected_effect_column:
					raise ValueError('Unexpected header in ' + effect_file + ': ' + '\t'.join(data))
				head_count = head_count + 1
				continue
			key = (data[0], data[1])
			if key in mapping:
				raise ValueError('Duplicate gene-variant pair in ' + effect_file + ': ' + data[0] + '\t' + data[1])
			mapping[key] = float(data[6])
	return mapping


def load_susie_file(susie_fine_mapping_file):
	mapping = {}
	with gzip.open(susie_fine_mapping_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			key = (data[0], data[1])
			if key in mapping:
				raise ValueError('Duplicate gene-variant pair in ' + susie_fine_mapping_file + ': ' + data[0] + '\t' + data[1])
			mapping[key] = {
				'pip': float(data[2]),
				'posterior_mean': float(data[3])
			}
	return mapping


def compute_correlation(x, y):
	valid_indices = np.isfinite(x) & np.isfinite(y)
	x = x[valid_indices]
	y = y[valid_indices]
	if len(x) < 2:
		return np.nan
	if np.std(x) == 0.0 or np.std(y) == 0.0:
		return np.nan
	return np.corrcoef(x, y)[0, 1]


def compute_regression_slope(x, y):
	valid_indices = np.isfinite(x) & np.isfinite(y)
	x = x[valid_indices]
	y = y[valid_indices]
	if len(x) < 2:
		return np.nan
	x_mean = np.mean(x)
	y_mean = np.mean(y)
	x_var = np.var(x)
	if x_var == 0.0:
		return np.nan
	return np.mean((x - x_mean)*(y - y_mean))/x_var


def compute_stats(borzoi_effects, eqtl_effects):
	return {
		'correlation': compute_correlation(borzoi_effects, eqtl_effects),
		'calibration_slope': compute_regression_slope(borzoi_effects, eqtl_effects)
	}


def summarize_bootstrap_distribution(observed_value, bootstrap_distribution):
	if np.all(np.isnan(bootstrap_distribution)):
		return np.nan, np.nan, np.nan, np.nan, np.nan
	bootstrap_mean = np.nanmean(bootstrap_distribution)
	bootstrap_se = np.nanstd(bootstrap_distribution, ddof=1)
	if bootstrap_se == 0.0 or np.isnan(bootstrap_se):
		gaussian_z_score = np.nan
	else:
		gaussian_z_score = observed_value/bootstrap_se
	empirical_ci = np.nanquantile(bootstrap_distribution, [.025, .975])
	return bootstrap_mean, bootstrap_se, gaussian_z_score, empirical_ci[0], empirical_ci[1]


def get_susie_info(susie_mapping, key):
	if key not in susie_mapping:
		return 0.0, 0.0
	return susie_mapping[key]['pip'], susie_mapping[key]['posterior_mean']


##########################
# Command line args
##########################
t1_est_borzoi_effect_size_file = sys.argv[1]
t2_est_borzoi_effect_size_file = sys.argv[2]
t1_susie_fine_mapping_file = sys.argv[3]
t2_susie_fine_mapping_file = sys.argv[4]
t1_causal_variant_gene_effect_size_file = sys.argv[5]
t2_causal_variant_gene_effect_size_file = sys.argv[6]
fm_corr_output_stem = sys.argv[7]
pip_thresh = float(sys.argv[8]) if len(sys.argv) > 8 else 0.9
n_bootstraps = int(sys.argv[9]) if len(sys.argv) > 9 else 100


t1_borzoi = load_effect_file(t1_est_borzoi_effect_size_file, 'borzoi_effect_size')
t2_borzoi = load_effect_file(t2_est_borzoi_effect_size_file, 'borzoi_effect_size')
t1_causal = load_effect_file(t1_causal_variant_gene_effect_size_file, 'effect_size')
t2_causal = load_effect_file(t2_causal_variant_gene_effect_size_file, 'effect_size')
t1_susie = load_susie_file(t1_susie_fine_mapping_file)
t2_susie = load_susie_file(t2_susie_fine_mapping_file)


selected_keys = []
for key in sorted(set(t1_susie.keys()).union(set(t2_susie.keys()))):
	t1_pip, _ = get_susie_info(t1_susie, key)
	t2_pip, _ = get_susie_info(t2_susie, key)
	if t1_pip >= pip_thresh or t2_pip >= pip_thresh:
		if key in t1_borzoi and key in t2_borzoi and key in t1_causal and key in t2_causal:
			selected_keys.append(key)


delta_borzoi = []
delta_fine_mapped_eqtl = []
delta_causal = []
for key in selected_keys:
	_, t1_posterior_mean = get_susie_info(t1_susie, key)
	_, t2_posterior_mean = get_susie_info(t2_susie, key)
	delta_borzoi.append(t1_borzoi[key] - t2_borzoi[key])
	delta_fine_mapped_eqtl.append(t1_posterior_mean - t2_posterior_mean)
	delta_causal.append(t1_causal[key] - t2_causal[key])

delta_borzoi = np.asarray(delta_borzoi)
delta_fine_mapped_eqtl = np.asarray(delta_fine_mapped_eqtl)
delta_causal = np.asarray(delta_causal)


n_selected = len(selected_keys)
output_names = ['correlation', 'calibration_slope']
observed_stats = {
	'fine_mapped_posterior_mean': compute_stats(delta_borzoi, delta_fine_mapped_eqtl),
	'oracle_selected_causal': compute_stats(delta_borzoi, delta_causal)
}

bootstrap_stats = {}
for estimator_name in observed_stats:
	bootstrap_stats[estimator_name] = {}
	for output_name in output_names:
		bootstrap_stats[estimator_name][output_name] = np.full(n_bootstraps, np.nan)

if n_selected >= 2:
	for bs_iter in range(n_bootstraps):
		bs_indices = np.random.choice(np.arange(n_selected), size=n_selected, replace=True)
		bs_observed_stats = {
			'fine_mapped_posterior_mean': compute_stats(delta_borzoi[bs_indices], delta_fine_mapped_eqtl[bs_indices]),
			'oracle_selected_causal': compute_stats(delta_borzoi[bs_indices], delta_causal[bs_indices])
		}
		for estimator_name in observed_stats:
			for output_name in output_names:
				bootstrap_stats[estimator_name][output_name][bs_iter] = bs_observed_stats[estimator_name][output_name]


summary_stats_output_file = fm_corr_output_stem + '_bootstrap_stats.txt'
with open(summary_stats_output_file, 'w') as t:
	t.write('estimator_name\toutput_name\tn_selected_variant_gene_pairs\tpip_threshold\tmean\tbootstrapped_mean\tbootstrap_se\tgaussian_z_score\tempirical_ci_lower\tempirical_ci_upper\n')
	for estimator_name in observed_stats:
		for output_name in output_names:
			observed_value = observed_stats[estimator_name][output_name]
			bootstrap_mean, bootstrap_se, gaussian_z_score, ci_lower, ci_upper = summarize_bootstrap_distribution(observed_value, bootstrap_stats[estimator_name][output_name])
			t.write(estimator_name + '\t' + output_name + '\t' + str(n_selected) + '\t' + str(pip_thresh) + '\t' + str(observed_value) + '\t' + str(bootstrap_mean) + '\t' + str(bootstrap_se) + '\t' + str(gaussian_z_score) + '\t' + str(ci_lower) + '\t' + str(ci_upper) + '\n')


selected_pairs_output_file = fm_corr_output_stem + '_selected_variant_gene_pairs.txt.gz'
with gzip.open(selected_pairs_output_file, 'wt') as t:
	t.write('gene\tvariant\tt1_pip\tt2_pip\tdelta_borzoi\tdelta_fine_mapped_eqtl\tdelta_causal\n')
	for key, cur_delta_borzoi, cur_delta_fm_eqtl, cur_delta_causal in zip(selected_keys, delta_borzoi, delta_fine_mapped_eqtl, delta_causal):
		t1_pip, _ = get_susie_info(t1_susie, key)
		t2_pip, _ = get_susie_info(t2_susie, key)
		t.write(key[0] + '\t' + key[1] + '\t' + str(t1_pip) + '\t' + str(t2_pip) + '\t' + str(cur_delta_borzoi) + '\t' + str(cur_delta_fm_eqtl) + '\t' + str(cur_delta_causal) + '\n')

print(summary_stats_output_file)
print(selected_pairs_output_file)
