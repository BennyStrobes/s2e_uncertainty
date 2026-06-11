import gzip
import numpy as np
import sys


def open_possibly_gzipped_file(file_name, mode):
	if file_name.endswith('.gz'):
		return gzip.open(file_name, mode)
	return open(file_name, mode)


def load_effect_size_file(effect_size_file, effect_size_column_name):
	records = []
	effect_sizes = []
	with open_possibly_gzipped_file(effect_size_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				if len(data) < 7 or data[6] != effect_size_column_name:
					raise ValueError('Unexpected header in ' + effect_size_file + ': ' + '\t'.join(data))
				head_count = head_count + 1
				continue
			records.append(data[:6])
			effect_sizes.append(float(data[6]))
	return records, np.asarray(effect_sizes)


def check_paired_records(record_sets, record_set_names):
	reference_records = record_sets[0]
	for record_set, record_set_name in zip(record_sets[1:], record_set_names[1:]):
		if len(reference_records) != len(record_set):
			raise ValueError(record_set_names[0] + ' and ' + record_set_name + ' have different numbers of rows.')
		for ii, (reference_record, record) in enumerate(zip(reference_records, record_set)):
			if reference_record != record:
				raise ValueError(record_set_names[0] + ' and ' + record_set_name + ' are not aligned at row ' + str(ii) + '.')


def compute_regression_slope_and_intercept(x, y):
	x_mean = np.mean(x)
	y_mean = np.mean(y)
	x_var = np.var(x)
	if x_var == 0.0:
		return np.nan, np.nan
	slope = np.mean((x - x_mean)*(y - y_mean))/x_var
	intercept = y_mean - (slope*x_mean)
	return slope, intercept


def compute_pearson_correlation(x, y):
	if len(x) < 2:
		return np.nan
	if np.std(x) == 0.0 or np.std(y) == 0.0:
		return np.nan
	return np.corrcoef(x, y)[0, 1]


def summarize_effects(effect_type, causal_effects, borzoi_effects):
	pearson_correlation = compute_pearson_correlation(causal_effects, borzoi_effects)
	regression_slope, regression_intercept = compute_regression_slope_and_intercept(borzoi_effects, causal_effects)

	return [
		effect_type,
		str(len(causal_effects)),
		str(np.sum(causal_effects != 0.0)),
		str(pearson_correlation),
		str(regression_slope),
		str(regression_intercept),
		str(np.mean(borzoi_effects)),
		str(np.mean(causal_effects)),
		str(np.var(borzoi_effects)),
		str(np.var(causal_effects)),
		str(np.mean((borzoi_effects - np.mean(borzoi_effects))*(causal_effects - np.mean(causal_effects))))
	]


####################
# Command line args
####################
t1_est_borzoi_effect_size_file = sys.argv[1]
t2_est_borzoi_effect_size_file = sys.argv[2]
t1_causal_variant_gene_effect_size_file = sys.argv[3]
t2_causal_variant_gene_effect_size_file = sys.argv[4]
simulation_parameter_summary_file = sys.argv[5]


t1_borzoi_records, t1_est_borzoi = load_effect_size_file(t1_est_borzoi_effect_size_file, 'borzoi_effect_size')
t2_borzoi_records, t2_est_borzoi = load_effect_size_file(t2_est_borzoi_effect_size_file, 'borzoi_effect_size')
t1_causal_records, t1_causal = load_effect_size_file(t1_causal_variant_gene_effect_size_file, 'effect_size')
t2_causal_records, t2_causal = load_effect_size_file(t2_causal_variant_gene_effect_size_file, 'effect_size')

check_paired_records(
	[t1_borzoi_records, t2_borzoi_records, t1_causal_records, t2_causal_records],
	['t1_borzoi', 't2_borzoi', 't1_causal', 't2_causal']
)

delta_est_borzoi = t1_est_borzoi - t2_est_borzoi
delta_causal = t1_causal - t2_causal

with open(simulation_parameter_summary_file, 'w') as t:
	t.write('effect_type\tn_variant_gene_pairs\tn_nonzero_causal\tpearson_correlation\tregression_slope_causal_on_borzoi\tregression_intercept_causal_on_borzoi\tmean_borzoi\tmean_causal\tvar_borzoi\tvar_causal\tcov_borzoi_causal\n')
	t.write('\t'.join(summarize_effects('t1', t1_causal, t1_est_borzoi)) + '\n')
	t.write('\t'.join(summarize_effects('t2', t2_causal, t2_est_borzoi)) + '\n')
	t.write('\t'.join(summarize_effects('delta', delta_causal, delta_est_borzoi)) + '\n')
