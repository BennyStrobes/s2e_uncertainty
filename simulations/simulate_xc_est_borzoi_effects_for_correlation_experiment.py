import numpy as np
import sys
import gzip


def load_causal_effect_file(causal_variant_gene_effect_size_file):
	records = []
	effect_sizes = []
	f = gzip.open(causal_variant_gene_effect_size_file, 'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		records.append(data[:6])
		effect_sizes.append(float(data[6]))
	f.close()
	return records, np.asarray(effect_sizes)


def check_paired_records(t1_records, t2_records):
	if len(t1_records) != len(t2_records):
		raise ValueError('Tissue 1 and tissue 2 causal effect files have different numbers of rows.')

	for ii, (t1_record, t2_record) in enumerate(zip(t1_records, t2_records)):
		if t1_record != t2_record:
			raise ValueError('Tissue 1 and tissue 2 causal effect files are not aligned at row ' + str(ii) + '.')


def simulate_orthogonal_noise(reference_values, noise_variance):
	if noise_variance == 0.0:
		return np.zeros(len(reference_values))

	reference_values = reference_values - np.mean(reference_values)
	noise = np.random.normal(loc=0.0, scale=1.0, size=len(reference_values))
	noise = noise - np.mean(noise)
	noise = noise - (reference_values*(np.dot(noise, reference_values)/np.dot(reference_values, reference_values)))
	noise = noise - np.mean(noise)
	noise = noise*np.sqrt(noise_variance/np.var(noise))
	return noise


def simulate_correlated_tissue_errors(delta_causal_effects, delta_error_variance, error_correlation):
	if delta_error_variance == 0.0:
		return np.zeros(len(delta_causal_effects)), np.zeros(len(delta_causal_effects))
	if error_correlation <= -1.0 or error_correlation >= 1.0:
		raise ValueError('borzoi_error_correlation must be greater than -1 and less than 1 when noise is nonzero.')

	delta_error = simulate_orthogonal_noise(delta_causal_effects, delta_error_variance)
	common_error_variance = delta_error_variance*(1.0 + error_correlation)/(4.0*(1.0 - error_correlation))
	common_error = simulate_orthogonal_noise(delta_error, common_error_variance)

	t1_errors = common_error + (delta_error/2.0)
	t2_errors = common_error - (delta_error/2.0)
	return t1_errors, t2_errors


def simulate_borzoi_effects(t1_causal_effects, t2_causal_effects, rho_delta, borzoi_error_correlation):
	if rho_delta < 0.0 or rho_delta > 1.0:
		raise ValueError('rho_delta must be between 0 and 1.')

	delta_causal_effects = t1_causal_effects - t2_causal_effects
	delta_causal_variance = np.var(delta_causal_effects)
	if delta_causal_variance <= 0.0:
		raise ValueError('Causal tissue-difference effects have zero variance.')

	if rho_delta == 0.0:
		signal_slope = 0.0
		delta_error_variance = delta_causal_variance
	else:
		signal_slope = 1.0
		delta_error_variance = delta_causal_variance*((1.0/np.square(rho_delta)) - 1.0)

	t1_errors, t2_errors = simulate_correlated_tissue_errors(delta_causal_effects, delta_error_variance, borzoi_error_correlation)

	t1_borzoi_effects = (signal_slope*t1_causal_effects) + t1_errors
	t2_borzoi_effects = (signal_slope*t2_causal_effects) + t2_errors

	return t1_borzoi_effects, t2_borzoi_effects


def write_borzoi_effect_file(records, borzoi_effects, output_file):
	t = gzip.open(output_file, 'wt')
	t.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\tborzoi_effect_size\n')
	for record, borzoi_effect in zip(records, borzoi_effects):
		t.write('\t'.join(record) + '\t' + str(borzoi_effect) + '\n')
	t.close()


#####################
# Command line args
#####################
t1_causal_variant_gene_effect_size_file = sys.argv[1]
t2_causal_variant_gene_effect_size_file = sys.argv[2]
t1_est_borzoi_effect_size_file = sys.argv[3]
t2_est_borzoi_effect_size_file = sys.argv[4]
simulation_iter = int(sys.argv[5])
rho_delta = float(sys.argv[6])

if len(sys.argv) > 7:
	borzoi_error_correlation = float(sys.argv[7])
else:
	borzoi_error_correlation = 0.3


np.random.seed(simulation_iter)


t1_records, t1_causal_effects = load_causal_effect_file(t1_causal_variant_gene_effect_size_file)
t2_records, t2_causal_effects = load_causal_effect_file(t2_causal_variant_gene_effect_size_file)
check_paired_records(t1_records, t2_records)

t1_borzoi_effects, t2_borzoi_effects = simulate_borzoi_effects(t1_causal_effects, t2_causal_effects, rho_delta, borzoi_error_correlation)

write_borzoi_effect_file(t1_records, t1_borzoi_effects, t1_est_borzoi_effect_size_file)
write_borzoi_effect_file(t2_records, t2_borzoi_effects, t2_est_borzoi_effect_size_file)

realized_rho_delta = np.corrcoef(t1_causal_effects - t2_causal_effects, t1_borzoi_effects - t2_borzoi_effects)[0, 1]
sys.stderr.write('Simulated cross-context Borzoi effects with target rho_delta=' + str(rho_delta) + ', realized rho_delta=' + str(realized_rho_delta) + ', and borzoi_error_correlation=' + str(borzoi_error_correlation) + '.\n')
