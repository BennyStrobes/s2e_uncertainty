import gzip
import math
import os
import sys

import numpy as np


def open_maybe_gzip(file_name, mode='rt'):
	if file_name.endswith('.gz'):
		return gzip.open(file_name, mode)
	return open(file_name, mode)


def create_variant_aliases(variant_id):
	aliases = [variant_id]
	if variant_id.endswith('_b38'):
		aliases.append(variant_id[:-4])
	else:
		aliases.append(variant_id + '_b38')
	return aliases


def normal_false_sign_rate(mean, variance):
	if variance < 0.0:
		variance = 0.0
	if variance == 0.0:
		if mean == 0.0:
			return 0.5
		return 0.0
	z_score = abs(mean)/math.sqrt(variance)
	return 0.5*math.erfc(z_score/math.sqrt(2.0))


def create_continuous_piecewise_linear_squared_basis(delta, knots):
	delta_sq = np.square(np.asarray(delta))
	basis = [np.ones(len(delta_sq)), delta_sq]
	for knot in knots:
		basis.append(np.maximum(delta_sq - knot, 0.0))
	return np.transpose(np.vstack(basis))


def load_ldscore_grid_squared_distribution(borzoi_based_prior_output_stem):
	summary_file = borzoi_based_prior_output_stem + '_ldscore_grid_squared_summary.txt'
	param_file = borzoi_based_prior_output_stem + '_ldscore_grid_squared_params.txt'

	grid_knots = None
	model_effect_scale = None
	with open(summary_file) as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if data[0] == 'grid_knots':
				grid_knots = np.asarray(data[1].split(',')).astype(float)
			elif data[0] == 'model_effect_scale':
				model_effect_scale = data[1]

	if grid_knots is None:
		raise ValueError('Could not find grid_knots in ' + summary_file)
	if model_effect_scale != 'allelic_grid_squared':
		raise ValueError('Expected model_effect_scale=allelic_grid_squared in ' + summary_file)

	a_prior = None
	a_prior_var = 0.0
	grid_coef_mapping = {}
	with open(param_file) as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			param_name = data[0]
			param_index = int(data[1])
			param_val = float(data[2])
			if param_name == 'a_prior':
				a_prior = param_val
			elif param_name == 'a_prior_var':
				a_prior_var = param_val
			elif param_name == 'grid_coef':
				grid_coef_mapping[param_index] = param_val

	if a_prior is None:
		raise ValueError('Could not find a_prior in ' + param_file)

	ordered_coef_indices = np.sort(np.asarray([*grid_coef_mapping]))
	expected_indices = np.arange(len(grid_knots) + 2)
	if np.array_equal(ordered_coef_indices, expected_indices) == False:
		raise ValueError('Grid coefficient indices do not match expected basis size in ' + param_file)

	return {
		'a_prior': a_prior,
		'a_prior_var': a_prior_var,
		'grid_knots': grid_knots,
		'grid_coefs': np.asarray([grid_coef_mapping[index] for index in ordered_coef_indices])
	}


def load_high_pip_finemapped_pairs(gtex_fm_file, target_tissue, fine_mapping_method, pip_threshold):
	fm_rows = []
	query_pairs = {}
	with open_maybe_gzip(gtex_fm_file) as f:
		header = f.readline().rstrip().split('\t')
		col = {name: ii for ii, name in enumerate(header)}
		required_cols = [
			'variant',
			'variant_hg38',
			'allele1',
			'allele2',
			'method',
			'tissue',
			'gene',
			'pip',
			'beta_posterior',
			'sd_posterior'
		]
		for required_col in required_cols:
			if required_col not in col:
				raise ValueError('Missing required column ' + required_col + ' in ' + gtex_fm_file)

		for line in f:
			line = line.rstrip()
			if line == '':
				continue
			data = line.split('\t')
			if len(data) != len(header):
				continue
			if data[col['tissue']] != target_tissue:
				continue
			if fine_mapping_method != 'all' and data[col['method']].lower() != fine_mapping_method.lower():
				continue
			try:
				pip = float(data[col['pip']])
			except ValueError:
				continue
			if pip <= pip_threshold:
				continue

			gene_id = data[col['gene']]
			gene_stable_id = gene_id.split('.')[0]
			variant_id = data[col['variant_hg38']]
			variant_aliases = create_variant_aliases(variant_id) + create_variant_aliases(data[col['variant']])

			row = {
				'chromosome': data[col['chromosome']] if 'chromosome' in col else 'NA',
				'start': data[col['start']] if 'start' in col else 'NA',
				'end': data[col['end']] if 'end' in col else 'NA',
				'variant': data[col['variant']],
				'variant_hg38': data[col['variant_hg38']],
				'allele1': data[col['allele1']],
				'allele2': data[col['allele2']],
				'method': data[col['method']],
				'tissue': data[col['tissue']],
				'gene': gene_id,
				'gene_stable_id': gene_stable_id,
				'pip': data[col['pip']],
				'beta_posterior': data[col['beta_posterior']],
				'sd_posterior': data[col['sd_posterior']],
				'variant_aliases': variant_aliases
			}
			fm_rows.append(row)

			if gene_stable_id not in query_pairs:
				query_pairs[gene_stable_id] = {}
			for variant_alias in variant_aliases:
				if variant_alias not in query_pairs[gene_stable_id]:
					query_pairs[gene_stable_id][variant_alias] = []
				query_pairs[gene_stable_id][variant_alias].append(len(fm_rows) - 1)

	return fm_rows, query_pairs


def annotate_rows_with_borzoi_effects(fm_rows, query_pairs, borzoi_effect_file):
	with open_maybe_gzip(borzoi_effect_file) as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if len(data) < 7:
				continue
			gene_id = data[0]
			variant_id = data[1]
			if gene_id not in query_pairs:
				continue
			if variant_id not in query_pairs[gene_id]:
				continue
			for row_index in query_pairs[gene_id][variant_id]:
				fm_rows[row_index]['borzoi_variant'] = variant_id
				fm_rows[row_index]['borzoi_a0'] = data[4]
				fm_rows[row_index]['borzoi_a1'] = data[5]
				fm_rows[row_index]['borzoi_effect_size'] = data[6]


def write_output(fm_rows, distribution_obj, output_file):
	output_dir = os.path.dirname(output_file)
	if output_dir != '':
		os.makedirs(output_dir, exist_ok=True)

	header = [
		'chromosome',
		'start',
		'end',
		'variant',
		'variant_hg38',
		'allele1',
		'allele2',
		'method',
		'tissue',
		'gene',
		'gene_stable_id',
		'pip',
		'fm_eqtl_effect_size',
		'fm_eqtl_effect_sd',
		'borzoi_variant',
		'borzoi_a0',
		'borzoi_a1',
		'borzoi_effect_size',
		'borzoi_prior_mean',
		'borzoi_prior_variance',
		'borzoi_prior_sd',
		'borzoi_prior_z',
		'borzoi_prior_false_sign_rate',
		'borzoi_prior_false_sign_rate_with_a_prior_uncertainty',
		'fm_eqtl_z',
		'fm_eqtl_false_sign_rate',
		'borzoi_fm_sign_agree'
	]

	with open(output_file, 'w') as t:
		t.write('\t'.join(header) + '\n')
		for row in fm_rows:
			borzoi_effect_string = row.get('borzoi_effect_size', 'nan')
			try:
				borzoi_effect = float(borzoi_effect_string)
				prior_mean = distribution_obj['a_prior']*borzoi_effect
				grid_basis = create_continuous_piecewise_linear_squared_basis(np.asarray([borzoi_effect]), distribution_obj['grid_knots'])
				prior_variance = np.dot(grid_basis, distribution_obj['grid_coefs'])[0]
				prior_variance = max(prior_variance, 0.0)
				fsr = normal_false_sign_rate(prior_mean, prior_variance)
				prior_sd = math.sqrt(prior_variance)
				prior_z = prior_mean/prior_sd if prior_sd > 0.0 else np.nan
				prior_variance_with_a_uncertainty = prior_variance + (np.square(borzoi_effect)*distribution_obj['a_prior_var'])
				fsr_with_a_uncertainty = normal_false_sign_rate(prior_mean, prior_variance_with_a_uncertainty)
			except ValueError:
				prior_mean = np.nan
				prior_variance = np.nan
				prior_sd = np.nan
				prior_z = np.nan
				fsr = np.nan
				fsr_with_a_uncertainty = np.nan

			try:
				fm_effect = float(row['beta_posterior'])
				fm_sd = float(row['sd_posterior'])
				fm_z = fm_effect/fm_sd if fm_sd > 0.0 else np.nan
				fm_fsr = normal_false_sign_rate(fm_effect, np.square(fm_sd))
				if np.isfinite(prior_mean) and np.isfinite(fm_effect) and prior_mean != 0.0 and fm_effect != 0.0:
					sign_agree = str(np.sign(prior_mean) == np.sign(fm_effect))
				else:
					sign_agree = 'nan'
			except ValueError:
				fm_z = np.nan
				fm_fsr = np.nan
				sign_agree = 'nan'

			out_data = [
				row['chromosome'],
				row['start'],
				row['end'],
				row['variant'],
				row['variant_hg38'],
				row['allele1'],
				row['allele2'],
				row['method'],
				row['tissue'],
				row['gene'],
				row['gene_stable_id'],
				row['pip'],
				row['beta_posterior'],
				row['sd_posterior'],
				row.get('borzoi_variant', 'NA'),
				row.get('borzoi_a0', 'NA'),
				row.get('borzoi_a1', 'NA'),
				borzoi_effect_string,
				str(prior_mean),
				str(prior_variance),
				str(prior_sd),
				str(prior_z),
				str(fsr),
				str(fsr_with_a_uncertainty),
				str(fm_z),
				str(fm_fsr),
				sign_agree
			]
			t.write('\t'.join(out_data) + '\n')


gtex_fm_file = sys.argv[1]
target_tissue = sys.argv[2]
fine_mapping_method = sys.argv[3]
pip_threshold = float(sys.argv[4])
borzoi_effect_file = sys.argv[5]
borzoi_based_prior_output_stem = sys.argv[6]
output_file = sys.argv[7]

distribution_obj = load_ldscore_grid_squared_distribution(borzoi_based_prior_output_stem)
fm_rows, query_pairs = load_high_pip_finemapped_pairs(gtex_fm_file, target_tissue, fine_mapping_method, pip_threshold)
annotate_rows_with_borzoi_effects(fm_rows, query_pairs, borzoi_effect_file)
write_output(fm_rows, distribution_obj, output_file)
