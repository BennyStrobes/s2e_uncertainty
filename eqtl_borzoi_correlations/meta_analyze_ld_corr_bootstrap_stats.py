#!/usr/bin/env python3

import argparse
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np


def normal_two_sided_pvalue(z_score):
	return math.erfc(abs(z_score)/math.sqrt(2.0))


def benjamini_hochberg(pvalues):
	n_tests = len(pvalues)
	if n_tests == 0:
		return []
	order = np.argsort(pvalues)
	adjusted = np.ones(n_tests)
	running_min = 1.0
	for reversed_rank, test_index in enumerate(order[::-1], start=1):
		rank = n_tests - reversed_rank + 1
		adjusted_pvalue = pvalues[test_index]*n_tests/float(rank)
		running_min = min(running_min, adjusted_pvalue)
		adjusted[test_index] = min(1.0, running_min)
	return adjusted.tolist()


def bonferroni(pvalues):
	n_tests = len(pvalues)
	return [min(1.0, pvalue*n_tests) for pvalue in pvalues]


def annotation_sort_key(annotation_name):
	match = re.search(r'(\d+)$', annotation_name)
	if match is None:
		return (annotation_name, -1)
	return (annotation_name[:match.start()], int(match.group(1)))


def extract_source_label(raw_source):
	base = os.path.basename(raw_source)
	match = re.search(r'ld_corr_results_(.+?)_GTEX-', base)
	if match is not None:
		return match.group(1)
	return base


def parse_stats_line(line, fallback_source):
	line = line.rstrip()
	if line == '' or line.startswith('#'):
		return None
	data = line.split('\t')
	if len(data) < 5:
		raise ValueError('Malformed line with fewer than 5 columns: ' + line)
	if data[0] == 'annotation_name' and data[1] == 'output_name':
		return None

	first_entry = data[0]
	if ':' in first_entry and not first_entry.startswith('magnitude_bin'):
		raw_source, annotation_name = first_entry.rsplit(':', 1)
	else:
		raw_source = fallback_source
		annotation_name = first_entry

	return {
		'source': extract_source_label(raw_source),
		'annotation_name': annotation_name,
		'output_name': data[1],
		'mean': float(data[2]),
		'bootstrap_mean': float(data[3]),
		'bootstrap_se': float(data[4]),
	}


def load_records(input_files, target_output_name):
	records = []
	for input_file in input_files:
		with open(input_file) as f:
			for line in f:
				record = parse_stats_line(line, input_file)
				if record is None:
					continue
				if record['output_name'] != target_output_name:
					continue
				if not np.isfinite(record['mean']) or not np.isfinite(record['bootstrap_se']):
					continue
				if record['bootstrap_se'] <= 0.0:
					continue
				records.append(record)
	return records


def meta_analyze_records(records):
	annotation_to_records = defaultdict(list)
	for record in records:
		annotation_to_records[record['annotation_name']].append(record)

	meta_results = []
	for annotation_name in sorted(annotation_to_records.keys(), key=annotation_sort_key):
		annotation_records = annotation_to_records[annotation_name]
		means = np.asarray([record['mean'] for record in annotation_records])
		standard_errors = np.asarray([record['bootstrap_se'] for record in annotation_records])
		weights = 1.0/np.square(standard_errors)
		meta_mean = np.sum(weights*means)/np.sum(weights)
		meta_se = np.sqrt(1.0/np.sum(weights))
		meta_z = meta_mean/meta_se
		meta_pvalue = normal_two_sided_pvalue(meta_z)
		meta_results.append({
			'annotation_name': annotation_name,
			'n_sources': len(annotation_records),
			'meta_mean': meta_mean,
			'meta_se': meta_se,
			'meta_z': meta_z,
			'meta_pvalue': meta_pvalue,
			'ci_lower': meta_mean - 1.96*meta_se,
			'ci_upper': meta_mean + 1.96*meta_se,
			'sources': ','.join(sorted(record['source'] for record in annotation_records)),
		})
	return meta_results


def pairwise_meta_differences(meta_results):
	pairwise_results = []
	for annotation_iter in range(len(meta_results)):
		for annotation_iter2 in range(annotation_iter + 1, len(meta_results)):
			result_1 = meta_results[annotation_iter]
			result_2 = meta_results[annotation_iter2]
			diff = result_1['meta_mean'] - result_2['meta_mean']
			diff_se = np.sqrt(np.square(result_1['meta_se']) + np.square(result_2['meta_se']))
			diff_z = diff/diff_se
			raw_pvalue = normal_two_sided_pvalue(diff_z)
			pairwise_results.append({
				'annotation_1': result_1['annotation_name'],
				'annotation_2': result_2['annotation_name'],
				'mean_diff': diff,
				'diff_se': diff_se,
				'diff_z': diff_z,
				'raw_pvalue': raw_pvalue,
			})

	raw_pvalues = [result['raw_pvalue'] for result in pairwise_results]
	bonferroni_pvalues = bonferroni(raw_pvalues)
	bh_pvalues = benjamini_hochberg(raw_pvalues)
	for result_iter, result in enumerate(pairwise_results):
		result['bonferroni_pvalue'] = bonferroni_pvalues[result_iter]
		result['bh_fdr_pvalue'] = bh_pvalues[result_iter]
		result['bonferroni_significant'] = bonferroni_pvalues[result_iter] < 0.05
		result['bh_fdr_significant'] = bh_pvalues[result_iter] < 0.05
	return pairwise_results


def print_or_write_results(meta_results, pairwise_results, output_stem):
	meta_header = 'annotation_name\tn_sources\tmeta_mean\tmeta_se\tmeta_z\tmeta_pvalue\tci_lower_95\tci_upper_95\tsources'
	pairwise_header = 'annotation_1\tannotation_2\tmean_diff\tdiff_se\tdiff_z\traw_pvalue\tbonferroni_pvalue\tbh_fdr_pvalue\tbonferroni_significant\tbh_fdr_significant'
	pairwise_note = '# Pairwise z-tests below treat the meta-analyzed bin estimates as independent. Because the same tissues contribute to multiple bins, these pairwise p-values are approximate.'

	if output_stem is None:
		print(meta_header)
		for result in meta_results:
			print(
				result['annotation_name'] + '\t' +
				str(result['n_sources']) + '\t' +
				str(result['meta_mean']) + '\t' +
				str(result['meta_se']) + '\t' +
				str(result['meta_z']) + '\t' +
				str(result['meta_pvalue']) + '\t' +
				str(result['ci_lower']) + '\t' +
				str(result['ci_upper']) + '\t' +
				result['sources']
			)
		print('')
		print(pairwise_note)
		print(pairwise_header)
		for result in pairwise_results:
			print(
				result['annotation_1'] + '\t' +
				result['annotation_2'] + '\t' +
				str(result['mean_diff']) + '\t' +
				str(result['diff_se']) + '\t' +
				str(result['diff_z']) + '\t' +
				str(result['raw_pvalue']) + '\t' +
				str(result['bonferroni_pvalue']) + '\t' +
				str(result['bh_fdr_pvalue']) + '\t' +
				str(result['bonferroni_significant']) + '\t' +
				str(result['bh_fdr_significant'])
			)
		return

	meta_output_file = output_stem + '_meta_analysis.txt'
	with open(meta_output_file, 'w') as t:
		t.write(meta_header + '\n')
		for result in meta_results:
			t.write(
				result['annotation_name'] + '\t' +
				str(result['n_sources']) + '\t' +
				str(result['meta_mean']) + '\t' +
				str(result['meta_se']) + '\t' +
				str(result['meta_z']) + '\t' +
				str(result['meta_pvalue']) + '\t' +
				str(result['ci_lower']) + '\t' +
				str(result['ci_upper']) + '\t' +
				result['sources'] + '\n'
			)

	pairwise_output_file = output_stem + '_pairwise_differences.txt'
	with open(pairwise_output_file, 'w') as t:
		t.write(pairwise_note + '\n')
		t.write(pairwise_header + '\n')
		for result in pairwise_results:
			t.write(
				result['annotation_1'] + '\t' +
				result['annotation_2'] + '\t' +
				str(result['mean_diff']) + '\t' +
				str(result['diff_se']) + '\t' +
				str(result['diff_z']) + '\t' +
				str(result['raw_pvalue']) + '\t' +
				str(result['bonferroni_pvalue']) + '\t' +
				str(result['bh_fdr_pvalue']) + '\t' +
				str(result['bonferroni_significant']) + '\t' +
				str(result['bh_fdr_significant']) + '\n'
			)

	print(meta_output_file)
	print(pairwise_output_file)


def main():
	parser = argparse.ArgumentParser(description='Fixed-effect inverse-variance meta-analysis of ld_corr bootstrap summary statistics across tissues.')
	parser.add_argument('input_files', nargs='+', help='One or more *_bootstrap_stats.txt files, or a concatenated grep-style file containing rows from many such files.')
	parser.add_argument('--output-name', default='calibration_slope', help='Output name to meta-analyze. Default: calibration_slope')
	parser.add_argument('--output-stem', default=None, help='If provided, write results to <output_stem>_meta_analysis.txt and <output_stem>_pairwise_differences.txt')
	args = parser.parse_args()

	records = load_records(args.input_files, args.output_name)
	if len(records) == 0:
		raise ValueError('No valid records found for output_name=' + args.output_name)

	meta_results = meta_analyze_records(records)
	pairwise_results = pairwise_meta_differences(meta_results)
	print_or_write_results(meta_results, pairwise_results, args.output_stem)


if __name__ == '__main__':
	main()
