import numpy as np
import sys
import pdb
import gzip
import pyarrow.parquet as pq


def write_annotation_header(data, t, anno_names):
	t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
	t.write('\t' + '\t'.join(anno_names) + '\n')


def create_one_hot_annotation_from_vg_pair_mapping(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, vg_pair_info_index, bins, bin_prefix):
	f = gzip.open(borzoi_effect_file,'rt')
	t = gzip.open(borzoi_annotation_file, 'wt')
	bin_counts = np.zeros(len(bins)-1, dtype=int)
	anno_names = []
	for bin_iter in range(len(bins)-1):
		anno_names.append(bin_prefix + str(bin_iter))

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			write_annotation_header(data, t, anno_names)
			continue

		vg_name = data[1] + ':' + data[0]
		if vg_name not in vg_pair_info:
			bin_membership = np.zeros(len(bins)-1)
		else:
			value = vg_pair_info[vg_name][vg_pair_info_index]
			bin_membership = np.zeros(len(bins)-1)
			for bin_iter in range(len(bins)-1):
				if value >= bins[bin_iter] and value < bins[bin_iter+1]:
					bin_membership[bin_iter] = 1.0
					bin_counts[bin_iter] = bin_counts[bin_iter] + 1
					break


		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
		t.write('\t' + '\t'.join(bin_membership.astype(str)) + '\n')

	f.close()
	t.close()

	for bin_iter in range(len(bins)-1):
		print(bin_prefix + str(bin_iter) + ' [' + str(bins[bin_iter]) + ', ' + str(bins[bin_iter+1]) + '): ' + str(bin_counts[bin_iter]))
	return


def create_annotation_using_borzoi_magnitude_bins(borzoi_effect_file, borzoi_annotation_file, bins=[0.0, 0.01, 0.075, 0.2, 10.0]):
	f = gzip.open(borzoi_effect_file,'rt')
	t = gzip.open(borzoi_annotation_file, 'wt')
	bin_counts = np.zeros(len(bins)-1, dtype=int)

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
			for bin_iter in range(len(bins)-1):
				t.write('\tmagnitude_bin' + str(bin_iter))
			t.write('\n')
			continue

		borzoi_effect_size = np.abs(float(data[6]))
		bin_membership = np.zeros(len(bins)-1)
		for bin_iter in range(len(bins)-1):
			if borzoi_effect_size >= bins[bin_iter] and borzoi_effect_size < bins[bin_iter+1]:
				bin_membership[bin_iter] = 1.0
				bin_counts[bin_iter] = bin_counts[bin_iter] + 1
				break

		if np.sum(bin_membership) != 1.0:
			print('assumption eroror')
			pdb.set_trace()

		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
		t.write('\t' + '\t'.join(bin_membership.astype(str)) + '\n')

	f.close()
	t.close()

	for bin_iter in range(len(bins)-1):
		print('magnitude_bin' + str(bin_iter) + ' [' + str(bins[bin_iter]) + ', ' + str(bins[bin_iter+1]) + '): ' + str(bin_counts[bin_iter]))
	return


def create_annotation_using_borzoi_magnitude_bins_x_maf_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, magnitude_bins=[0.0, 0.01, 0.075, 0.2, 10.0], maf_bins=[0.0, 0.05, 0.1, 0.2, 0.3, 0.500001]):
	f = gzip.open(borzoi_effect_file,'rt')
	t = gzip.open(borzoi_annotation_file, 'wt')
	n_magnitude_bins = len(magnitude_bins) - 1
	n_maf_bins = len(maf_bins) - 1
	bin_counts = np.zeros((n_magnitude_bins, n_maf_bins), dtype=int)
	anno_names = []
	for magnitude_bin_iter in range(n_magnitude_bins):
		for maf_bin_iter in range(n_maf_bins):
			anno_names.append('magnitude_bin' + str(magnitude_bin_iter) + 'Xmaf_bin' + str(maf_bin_iter))

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			write_annotation_header(data, t, anno_names)
			continue

		bin_membership = np.zeros(n_magnitude_bins*n_maf_bins)
		vg_name = data[1] + ':' + data[0]
		if vg_name in vg_pair_info:
			borzoi_effect_size = np.abs(float(data[6]))
			maf = vg_pair_info[vg_name][0]

			magnitude_bin_index = -1
			for magnitude_bin_iter in range(n_magnitude_bins):
				if borzoi_effect_size >= magnitude_bins[magnitude_bin_iter] and borzoi_effect_size < magnitude_bins[magnitude_bin_iter+1]:
					magnitude_bin_index = magnitude_bin_iter
					break

			maf_bin_index = -1
			for maf_bin_iter in range(n_maf_bins):
				if maf >= maf_bins[maf_bin_iter] and maf < maf_bins[maf_bin_iter+1]:
					maf_bin_index = maf_bin_iter
					break

			if magnitude_bin_index != -1 and maf_bin_index != -1:
				interaction_bin_index = (magnitude_bin_index*n_maf_bins) + maf_bin_index
				bin_membership[interaction_bin_index] = 1.0
				bin_counts[magnitude_bin_index, maf_bin_index] = bin_counts[magnitude_bin_index, maf_bin_index] + 1

		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
		t.write('\t' + '\t'.join(bin_membership.astype(str)) + '\n')

	f.close()
	t.close()

	for magnitude_bin_iter in range(n_magnitude_bins):
		for maf_bin_iter in range(n_maf_bins):
			print('magnitude_bin' + str(magnitude_bin_iter) + 'Xmaf_bin' + str(maf_bin_iter) + ': ' + str(bin_counts[magnitude_bin_iter, maf_bin_iter]))
	return


def create_annotation_using_borzoi_magnitude_bins_x_dist_to_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, magnitude_bins=[0.0, 0.01, 0.075, 0.2, 10.0], dist_bins=[0.0, 1000.0, 5000.0, 20000.0, 50000.0, 100001.0]):
	f = gzip.open(borzoi_effect_file,'rt')
	t = gzip.open(borzoi_annotation_file, 'wt')
	n_magnitude_bins = len(magnitude_bins) - 1
	n_dist_bins = len(dist_bins) - 1
	bin_counts = np.zeros((n_magnitude_bins, n_dist_bins), dtype=int)
	anno_names = []
	for magnitude_bin_iter in range(n_magnitude_bins):
		for dist_bin_iter in range(n_dist_bins):
			anno_names.append('magnitude_bin' + str(magnitude_bin_iter) + 'Xdist_to_tss_bin' + str(dist_bin_iter))

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			write_annotation_header(data, t, anno_names)
			continue

		bin_membership = np.zeros(n_magnitude_bins*n_dist_bins)
		vg_name = data[1] + ':' + data[0]
		if vg_name in vg_pair_info:
			borzoi_effect_size = np.abs(float(data[6]))
			dist_to_tss = vg_pair_info[vg_name][1]

			magnitude_bin_index = -1
			for magnitude_bin_iter in range(n_magnitude_bins):
				if borzoi_effect_size >= magnitude_bins[magnitude_bin_iter] and borzoi_effect_size < magnitude_bins[magnitude_bin_iter+1]:
					magnitude_bin_index = magnitude_bin_iter
					break

			dist_bin_index = -1
			for dist_bin_iter in range(n_dist_bins):
				if dist_to_tss >= dist_bins[dist_bin_iter] and dist_to_tss < dist_bins[dist_bin_iter+1]:
					dist_bin_index = dist_bin_iter
					break

			if magnitude_bin_index != -1 and dist_bin_index != -1:
				interaction_bin_index = (magnitude_bin_index*n_dist_bins) + dist_bin_index
				bin_membership[interaction_bin_index] = 1.0
				bin_counts[magnitude_bin_index, dist_bin_index] = bin_counts[magnitude_bin_index, dist_bin_index] + 1

		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
		t.write('\t' + '\t'.join(bin_membership.astype(str)) + '\n')

	f.close()
	t.close()

	for magnitude_bin_iter in range(n_magnitude_bins):
		for dist_bin_iter in range(n_dist_bins):
			print('magnitude_bin' + str(magnitude_bin_iter) + 'Xdist_to_tss_bin' + str(dist_bin_iter) + ': ' + str(bin_counts[magnitude_bin_iter, dist_bin_iter]))
	return


def create_annotation_using_distance_to_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, bins=[0.0, 5000.0, 50000.0, 100001.0]):
	return create_one_hot_annotation_from_vg_pair_mapping(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, 1, bins, 'dist_to_tss_bin')

def create_annotation_using_distance_to_strand_based_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, bins=[-100001.0, -50000.0, -5000.0, 0.0, 5000.0, 50000.0, 100001.0]):
	return create_one_hot_annotation_from_vg_pair_mapping(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, 2, bins, 'strand_based_dist_to_tss_bin')


def create_annotation_using_maf_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, bins=[0.0, 0.05, 0.1, 0.2, 0.3, 0.500001]):
	return create_one_hot_annotation_from_vg_pair_mapping(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, 0, bins, 'maf_bin')


def extract_info_on_each_variant_gene_pair(eqtl_sumstats_file):
	dicti = {}

	f = gzip.open(eqtl_sumstats_file,'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		gene_id = data[0]
		var_id = data[1]
		maf = float(data[9])
		if maf > 1.0 or maf < 0.0:
			print('assumption errorro')
			pdb.set_trace()
		if maf > 0.5:
			maf = 1.0 - maf
		directional_distance = float(data[10])
		abs_tss_dist = np.abs(float(data[10]))
		vg_name = var_id + ':' + gene_id
		dicti[vg_name] = (maf, abs_tss_dist, directional_distance)
	f.close()
	return dicti


#####################
# Command line args
#####################
borzoi_effect_file = sys.argv[1]
anno_method = sys.argv[2]
borzoi_annotation_file = sys.argv[3]
eqtl_sumstats_file = sys.argv[4]
tissue_name = sys.argv[5]

# Create mapping from variant-gene pairs to auxiliary info
vg_pair_info = extract_info_on_each_variant_gene_pair(eqtl_sumstats_file)


if anno_method == 'borzoi_magnitude_bins':
	create_annotation_using_borzoi_magnitude_bins(borzoi_effect_file, borzoi_annotation_file, bins=[0.0, 0.001, 0.01, 0.075, 0.2, 10.0])
elif anno_method == 'dist_to_tss_bins':
	create_annotation_using_distance_to_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, bins=[0.0, 5000.0, 50000.0, 100001.0])
elif anno_method == 'strand_dist_to_tss_bins':
	create_annotation_using_distance_to_strand_based_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info,bins=[-100001.0, -50000.0, -5000.0, 0.0, 5000.0, 50000.0, 100001.0])
elif anno_method == 'af_bins':
	create_annotation_using_maf_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, bins=[0.0, 0.05, 0.1, 0.2, 0.3, 0.500001])
elif anno_method == 'borzoi_magnitude_binsXaf_bins':
	create_annotation_using_borzoi_magnitude_bins_x_maf_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, magnitude_bins=[0.0, 0.001, 0.01, 0.075, 0.2, 10.0], maf_bins=[0.0, 0.05, 0.1, 0.2, 0.3, 0.500001])
elif anno_method == 'borzoi_magnitude_binsXdist_to_tss_bins':
	create_annotation_using_borzoi_magnitude_bins_x_dist_to_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, magnitude_bins=[0.0, 0.001, 0.01, 0.075, 0.2, 10.0], dist_bins=[0.0, 5000.0, 50000.0, 100001.0])
