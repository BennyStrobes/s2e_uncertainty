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
			print('assumption erororo')
			pdb.set_trace()

		value = vg_pair_info[vg_name][vg_pair_info_index]
		bin_membership = np.zeros(len(bins)-1)
		for bin_iter in range(len(bins)-1):
			if value >= bins[bin_iter] and value < bins[bin_iter+1]:
				bin_membership[bin_iter] = 1.0
				bin_counts[bin_iter] = bin_counts[bin_iter] + 1
				break

		if np.sum(bin_membership) != 1.0:
			print('assumption eroror')
			print(vg_name)
			pdb.set_trace()

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
			t.write('\tmagnitude_bin0\tmagnitude_bin1\tmagnitude_bin2\tmagnitude_bin3\n')
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


def create_annotation_using_distance_to_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, bins=[0.0, 1000.0, 5000.0, 10000.0, 50000.0, 100001.0]):
	return create_one_hot_annotation_from_vg_pair_mapping(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, 1, bins, 'dist_to_tss_bin')


def create_annotation_using_maf_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, bins=[0.0, 0.05, 0.1, 0.2, 0.3, 0.500001]):
	return create_one_hot_annotation_from_vg_pair_mapping(borzoi_effect_file, borzoi_annotation_file, vg_pair_info, 0, bins, 'maf_bin')


def extract_info_on_each_variant_gene_pair(eqtl_sumstats_dir, tissue_name):
	dicti = {}
	for chrom_num in range(1,23):
		filer = eqtl_sumstats_dir + tissue_name + '.v10.allpairs.chr' + str(chrom_num) + '.parquet'
		pf = pq.ParquetFile(filer)

		for rg in range(pf.num_row_groups):
			table = pf.read_row_group(
				rg,
				columns=['gene_id', 'variant_id', 'tss_distance', 'af', 'slope', 'slope_se', 'ma_count']
			)

			if table.num_columns == 0:
				continue

			gene_ids = table['gene_id'].to_pylist()
			variant_ids = table['variant_id'].to_pylist()
			tss_dists = table['tss_distance'].to_pylist()
			afs = table['af'].to_pylist()

			for ii, gene_id in enumerate(gene_ids):
				maf = afs[ii]
				if maf > 0.5:
					maf = 1.0 - maf
				if maf == 0.0:
					continue

				abs_tss_dist = np.abs(tss_dists[ii])
				if abs_tss_dist > 100000:
					continue

				variant_id = variant_ids[ii]
				vg_name = variant_id + ':' + gene_id.split('.')[0]
				dicti[vg_name] = (maf, abs_tss_dist)
	return dicti


#####################
# Command line args
#####################
borzoi_effect_file = sys.argv[1]
anno_method = sys.argv[2]
borzoi_annotation_file = sys.argv[3]
eqtl_sumstats_dir = sys.argv[4]
tissue_name = sys.argv[5]

# Create mapping from variant-gene pairs to auxiliary info
vg_pair_info = extract_info_on_each_variant_gene_pair(eqtl_sumstats_dir, tissue_name)


if anno_method == 'borzoi_magnitude_bins':
	create_annotation_using_borzoi_magnitude_bins(borzoi_effect_file, borzoi_annotation_file)
elif anno_method == 'dist_to_tss_bins':
	create_annotation_using_distance_to_tss_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info)
elif anno_method == 'af_bins':
	create_annotation_using_maf_bins(borzoi_effect_file, borzoi_annotation_file, vg_pair_info)
