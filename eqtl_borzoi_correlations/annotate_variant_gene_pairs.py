import numpy as np
import sys
import pdb
import gzip


def write_annotation_header(data, t, anno_names):
	t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
	t.write('\t' + '\t'.join(anno_names) + '\n')


def extract_annotation_names_and_sources(annotation_name_file):
	# Columns (tab-separated, with header): anno_name  source
	# source is either 'custom' or 'sldsc'
	anno_names = []
	anno_sources = []
	f = open(annotation_name_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_name = data[0]
		anno_source = data[1]
		if anno_source != 'custom' and anno_source != 'sldsc':
			print('assumption eroror: unknown annotation source ' + anno_source)
			pdb.set_trace()
		if anno_name in anno_names:
			print('assumption eroror: repeat annotation name ' + anno_name)
			pdb.set_trace()
		anno_names.append(anno_name)
		anno_sources.append(anno_source)
	f.close()
	return anno_names, anno_sources


def create_one_hot_category_names(bins, bin_prefix):
	category_names = []
	for bin_iter in range(len(bins)-1):
		category_names.append(bin_prefix + str(bin_iter))
	return category_names


def create_interaction_category_names(magnitude_bins, other_bins, other_bin_prefix):
	category_names = []
	for magnitude_bin_iter in range(len(magnitude_bins)-1):
		for other_bin_iter in range(len(other_bins)-1):
			category_names.append('magnitude_bin' + str(magnitude_bin_iter) + 'X' + other_bin_prefix + str(other_bin_iter))
	return category_names


def create_custom_annotation_config(anno_name):
	# A custom annotation is a named (binning function, bin edges) pair.
	# 'kind' tells the row loop how to turn a variant-gene pair into a category index.
	anno_config = {}
	anno_config['anno_name'] = anno_name
	anno_config['source'] = 'custom'

	if anno_name == 'intercept':
		anno_config['kind'] = 'intercept'
		anno_config['category_names'] = ['intercept']
	elif anno_name == 'borzoi_magnitude_bins':
		anno_config['kind'] = 'borzoi_magnitude'
		anno_config['bins'] = [0.0, 0.001, 0.01, 0.075, 0.2, 10.0]
		anno_config['category_names'] = create_one_hot_category_names(anno_config['bins'], 'magnitude_bin')
	elif anno_name == 'borzoi_effect_size_bins':
		anno_config['kind'] = 'borzoi_effect_size'
		anno_config['bins'] = [-10.0, -0.2, -0.075, -0.01, -0.001, 0.0, 0.001, 0.01, 0.075, 0.2, 10.0]
		anno_config['category_names'] = create_one_hot_category_names(anno_config['bins'], 'effect_size_bin')
	elif anno_name == 'borzoi_finer_effect_size_bins':
		anno_config['kind'] = 'borzoi_effect_size'
		anno_config['bins'] = [-10.0, -0.25, -.1, -0.04, -0.01, -0.001, 0.0, 0.001, 0.01, 0.04, .1, 0.25, 10.0]
		anno_config['category_names'] = create_one_hot_category_names(anno_config['bins'], 'effect_size_bin')
	elif anno_name == 'dist_to_tss_bins':
		anno_config['kind'] = 'vg_pair_info'
		anno_config['vg_pair_info_index'] = 1
		anno_config['bins'] = [0.0, 1000.0, 5000.0, 25000.0, 50000.0, 100001.0]
		anno_config['category_names'] = create_one_hot_category_names(anno_config['bins'], 'dist_to_tss_bin')
	elif anno_name == 'strand_dist_to_tss_bins':
		anno_config['kind'] = 'vg_pair_info'
		anno_config['vg_pair_info_index'] = 2
		anno_config['bins'] = [-100001.0, -50000.0, -25000.0, -5000.0, -1000.0, 0.0, 1000.0, 5000.0, 25000.0, 50000.0, 100001.0]
		anno_config['category_names'] = create_one_hot_category_names(anno_config['bins'], 'strand_based_dist_to_tss_bin')
	elif anno_name == 'af_bins':
		anno_config['kind'] = 'vg_pair_info'
		anno_config['vg_pair_info_index'] = 0
		anno_config['bins'] = [0.0, 0.05, 0.1, 0.2, 0.3, 0.500001]
		anno_config['category_names'] = create_one_hot_category_names(anno_config['bins'], 'maf_bin')
	elif anno_name == 'borzoi_magnitude_binsXaf_bins':
		anno_config['kind'] = 'borzoi_magnitudeXvg_pair_info'
		anno_config['vg_pair_info_index'] = 0
		anno_config['bins'] = [0.0, 0.001, 0.01, 0.075, 0.2, 10.0]
		anno_config['other_bins'] = [0.0, 0.05, 0.1, 0.2, 0.3, 0.500001]
		anno_config['category_names'] = create_interaction_category_names(anno_config['bins'], anno_config['other_bins'], 'maf_bin')
	elif anno_name == 'borzoi_magnitude_binsXdist_to_tss_bins':
		anno_config['kind'] = 'borzoi_magnitudeXvg_pair_info'
		anno_config['vg_pair_info_index'] = 1
		anno_config['bins'] = [0.0, 0.001, 0.01, 0.075, 0.2, 10.0]
		anno_config['other_bins'] = [0.0, 5000.0, 50000.0, 100001.0]
		anno_config['category_names'] = create_interaction_category_names(anno_config['bins'], anno_config['other_bins'], 'dist_to_tss_bin')
	else:
		print('assumption eroror: unknown custom annotation ' + anno_name)
		pdb.set_trace()

	return anno_config


def create_sldsc_annotation_config(anno_name, sldsc_index):
	# An s-ldsc annotation is binary, so it becomes a two category partition
	anno_config = {}
	anno_config['anno_name'] = anno_name
	anno_config['source'] = 'sldsc'
	anno_config['kind'] = 'sldsc'
	anno_config['sldsc_index'] = sldsc_index
	anno_config['category_names'] = ['absent', 'present']
	return anno_config


def create_annotation_configs(anno_names, anno_sources):
	anno_configs = []
	sldsc_anno_names = []
	for anno_iter, anno_name in enumerate(anno_names):
		if anno_sources[anno_iter] == 'custom':
			anno_config = create_custom_annotation_config(anno_name)
		else:
			anno_config = create_sldsc_annotation_config(anno_name, len(sldsc_anno_names))
			sldsc_anno_names.append(anno_name)
		anno_config['category_counts'] = np.zeros(len(anno_config['category_names']), dtype=int)
		anno_config['missing_count'] = 0
		anno_configs.append(anno_config)
	return anno_configs, sldsc_anno_names


def extract_bin_index(value, bins):
	for bin_iter in range(len(bins)-1):
		if value >= bins[bin_iter] and value < bins[bin_iter+1]:
			return bin_iter
	return -1


def extract_category_index(anno_config, data, vg_pair_info, variant_to_sldsc_index, sldsc_anno_mat):
	# Returns the index of the category this variant-gene pair falls in (-1 if the pair has no category)
	anno_kind = anno_config['kind']

	if anno_kind == 'intercept':
		return 0

	if anno_kind == 'borzoi_magnitude':
		category_index = extract_bin_index(np.abs(float(data[6])), anno_config['bins'])
		if category_index == -1:
			print('assumption eroror')
			pdb.set_trace()
		return category_index

	if anno_kind == 'borzoi_effect_size':
		category_index = extract_bin_index(float(data[6]), anno_config['bins'])
		if category_index == -1:
			print('assumption eroror')
			pdb.set_trace()
		return category_index

	if anno_kind == 'vg_pair_info':
		vg_name = data[1] + ':' + data[0]
		if vg_name not in vg_pair_info:
			return -1
		return extract_bin_index(vg_pair_info[vg_name][anno_config['vg_pair_info_index']], anno_config['bins'])

	if anno_kind == 'borzoi_magnitudeXvg_pair_info':
		vg_name = data[1] + ':' + data[0]
		if vg_name not in vg_pair_info:
			return -1
		magnitude_bin_index = extract_bin_index(np.abs(float(data[6])), anno_config['bins'])
		other_bin_index = extract_bin_index(vg_pair_info[vg_name][anno_config['vg_pair_info_index']], anno_config['other_bins'])
		if magnitude_bin_index == -1 or other_bin_index == -1:
			return -1
		return (magnitude_bin_index*(len(anno_config['other_bins'])-1)) + other_bin_index

	if anno_kind == 'sldsc':
		var_id = data[2] + ':' + data[3]
		if var_id not in variant_to_sldsc_index:
			return 0
		return int(sldsc_anno_mat[variant_to_sldsc_index[var_id], anno_config['sldsc_index']])

	print('assumption eroror: unknown annotation kind ' + anno_kind)
	pdb.set_trace()


def create_mapping_from_variant_to_sldsc_annos(baselineLD_anno_dir, sldsc_anno_names):
	# Pull every requested s-ldsc annotation column in a single pass through the baselineLD files
	variant_to_sldsc_index = {}
	sldsc_anno_mat = []

	if len(sldsc_anno_names) == 0:
		return variant_to_sldsc_index, None

	for chrom_num in range(1,23):
		print('sldsc chrom' + str(chrom_num))
		anno_file = baselineLD_anno_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'
		f = gzip.open(anno_file,'rt')
		head_count = 0
		chrom_anno_mat = []
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				anno_indices = []
				for sldsc_anno_name in sldsc_anno_names:
					tmp = np.where(np.asarray(data) == sldsc_anno_name)[0]
					if len(tmp) != 1:
						print('assumption oeoror: ' + sldsc_anno_name + ' not found in ' + anno_file)
						pdb.set_trace()
					anno_indices.append(tmp[0])
				continue

			var_id = data[0] + ':' + data[1]
			if var_id in variant_to_sldsc_index:
				continue

			anno_values = []
			for anno_index in anno_indices:
				anno_value = data[anno_index]
				if anno_value != '0' and anno_value != '1':
					print('assumptioneroorr: assume binary anno')
					pdb.set_trace()
				anno_values.append(int(anno_value))

			variant_to_sldsc_index[var_id] = len(variant_to_sldsc_index)
			chrom_anno_mat.append(anno_values)
		f.close()
		sldsc_anno_mat.append(np.asarray(chrom_anno_mat, dtype=np.int8).reshape(len(chrom_anno_mat), len(sldsc_anno_names)))

	sldsc_anno_mat = np.vstack(sldsc_anno_mat)
	print(str(sldsc_anno_mat.shape[0]) + ' variants with s-ldsc annotations')
	return variant_to_sldsc_index, sldsc_anno_mat


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


def create_annotation_file(borzoi_effect_file, borzoi_annotation_file, anno_configs, vg_pair_info, variant_to_sldsc_index, sldsc_anno_mat):
	# One column per annotation, holding the index of the category the variant-gene pair falls in.
	# -1 means the variant-gene pair falls in no category of that annotation.
	anno_names = []
	for anno_config in anno_configs:
		anno_names.append(anno_config['anno_name'])

	f = gzip.open(borzoi_effect_file,'rt')
	t = gzip.open(borzoi_annotation_file, 'wt')
	head_count = 0
	vg_pair_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			write_annotation_header(data, t, anno_names)
			continue

		if data[4] == data[5]:
			print('assumption eroroor')
			pdb.set_trace()

		category_indices = []
		for anno_config in anno_configs:
			category_index = extract_category_index(anno_config, data, vg_pair_info, variant_to_sldsc_index, sldsc_anno_mat)
			if category_index == -1:
				anno_config['missing_count'] = anno_config['missing_count'] + 1
			else:
				anno_config['category_counts'][category_index] = anno_config['category_counts'][category_index] + 1
			category_indices.append(str(category_index))

		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
		t.write('\t' + '\t'.join(category_indices) + '\n')
		vg_pair_count = vg_pair_count + 1

	f.close()
	t.close()

	print(str(vg_pair_count) + ' variant-gene pairs annotated')
	return


def create_annotation_category_file(annotation_category_file, anno_configs):
	# The (annotation, category) pairs found in the annotation file
	t = open(annotation_category_file, 'w')
	t.write('anno_name\tsource\tcategory_index\tcategory_name\n')
	for anno_config in anno_configs:
		for category_iter, category_name in enumerate(anno_config['category_names']):
			t.write(anno_config['anno_name'] + '\t' + anno_config['source'] + '\t' + str(category_iter) + '\t' + category_name + '\n')
	t.close()
	return


def print_annotation_category_counts(anno_configs):
	for anno_config in anno_configs:
		print(anno_config['anno_name'])
		for category_iter, category_name in enumerate(anno_config['category_names']):
			print('\t' + category_name + ': ' + str(anno_config['category_counts'][category_iter]))
		if anno_config['missing_count'] > 0:
			print('\tno category: ' + str(anno_config['missing_count']))
	return


#####################
# Command line args
#####################
borzoi_effect_file = sys.argv[1]
annotation_name_file = sys.argv[2]
borzoi_annotation_file = sys.argv[3]
eqtl_sumstats_file = sys.argv[4]
tissue_name = sys.argv[5]
baselineLD_anno_dir = sys.argv[6]

# Describes the (annotation, category) pairs making up the annotation file
annotation_category_file = borzoi_annotation_file.split('.txt.gz')[0] + '_categories.txt'


# Extract the annotations we want to include
anno_names, anno_sources = extract_annotation_names_and_sources(annotation_name_file)
anno_configs, sldsc_anno_names = create_annotation_configs(anno_names, anno_sources)
print(str(len(anno_configs)) + ' annotations (' + str(len(sldsc_anno_names)) + ' of them from s-ldsc)')

# Create mapping from variant-gene pairs to auxiliary info
vg_pair_info = extract_info_on_each_variant_gene_pair(eqtl_sumstats_file)



# Create mapping from variant to each of its s-ldsc annotations
variant_to_sldsc_index, sldsc_anno_mat = create_mapping_from_variant_to_sldsc_annos(baselineLD_anno_dir, sldsc_anno_names)

# Print annotation file
create_annotation_file(borzoi_effect_file, borzoi_annotation_file, anno_configs, vg_pair_info, variant_to_sldsc_index, sldsc_anno_mat)

# Print (annotation, category) pair file
create_annotation_category_file(annotation_category_file, anno_configs)

print_annotation_category_counts(anno_configs)
