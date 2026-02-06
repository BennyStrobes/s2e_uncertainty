import numpy as np
import os
import sys
import pdb
from bgen.reader import BgenFile
import gzip



def extract_list_of_trait_names(trait_names_file, chrom_num):
	f = open(trait_names_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[3] != 'chr' + chrom_num:
			continue
		arr.append(data[0])
	f.close()

	return np.unique(np.asarray(arr))



def generate_dictionary_list_of_wb_unrelated_samples(wb_unrelated_samples_file):
	f = open(wb_unrelated_samples_file)
	dicti = {}
	arr = []
	ii = 0
	for line in f:
		line = line.rstrip()
		data = line.split(' ')
		word = data[0] + '\t' + data[1]
		if word in dicti:
			print('repeat erororro')
			pdb.set_trace()
		dicti[word] = ii
		ii = ii + 1
		arr.append(word)
	f.close()

	return dicti, np.asarray(arr)


def create_mapping_from_genotyped_samples_to_unrelated_wb_ancestry_samples(ordered_wb_unrelated_samples, genotype_sample_names):
	# First create mapping from genotype sample to genotype position
	g_sample_to_g_position = {}
	g_samples_arr = []
	head_count = 0
	ii = 0
	f = open(genotype_sample_names)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count < 2:
			head_count = head_count + 1
			continue
		word = data[0] + '\t' + data[1]
		if word in g_sample_to_g_position:
			print('repeate rororor')
			pdb.set_trace()
		g_sample_to_g_position[word] = ii
		g_samples_arr.append(word)
		ii = ii + 1
	f.close()
	g_samples_arr = np.asarray(g_samples_arr)


	genotype_sample_indexer = []
	for sample_name in ordered_wb_unrelated_samples:
		genotype_sample_indexer.append(g_sample_to_g_position[sample_name])
	genotype_sample_indexer = np.asarray(genotype_sample_indexer)


	# Quick error check
	if np.array_equal(g_samples_arr[genotype_sample_indexer], ordered_wb_unrelated_samples) == False:
		print('assumptione ornorr')
		pdb.set_trace()

	return genotype_sample_indexer


def create_mapping_from_trait_name_to_pheno_file(trait_names_file):
	dicti = {}
	head_count = 0
	f = open(trait_names_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_name_orig = data[1]
		pheno_file = data[2]
		trait_name_pheno_file = data[3]

		if trait_name_orig in dicti:
			print('trait repeate erooror')
			pdb.set_trace()

		dicti[trait_name_orig] = (pheno_file, trait_name_pheno_file)

	f.close()
	return dicti


def extract_pheno_vec_v4(pheno_file_trait_name, ordered_wb_unrelated_samples, ukbb_pheno_file_v4):
	sample_names = []
	trait_vec = []
	f = open(ukbb_pheno_file_v4)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 32:
			print('assumption eorroorro')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			indices = np.where(np.asarray(data) == pheno_file_trait_name)[0]
			if len(indices) != 1:
				print('assumption error')
				pdb.set_trace()
			trait_index = indices[0]
			continue
		sample_name = data[0] + '\t' + data[1]
		sample_names.append(sample_name)
		if data[trait_index] == 'NA':
			trait_vec.append(np.nan)
		else:
			trait_vec.append(float(data[trait_index]))
	f.close()
	sample_names = np.asarray(sample_names)
	trait_vec = np.asarray(trait_vec)

	mapping = {}
	for ii, sample_name in enumerate(sample_names):
		mapping[sample_name] = ii


	reordered_trait = []
	for sample in ordered_wb_unrelated_samples:
		if sample in mapping:
			index = mapping[sample]
			reordered_trait.append(trait_vec[index])
		else:
			reordered_trait.append(np.nan)
	reordered_trait = np.asarray(reordered_trait)

	return reordered_trait


def extract_pheno_vec_v2(pheno_file_trait_name, ordered_wb_unrelated_samples, ukbb_pheno_file_v2):
	sample_names = []
	trait_vec = []
	f = open(ukbb_pheno_file_v2)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 45:
			print('assumption eorroorro')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			indices = np.where(np.asarray(data) == pheno_file_trait_name)[0]
			if len(indices) != 1:
				print('assumption error')
				pdb.set_trace()
			trait_index = indices[0]
			continue
		sample_name = data[0] + '\t' + data[1]
		sample_names.append(sample_name)
		if data[trait_index] == 'NA':
			trait_vec.append(np.nan)
		else:
			trait_vec.append(float(data[trait_index]))
	f.close()
	sample_names = np.asarray(sample_names)
	trait_vec = np.asarray(trait_vec)

	mapping = {}
	for ii, sample_name in enumerate(sample_names):
		mapping[sample_name] = ii


	reordered_trait = []
	for sample in ordered_wb_unrelated_samples:
		if sample in mapping:
			index = mapping[sample]
			reordered_trait.append(trait_vec[index])
		else:
			reordered_trait.append(np.nan)
	reordered_trait = np.asarray(reordered_trait)

	return reordered_trait


def extract_pheno_vec_v3(pheno_file_trait_name, ordered_wb_unrelated_samples, ukbb_pheno_file_v3):
	sample_names = []
	trait_vec = []
	f = gzip.open(ukbb_pheno_file_v3,'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 112:
			print('assumption eorroorro')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			indices = np.where(np.asarray(data) == pheno_file_trait_name)[0]
			if len(indices) != 1:
				print('assumption error')
				pdb.set_trace()
			trait_index = indices[0]
			continue
		sample_name = data[0] + '\t' + data[1]
		sample_names.append(sample_name)
		if data[trait_index] == 'NA':
			trait_vec.append(np.nan)
		else:
			trait_vec.append(float(data[trait_index]))
	f.close()
	sample_names = np.asarray(sample_names)
	trait_vec = np.asarray(trait_vec)

	mapping = {}
	for ii, sample_name in enumerate(sample_names):
		mapping[sample_name] = ii


	reordered_trait = []
	for sample in ordered_wb_unrelated_samples:
		if sample in mapping:
			index = mapping[sample]
			reordered_trait.append(trait_vec[index])
		else:
			reordered_trait.append(np.nan)
	reordered_trait = np.asarray(reordered_trait)

	return reordered_trait



##########################
# Command line args
##########################
ukbb_pheno_file_v2 = sys.argv[1]
ukbb_pheno_file_v3 = sys.argv[2]
ukbb_pheno_file_v4 = sys.argv[3]
organized_pred_results_file = sys.argv[4]
burden_genes_summary_file = sys.argv[5]
genotype_dir = sys.argv[6]
trait_names_file = sys.argv[7]
chrom_num = sys.argv[8]
genotype_sample_names = sys.argv[9]
wb_unrelated_samples_file = sys.argv[10]


##############
# First extract list of trait names (on this chromosome)
trait_names = extract_list_of_trait_names(burden_genes_summary_file, chrom_num)

# Create mapping from trait name to pheno file
trait_name_to_pheno_file_class = create_mapping_from_trait_name_to_pheno_file(trait_names_file)

# Create dictionary list of wb unrelated samples (need to subset to these)
wb_unrelated_samples_to_pos, ordered_wb_unrelated_samples = generate_dictionary_list_of_wb_unrelated_samples(wb_unrelated_samples_file)

# Create mapping from genotyped samples to unrelated wb ancestry samples
genotype_sample_indexer = create_mapping_from_genotyped_samples_to_unrelated_wb_ancestry_samples(ordered_wb_unrelated_samples, genotype_sample_names)


##################
# Create mapping from trait name to ordered phenotype vector (ok if nans exist)
trait_to_pheno_vec = {}
for trait_name in trait_names:

	# Where is data stored
	pheno_file_class, pheno_file_trait_name = trait_name_to_pheno_file_class[trait_name]


	if pheno_file_class == 'v2':
		pheno_vec = extract_pheno_vec_v2(pheno_file_trait_name, ordered_wb_unrelated_samples, ukbb_pheno_file_v2)
	elif pheno_file_class == 'v3':
		pheno_vec = extract_pheno_vec_v3(pheno_file_trait_name, ordered_wb_unrelated_samples, ukbb_pheno_file_v3)
	elif pheno_file_class == 'v4':
		pheno_vec = extract_pheno_vec_v4(pheno_file_trait_name, ordered_wb_unrelated_samples, ukbb_pheno_file_v4)
	else:
		print('should not be here')
		pdb.set_trace()
	

	if trait_name in trait_to_pheno_vec:
		print('erroro')
		pdb.set_trace()
	trait_to_pheno_vec[trait_name] = pheno_vec




# Stream variants









