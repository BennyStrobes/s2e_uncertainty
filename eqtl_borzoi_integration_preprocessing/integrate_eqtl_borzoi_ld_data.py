import numpy as np
import os
import sys
import pdb
import pyarrow.parquet as pq
import pandas as pd
import pyarrow.compute as pc
from pgenlib import PgenReader
import h5py



def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf, chrom_num):
	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != 'chr' + str(chrom_num):
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = 1

	f.close()

	return dicti





def extract_ordered_gtex_tissues(gtex_sumstats_dir, chrom_num):
	arr = []
	for file_name in os.listdir(gtex_sumstats_dir):
		if not file_name.endswith('chr' + str(chrom_num) + '.parquet'):
			continue
		tissue_name = file_name.split('.v10')[0]
		arr.append(tissue_name)
	return np.sort(np.unique(arr))



def load_in_eqtl_sumstats(gtex_sumstats_dir, chrom_num, ordered_gtex_tissues, pc_genes):

	tissue_mapping = {}
	for ii, tissue in enumerate(ordered_gtex_tissues):
		tissue_mapping[tissue] = ii

	qtl_obj = {}
	for file_name in os.listdir(gtex_sumstats_dir):
		if not file_name.endswith('chr' + str(chrom_num) + '.parquet'):
			continue
		tissue_name = file_name.split('.v10')[0]
		print(tissue_name)
		pf = pq.ParquetFile(gtex_sumstats_dir + file_name)

		for rg in range(pf.num_row_groups):
			table = pf.read_row_group(
				rg,
				columns=['gene_id', 'variant_id', 'tss_distance', 'af', 'slope', 'slope_se', 'ma_count']
			)

			if table.num_columns == 0:
				continue


			gene_col = table['gene_id']
			variant_col = table['variant_id']
			dist_col = table['tss_distance']
			af_col = table['af']
			slope_col = table['slope']
			slope_se_col = table['slope_se']


			gene_ids = table['gene_id'].to_pylist()
			variant_ids = table['variant_id'].to_pylist()
			tss_dists = table['tss_distance'].to_pylist()
			afs = table['af'].to_pylist()
			ma_counts = table['ma_count'].to_pylist()
			slopes = table['slope'].to_pylist()
			slope_ses = table['slope_se'].to_pylist()


			for ii, gene_id in enumerate(gene_ids):
				if gene_id.split('.')[0] not in pc_genes:
					continue

				maf = afs[ii]
				if maf > 0.5:
					maf = 1.0 - maf
				if maf == 0.0:
					continue
				pred_N = (ma_counts[ii]/maf)/2.0

				round_pred_N = np.round(pred_N)

				line_tss_dist = tss_dists[ii]
				if np.abs(line_tss_dist) > 100000:
					continue
				variant_id = variant_ids[ii]
				var_info = variant_id.split('_')
				if len(var_info[2]) != 1 or len(var_info[3]) != 1:
					continue

				af = afs[ii]
				if slopes[ii] is None or slope_ses[ii] is None:
					continue
				zed = slopes[ii]/slope_ses[ii]

				if gene_id not in qtl_obj:
					qtl_obj[gene_id] = {}
				if variant_id not in qtl_obj[gene_id]:
					qtl_obj[gene_id][variant_id] = {}
					qtl_obj[gene_id][variant_id]['zeds']= np.asarray([np.nan]*len(ordered_gtex_tissues))
					qtl_obj[gene_id][variant_id]['N_eff']= np.asarray([np.nan]*len(ordered_gtex_tissues))
					qtl_obj[gene_id][variant_id]['distance'] = line_tss_dist
					qtl_obj[gene_id][variant_id]['af'] = af

				qtl_obj[gene_id][variant_id]['zeds'][tissue_mapping[tissue_name]] = zed
				qtl_obj[gene_id][variant_id]['N_eff'][tissue_mapping[tissue_name]] = round_pred_N

	return qtl_obj


def print_eqtl_summary(gene_eqtl_summary_file, eqtl_ss_obj, eqtl_ss_output_dir, chrom_num):
	t = open(gene_eqtl_summary_file,'w')
	t.write('gene_name\tvariant_info_file\teqtl_ss_matrix_npy_file\teQTL_N_matrix_npy_file\n')

	for gene_id in np.sort([*eqtl_ss_obj]):
		# Load in data for this gene
		gene_variants = []
		gene_afs = []
		gene_distances = []
		gene_zeds = []
		gene_Ns = []
		gene_snp_poss = []
		for var_id in np.sort([*eqtl_ss_obj[gene_id]]):
			gene_variants.append(var_id)
			gene_afs.append(eqtl_ss_obj[gene_id][var_id]['af'])
			gene_distances.append(eqtl_ss_obj[gene_id][var_id]['distance'])
			gene_zeds.append(eqtl_ss_obj[gene_id][var_id]['zeds'])
			gene_Ns.append(eqtl_ss_obj[gene_id][var_id]['N_eff'])
			gene_snp_poss.append(int(var_id.split('_')[1]))
		gene_variants = np.asarray(gene_variants)
		gene_afs = np.asarray(gene_afs)
		gene_distances = np.asarray(gene_distances)
		gene_zeds = np.asarray(gene_zeds)
		gene_Ns = np.asarray(gene_Ns)
		gene_snp_poss = np.asarray(gene_snp_poss)

		# Reorder
		snp_ordering = np.argsort(gene_snp_poss)
		gene_variants = gene_variants[snp_ordering]
		gene_afs = gene_afs[snp_ordering]
		gene_distances = gene_distances[snp_ordering]
		gene_zeds = gene_zeds[snp_ordering, :]
		gene_Ns = gene_Ns[snp_ordering, :]

		# Print snp summary file for gene
		gene_snp_summary_file = eqtl_ss_output_dir + 'chr' + str(chrom_num) + '_' + gene_id + '_snp_summary.txt'
		t2 = open(gene_snp_summary_file,'w')
		t2.write('chr\tsnp_id\tsnp_pos\ta1\ta2\ttss_dist\taf\n')
		for ii, var_id in enumerate(gene_variants):
			var_info = var_id.split('_')
			if var_info[4] != 'b38':
				print('assumption erororro')
				pdb.set_trace()
			if var_info[0] != 'chr' + chrom_num:
				print('assumptioneornron')
				pdb.set_trace()
			t2.write(var_info[0] + '\t' + var_id + '\t' + var_info[1] + '\t' + var_info[2] + '\t' + var_info[3] + '\t' + str(gene_distances[ii]) + '\t' + str(gene_afs[ii]) + '\n')
		t2.close()

		# Print zeds to output file
		gene_zed_npy_file = eqtl_ss_output_dir + 'chr' + str(chrom_num) + '_' + gene_id + '_eqtl_zed.npy'
		np.save(gene_zed_npy_file, gene_zeds)

		# Print Ns to output file
		gene_N_npy_file = eqtl_ss_output_dir + 'chr' + str(chrom_num) + '_' + gene_id + '_eqtl_N_eff.npy'
		np.save(gene_N_npy_file, gene_Ns)

		# Print to global output file
		t.write(gene_id + '\t' + gene_snp_summary_file + '\t' + gene_zed_npy_file + '\t' + gene_N_npy_file + '\n')

	t.close()
	return


def extract_plink_prefix(plink2_genotype_data_dir, chrom_num):
	default_prefixes = [
		os.path.join(plink2_genotype_data_dir, 'gtex_v9_eqtl_chr' + str(chrom_num)),
		os.path.join(plink2_genotype_data_dir, 'chr' + str(chrom_num)),
	]
	for default_prefix in default_prefixes:
		if os.path.isfile(default_prefix + '.pgen'):
			return default_prefix
	print('assumption error: unable to infer genotype prefix')
	pdb.set_trace()


def load_pvar_variant_info(genotype_data_prefix):
	pvar_file = genotype_data_prefix + '.pvar'
	header_line = None
	with open(pvar_file) as f:
		for ii, line in enumerate(f):
			if line.startswith('#CHROM'):
				header_line = ii
				break
	if header_line is None:
		print('assumption error: pvar header not found')
		pdb.set_trace()

	pvar_df = pd.read_csv(pvar_file, sep='\t', header=header_line)
	required_cols = np.asarray(['#CHROM', 'POS', 'REF', 'ALT'])
	if np.array_equal(np.isin(required_cols, pvar_df.columns), np.asarray([True, True, True, True])) == False:
		print('assumption error: pvar missing required columns')
		pdb.set_trace()


	variant_info = {}
	for ii, row in pvar_df.iterrows():
		chrom = str(row['#CHROM'])
		pos = str(row['POS'])
		ref = str(row['REF'])
		alts = str(row['ALT']).split(',')
		var_id = str(row['ID'])
		if var_id.split(':')[2] != ref:
			print('assumption erooror')
			pdb.set_trace()
		if var_id.split(':')[3] != alts[0]:
			print('assumption eroroor')
			pdb.set_trace()
		if len(alts) != 1:
			print('assumptino erororor')
			pdb.set_trace()
		for alt in alts:
			variant_key = ('chr' + chrom, pos, ref, alt)
			if variant_key not in variant_info:
				variant_info[variant_key] = []
			variant_info[variant_key].append(ii)
	return variant_info


def map_gtex_variant_to_pgen_index(gtex_variant_id, pvar_variant_info):
	var_info = gtex_variant_id.split('_')
	if len(var_info) < 5:
		print('assumption error: malformed gtex variant id')
		pdb.set_trace()

	chrom = var_info[0]
	pos = var_info[1]
	gtex_a1 = var_info[2]
	gtex_a2 = var_info[3]
	forward_key = (chrom, pos, gtex_a1, gtex_a2)
	reverse_key = (chrom, pos, gtex_a2, gtex_a1)

	if forward_key in pvar_variant_info:
		matches = pvar_variant_info[forward_key]
		if len(matches) != 1:
			print('assumption error: multiple forward matches in genotype data: ' + gtex_variant_id)
			pdb.set_trace()
		return matches[0], 1.0

	if reverse_key in pvar_variant_info:
		print('reverse')
		matches = pvar_variant_info[reverse_key]
		if len(matches) != 1:
			print('assumption error: multiple reverse matches in genotype data: ' + gtex_variant_id)
			pdb.set_trace()
		return matches[0], -1.0

	print('assumption error: allele mismatch between gtex and genotype data: ' + gtex_variant_id)
	pdb.set_trace()


def load_gene_genotype_from_pgen(genotype_data_prefix, gene_variant_indices):
	variant_indices = np.asarray(gene_variant_indices, dtype=np.uint32)
	reader = PgenReader((genotype_data_prefix + '.pgen').encode('utf-8'))
	sample_ct = reader.get_raw_sample_ct()
	genotype = np.empty((len(variant_indices), sample_ct), dtype=np.int8)
	try:
		reader.read_list(variant_indices, genotype)
	finally:
		reader.close()
	return genotype


def compute_ld_from_genotype_matrix(raw_genotype):
	genotype = np.asarray(raw_genotype, dtype=float)
	if len(genotype.shape) != 2:
		print('assumption error: genotype matrix is not 2d')
		pdb.set_trace()

	for ii in range(genotype.shape[0]):
		variant_values = genotype[ii, :]
		missing = variant_values == -9
		if np.sum(missing) == len(variant_values):
			print('assumption error: all genotype values missing')
			pdb.set_trace()
		if np.any(missing):
			var_missing_rate = np.sum(missing)/len(missing)
			if var_missing_rate > .17:
				print('high variant missing rate')
			variant_mean = np.mean(variant_values[missing == False])
			variant_values[missing] = variant_mean
		variant_sd = np.std(variant_values)
		if variant_sd == 0.0:
			print('assumption error: monomorphic variant encountered')
			pdb.set_trace()
		genotype[ii, :] = (variant_values - np.mean(variant_values))/variant_sd

	ld_mat = np.dot(genotype, np.transpose(genotype))/genotype.shape[1]
	return ld_mat


def compute_diagonal_padded_inverse_ld(ld_mat, diag_padding=1e-4):
	if ld_mat.shape[0] != ld_mat.shape[1]:
		print('assumption error: ld matrix is not square')
		pdb.set_trace()
	regularized_ld = np.copy(ld_mat)
	regularized_ld[np.diag_indices(regularized_ld.shape[0])] = regularized_ld[np.diag_indices(regularized_ld.shape[0])] + diag_padding
	return np.linalg.inv(regularized_ld)


def generate_gene_ld_summary(gene_eqtl_summary_file, plink2_genotype_data_dir, LD_output_dir, chrom_num):
	os.makedirs(LD_output_dir, exist_ok=True)
	print('start')
	gene_eqtl_summary_df = pd.read_csv(gene_eqtl_summary_file, sep='\t')
	genotype_data_prefix = extract_plink_prefix(plink2_genotype_data_dir, chrom_num)

	pvar_variant_info = load_pvar_variant_info(genotype_data_prefix)

	ld_summary_file = os.path.join(LD_output_dir, 'gene_ld_summary_chr' + str(chrom_num) + '.txt')
	t = open(ld_summary_file, 'w')
	t.write('gene_name\tvariant_info_file\teqtl_ss_matrix_npy_file\teqtl_N_npy_file\tld_matrix_npy_file\tn_variants\n')
	for _, row in gene_eqtl_summary_df.iterrows():
		gene_name = row['gene_name']
		print(gene_name)
		variant_info_file = row['variant_info_file']
		eqtl_ss_matrix_npy_file = row['eqtl_ss_matrix_npy_file']
		eqtl_N_eff_npy_file = row['eQTL_N_matrix_npy_file']
		variant_df = pd.read_csv(variant_info_file, sep='\t')
		gene_variant_ids = variant_df['snp_id'].astype(str).to_numpy()
		gene_variant_afs = variant_df['af'].astype(float).to_numpy()
		gene_variant_mafs = np.minimum(gene_variant_afs, 1.0-gene_variant_afs)
		gene_variant_indices = []
		gene_variant_signs = []
		for variant_id in gene_variant_ids:
			variant_index, variant_sign = map_gtex_variant_to_pgen_index(variant_id, pvar_variant_info)
			gene_variant_indices.append(variant_index)
			gene_variant_signs.append(variant_sign)
		gene_variant_indices = np.asarray(gene_variant_indices)
		gene_variant_signs = np.asarray(gene_variant_signs)

		common_indices = gene_variant_mafs >= 0.05
		if np.sum(common_indices) < 20:
			continue

		gene_genotype = load_gene_genotype_from_pgen(genotype_data_prefix, gene_variant_indices)
		ld_mat = compute_ld_from_genotype_matrix(gene_genotype)
		ld_mat = (gene_variant_signs[:, None]*ld_mat)*gene_variant_signs[None, :]
		ld_matrix_npy_file = os.path.join(LD_output_dir, gene_name + '_ld.npy')
		#ld_diagonal_padded_inverse_npy_file = os.path.join(LD_output_dir, gene_name + '_ld_diagonal_padded_inverse.npy')
		#ld_diagonal_padded_inverse = compute_diagonal_padded_inverse_ld(ld_mat[common_indices,:][:, common_indices], diag_padding=1e-3)
		np.save(ld_matrix_npy_file, ld_mat)
		#np.save(ld_diagonal_padded_inverse_npy_file, ld_diagonal_padded_inverse)
		t.write(gene_name + '\t' + variant_info_file + '\t' + eqtl_ss_matrix_npy_file + '\t' + eqtl_N_eff_npy_file + '\t' + ld_matrix_npy_file + '\t' + str(len(gene_variant_ids)) + '\n')
	t.close()
	return


def extract_borzoi_h5_files(borzoi_input_dir):
	arr = []

	for chunk_iter in range(20):
		filer = borzoi_input_dir + 'model_0_chunk_' + str(chunk_iter) + '_borzoi_results.h5'
		arr.append(filer)
	return arr


def collect_gene_borzoi_row_indices(borzoi_h5_files, valid_genes, chunk_size=10000):
	gene_to_rows = {}
	for gene_name in valid_genes:
		gene_to_rows[gene_name] = []

	for file_index, h5_file in enumerate(borzoi_h5_files):
		with h5py.File(h5_file, "r") as f:
			n_rows = f["gene"].shape[0]

			for start in range(0, n_rows, chunk_size):
				end = np.min([start + chunk_size, n_rows])
				gene_chunk = f["gene"][start:end].astype(str)
				snp_chunk = f["snp"][start:end].astype(str)

				for ii, gene_name in enumerate(gene_chunk):
					small_gene_name = gene_name.split('.')[0]
					if small_gene_name not in gene_to_rows:
						continue
					gene_to_rows[small_gene_name].append((file_index, start + ii, snp_chunk[ii]))

	return gene_to_rows


def generate_borzoi_summary_file(ld_summary_file, chrom_num, borzoi_input_dir, borzoi_output_dir):
	os.makedirs(borzoi_output_dir, exist_ok=True)
	ld_df = pd.read_csv(ld_summary_file, sep='\t')
	borzoi_h5_files = extract_borzoi_h5_files(borzoi_input_dir)
	if len(borzoi_h5_files) == 0:
		print('assumption error: no borzoi h5 files found')
		pdb.set_trace()

	gene_to_variant_index = {}
	gene_to_n_variants = {}
	for _, row in ld_df.iterrows():
		gene_name = row['gene_name'].split('.')[0]
		variant_df = pd.read_csv(row['variant_info_file'], sep='\t')
		ordered_variants = variant_df['snp_id'].astype(str).to_numpy()
		gene_to_variant_index[gene_name] = {}
		for ii, variant_id in enumerate(ordered_variants):
			gene_to_variant_index[gene_name][variant_id] = (ii, 1)
			var_info = variant_id.split('_')
			alt_variant_id = var_info[0] + '_' + var_info[1] + '_' + var_info[3] + '_' + var_info[2] + '_' + var_info[4]
			gene_to_variant_index[gene_name][alt_variant_id] = (ii, -1)
		gene_to_n_variants[gene_name] = len(ordered_variants)

	gene_to_rows = collect_gene_borzoi_row_indices(borzoi_h5_files, gene_to_variant_index)

	output_summary_file = os.path.join(borzoi_output_dir, 'gene_borzoi_summary_chr' + str(chrom_num) + '.txt')
	t = open(output_summary_file, 'w')
	t.write('gene_name\tvariant_info_file\teqtl_ss_matrix_npy_file\teqtl_N_npy_file\tld_matrix_npy_file\tborzoi_matrix_npy_file\tn_variants\tn_borzoi_matches\n')

	for _, row in ld_df.iterrows():
		gene_name = row['gene_name'].split('.')[0]
		print(gene_name)
		gene_rows = gene_to_rows[gene_name]
		grouped_rows = {}
		for file_index, row_index, snp_id in gene_rows:
			if file_index not in grouped_rows:
				grouped_rows[file_index] = []
			grouped_rows[file_index].append((row_index, snp_id))

		borzoi_mat = None
		n_borzoi_matches = 0
		for file_index in grouped_rows:
			h5_file = borzoi_h5_files[file_index]
			file_rows = np.asarray([*grouped_rows[file_index]])
			row_indices = file_rows[:, 0].astype(int)
			snp_ids = file_rows[:, 1].astype(str)

			with h5py.File(h5_file, "r") as f:
				log_ref = f["logRef"][row_indices, :]
				log_alt = f["logAlt"][row_indices, :]
				delta = log_alt - log_ref

			if borzoi_mat is None:
				n_targets = delta.shape[1]
				borzoi_mat = np.asarray(np.nan*np.ones((gene_to_n_variants[gene_name], n_targets)), dtype=np.float32)

			for ii, snp_id in enumerate(snp_ids):
				if snp_id not in gene_to_variant_index[gene_name]:
					print('missing snp error!')
					continue
				variant_row, effect_direction = gene_to_variant_index[gene_name][snp_id]
				if effect_direction == -1:
					print('negative effect direction')
				borzoi_mat[variant_row, :] = effect_direction*delta[ii, :]
				n_borzoi_matches = n_borzoi_matches + 1

		if borzoi_mat is None:
			print('missing')
			continue

		borzoi_matrix_npy_file = os.path.join(borzoi_output_dir, gene_name + '_borzoi.npy')
		np.save(borzoi_matrix_npy_file, borzoi_mat)
		t.write(
			gene_name + '\t' +
			row['variant_info_file'] + '\t' +
			row['eqtl_ss_matrix_npy_file'] + '\t' +
			row['eqtl_N_npy_file'] + '\t' +
			row['ld_matrix_npy_file'] + '\t' +
			borzoi_matrix_npy_file + '\t' +
			str(gene_to_n_variants[gene_name]) + '\t' +
			str(n_borzoi_matches) + '\n'
		)
	t.close()
	return

def generate_summary_file_ready_for_inv_ld(borzoi_summary_file, inverse_ld_extension_summary_file, LD_output_dir, eqtl_ss_output_dir, filter_LD=True):
	head_count = 0
	f = open(borzoi_summary_file)
	t = open(inverse_ld_extension_summary_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' 'inv_ld_matrix_npy_file' + '\t' + data[5] + '\t' + data[6] + '\t' + data[7] + '\n')
			print(data)
			continue
		gene_name = data[0]
		snp_summary_file = data[1]
		zed_file = data[2]
		N_eff_file = data[3]
		ld_file = data[4]
		borzoi_file = data[5]
		tot_snps = int(data[6])
		analysis_snps = int(data[7])

		# Load in data
		ld_mat = np.load(ld_file)
		zeds = np.load(zed_file)
		N_eff_file = np.load(N_eff_file)
		snp_afs = (np.loadtxt(snp_summary_file,dtype=str,delimiter='\t')[1:,-1]).astype(float)
		snp_mafs = np.minimum(snp_afs, 1.0-snp_afs)
		borzoi_mat = np.load(borzoi_file)
		valid_indices = np.where(~np.isnan(borzoi_mat).any(axis=1))[0]
		small_zeds = zeds[valid_indices, :]

		valid_gene = True
		if np.sum(np.isnan(small_zeds)) > 0.0:
			for tiss_iter in range(small_zeds.shape[1]):
				if np.sum(np.isnan(small_zeds[:,tiss_iter])) > 0.0:
					if np.sum(np.isnan(small_zeds[:,tiss_iter])) != len(small_zeds[:,tiss_iter]):
						valid_gene = False
		if valid_gene == False:
			print('skipped gene')
			continue

		if filter_LD:
			small_ld = ld_mat[valid_indices, :][:, valid_indices]
			inv_ld_file = LD_output_dir + gene_name + '_inv_ld_padded_1e-3.npy'
			if small_ld.shape[0] != analysis_snps:
				print('assumption error: mismatch between LD matrix and analysis SNPs')
				pdb.set_trace()
		else:
			small_ld = np.copy(ld_mat)
			analysis_snps = small_ld.shape[0]
			inv_ld_file = LD_output_dir + gene_name + '_big_inv_ld_padded_1e-3.npy'
		small_ld_inv = compute_diagonal_padded_inverse_ld(small_ld, diag_padding=1e-3)


		np.save(inv_ld_file, small_ld_inv)

		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + inv_ld_file + '\t' + data[5] + '\t' + data[6] + '\t' + str(analysis_snps) + '\n')


	f.close()
	t.close()
	return

########################
# Command line args
########################
borzoi_results_dir = sys.argv[1]
plink2_genotype_data_dir = sys.argv[2]
eqtl_sumstats_dir = sys.argv[3]
borzoi_output_dir = sys.argv[4] # outputdir
LD_output_dir = sys.argv[5] # outputdir
eqtl_ss_output_dir = sys.argv[6] # outputdir
chrom_num = sys.argv[7]
gtex_v10_pc_genes_gtf = sys.argv[8]


# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf, chrom_num)

# Extract ordered list of tissues
ordered_gtex_tissues = extract_ordered_gtex_tissues(eqtl_sumstats_dir, chrom_num)

##############
gtex_tissue_names_file=borzoi_output_dir + 'ordered_gtex_tissues_chr' + str(chrom_num) + '.txt'
# Print to output
t = open(gtex_tissue_names_file, 'w')
t.write('tissue_name\n')
for tissue_name in ordered_gtex_tissues:
	t.write(tissue_name + '\n')
t.close()


# Load in eqtl sumstats object
eqtl_ss_obj = load_in_eqtl_sumstats(eqtl_sumstats_dir, chrom_num, ordered_gtex_tissues, pc_genes)


# Print eQTL sumstat object to output
gene_eqtl_summary_file = eqtl_ss_output_dir + 'gene_summary_chr' + str(chrom_num) + '.txt'
print_eqtl_summary(gene_eqtl_summary_file, eqtl_ss_obj, eqtl_ss_output_dir, chrom_num)

# Generate LD summary file
ld_summary_file = os.path.join(LD_output_dir, 'gene_ld_summary_chr' + str(chrom_num) + '.txt')
generate_gene_ld_summary(gene_eqtl_summary_file, plink2_genotype_data_dir, LD_output_dir, chrom_num)


# Finally generate borzoi summary file
generate_borzoi_summary_file(ld_summary_file, chrom_num, borzoi_results_dir, borzoi_output_dir)
borzoi_summary_file = os.path.join(borzoi_output_dir, 'gene_borzoi_summary_chr' + str(chrom_num) + '.txt')



# Get inverse LD Ready
inverse_ld_extension_summary_file = os.path.join(borzoi_output_dir, 'gene_borzoi_summary_inv_ld_ready_chr' + str(chrom_num) + '.txt')
generate_summary_file_ready_for_inv_ld(borzoi_summary_file, inverse_ld_extension_summary_file, LD_output_dir, eqtl_ss_output_dir)



# Get inverse LD Ready
inverse_ld_extension_summary_file = os.path.join(borzoi_output_dir, 'gene_borzoi_summary_big_inv_ld_ready_chr' + str(chrom_num) + '.txt')
generate_summary_file_ready_for_inv_ld(borzoi_summary_file, inverse_ld_extension_summary_file, LD_output_dir, eqtl_ss_output_dir, filter_LD=False)



