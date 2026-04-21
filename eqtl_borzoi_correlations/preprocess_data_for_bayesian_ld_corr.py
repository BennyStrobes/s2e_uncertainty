import numpy as np
import os
import sys
import pdb
import pickle
import glob
from pandas_plink import read_plink
import pyarrow.parquet as pq
import pandas as pd
import pyarrow.compute as pc
import h5py





def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf):
	valid_chroms = {}
	for chrom_num in range(1,23):
		valid_chroms['chr' + str(chrom_num)] = 1


	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in valid_chroms:
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = 1

	f.close()

	return dicti





def load_in_eqtl_data(ordered_gtex_tissues, pc_genes, eqtl_sumstats_dir, chrom_num):
	mapping = {}
	for tissue_name in ordered_gtex_tissues:
		print(tissue_name)
		#for chrom_num in range(1,23):
		if True:
			filer = eqtl_sumstats_dir + tissue_name + '.v10.allpairs.chr' + str(chrom_num) + '.parquet'
			pf = pq.ParquetFile(filer)

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


				for ii, gene_id_full in enumerate(gene_ids):
					gene_id = gene_id_full.split('.')[0]
					if gene_id not in pc_genes:
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
					if slopes[ii] is None or slope_ses[ii] is None:
						continue
					if slope_ses[ii] == 0:
						continue

					desired_se = np.sqrt(1.0/round_pred_N)
					if gene_id not in mapping:
						mapping[gene_id] = {}
					if tissue_name not in mapping[gene_id]:
						mapping[gene_id][tissue_name] = {}
					zed = slopes[ii]/slope_ses[ii]
					std_effect = zed*desired_se
					if variant_id in mapping[gene_id][tissue_name]:
						print('repeat snp error')
						pdb.set_trace()
					mapping[gene_id][tissue_name][variant_id] = (gene_id, variant_id, chrom_num, var_info[1], var_info[2], var_info[3], std_effect, desired_se, maf)
					#mapping[gene_id][tissue_name][variant_id] = (gene_id, variant_id, chrom_num, var_info[1], var_info[2], var_info[3], std_effect, desired_se, np.copy(afs[ii]))

	return mapping


def save_pickle_object(pickle_file, data_obj):
	with open(pickle_file, 'wb') as f:
		pickle.dump(data_obj, f)


def load_pickle_object(pickle_file):
	with open(pickle_file, 'rb') as f:
		data_obj = pickle.load(f)
	return data_obj



def decode_if_bytes(arr):
	if len(arr) == 0:
		return arr
	if isinstance(arr[0], bytes):
		return np.asarray([x.decode('utf-8') for x in arr])
	return arr


def extract_chunk_num_from_borzoi_file(borzoi_input_file):
	file_name = os.path.basename(borzoi_input_file)
	return int(file_name.split('_chunk_')[1].split('_')[0])


def stream_borzoi_results(borzoi_results_dir, borzoi_target_indices, chunk_size=50000):
	file_pattern = os.path.join(borzoi_results_dir, 'model_0_chunk_*_borzoi_gtex_only_results.h5')
	borzoi_input_files = sorted(glob.glob(file_pattern), key=extract_chunk_num_from_borzoi_file)
	target_indices = np.asarray(borzoi_target_indices).astype(int)

	for borzoi_input_file in borzoi_input_files:
		print(borzoi_input_file)
		with h5py.File(borzoi_input_file, 'r') as f:
			snp_chrom_ds = f['snp_chrom']
			snp_pos_ds = f['snp_pos']
			snp_ds = f['snp']
			gene_ds = f['gene']
			logRef_ds = f['logRef']
			logAlt_ds = f['logAlt']

			n_rows = gene_ds.shape[0]
			for start in range(0, n_rows, chunk_size):
				end = min(start + chunk_size, n_rows)

				snp_chrom = decode_if_bytes(snp_chrom_ds[start:end])
				snp_pos = snp_pos_ds[start:end]
				snp_ids = decode_if_bytes(snp_ds[start:end])
				gene_ids = decode_if_bytes(gene_ds[start:end])
				logRef = logRef_ds[start:end][:, target_indices]
				logAlt = logAlt_ds[start:end][:, target_indices]
				borzoi_effects = logAlt - logRef

				yield snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects


def load_in_borzoi_data(eqtl_data_obj, borzoi_results_dir, borzoi_target_indices):
	mapping = {}
	for snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects in stream_borzoi_results(borzoi_results_dir, borzoi_target_indices):
		for chrom, pos, variant_id, gene_id, borzoi_effect_size_arr in zip(snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects):
			gene_name = gene_id.split('.')[0]
			if gene_name not in eqtl_data_obj:
				continue

			if gene_name not in mapping:
				mapping[gene_name] = {}

			var_info = variant_id.split('_')
			mapping[gene_name][variant_id] = (gene_name, variant_id, chrom.split('hr')[1], pos, var_info[2], var_info[3], borzoi_effect_size_arr)

	return mapping



def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	mapping = {}

	n_snps = len(ordered_snps)
	for snp_iter in range(n_snps):
		snp_name = ordered_snps[snp_iter]

		if snp_name in mapping:
			print('asssumption erororo')
			pdb.set_trace()
		mapping[snp_name] = snp_iter

	return mapping


def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
	if len(snp_array) != len(a0_arr):
		print('assumption eorroro')
		pdb.set_trace()
	if len(snp_array) != len(a1_arr):
		print('assumption eorroro')
		pdb.set_trace()

	dicti = {}

	for ii, snp_id in enumerate(snp_array):
		if snp_id in dicti:
			print('assumpationoenroer')
			pdb.set_trace()

		snp_info = snp_id.split('_')
		if snp_info[2] != a1_arr[ii]:
			print('error')
			pdb.set_trace()
		if snp_info[3] != a0_arr[ii]:
			print('errror')
			pdb.set_trace()

		dicti[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
	return dicti

def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num


def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_borzoi_effects, var_to_est_eqtl_effects, ordered_gtex_tissues):
	unique_vars = np.unique(np.hstack(([*var_to_est_borzoi_effects],[*var_to_est_eqtl_effects])))
	final_vars = []
	for var in unique_vars:
		if var not in rsid_to_genotype_index:
			continue
		geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])

		passing = False
		if var in var_to_est_borzoi_effects:
			passing = True
			borzoi_alleles = var_to_est_borzoi_effects[var][4:6]
			if set(geno_alleles) != set(borzoi_alleles):
				passing = False
		passing2 = False
		for tissue in ordered_gtex_tissues:
			if tissue in var_to_est_eqtl_effects and var in var_to_est_eqtl_effects[tissue]:
				passing2 = True
				eqtl_alleles = var_to_est_eqtl_effects[tissue][var][4:6]
				if set(geno_alleles) != set(eqtl_alleles):
					passing = False

		if passing == False:
			continue
		if passing2 == False:
			continue
		final_vars.append(var)
	return np.asarray(final_vars)

def load_in_snp_gene_eqtl_data(ordered_cis_variants, var_to_est_eqtl_effects, ordered_gtex_tissues):
	effects = []
	alleles = []
	effect_ses = []
	mafs = []

	for variant_id in ordered_cis_variants:
		var_effects = []
		var_effects_ses = []
		var_alleles = []
		var_mafs = []
		for tiss_iter, tissue in enumerate(ordered_gtex_tissues):
			if tissue not in var_to_est_eqtl_effects or variant_id not in var_to_est_eqtl_effects[tissue]:
				var_effects.append(np.nan)
				var_effects_ses.append(np.nan)
				var_mafs.append(np.nan)
				var_alleles.append('nan_nan')
			else:
				var_info = var_to_est_eqtl_effects[tissue][variant_id]
				var_effects.append(var_info[6])
				var_effects_ses.append(var_info[7])
				var_alleles.append(var_info[4] + '_' + var_info[5])
				var_mafs.append(var_info[8])
		var_effects = np.asarray(var_effects)
		var_effects_ses = np.asarray(var_effects_ses)
		var_alleles = np.asarray(var_alleles)
		var_mafs = np.asarray(var_mafs)

		effects.append(var_effects)
		effect_ses.append(var_effects_ses)
		alleles.append(var_alleles)
		mafs.append(var_mafs)
	return np.asarray(effects), np.asarray(alleles), np.asarray(effect_ses), np.asarray(mafs)


def load_in_borzoi_snp_gene_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []
	n_tracks = len(var_to_est_eqtl_effects[[*var_to_est_eqtl_effects][0]][6])

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.full(n_tracks, np.nan))
			alleles.append('nan_nan')
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			alleles.append(var_info[4] + '_' + var_info[5])
	return np.vstack(effects), np.asarray(alleles)

def generate_bayesian_ld_corr_input_data(gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, onek_genomes_plink_filestem, output_stem, ordered_gtex_tissues, chrom_num, t):

	# Loop through chromsomes
	#for chrom_num in range(1,23):
	if True:
		print(chrom_num)

		##################################
		# Load in per-chrom-genotype data
		##################################
		# string of chromosome name
		chrom_string = 'chr' + str(chrom_num)
		# Load in chromosome plink data
		(bim, fam, G) = read_plink(onek_genomes_plink_filestem + str(chrom_num))
		# Create mapping from variant id to index
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		# Create mapping from rsid to a0, a1
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))


		##################################
		# Loop through genes on this chromosome
		# (Analysis done seperately for each gene)
		##################################
		for gene_id in [*gene_id_to_est_borzoi_effects]:

			# Limit to genes on this chromosome
			gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
			if str(gene_chrom_num) != str(chrom_num):
				continue

			# Gene needs both borzoi effects AND eQTLs
			if gene_id not in gene_id_to_est_eqtl_effects:
				continue

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id], gene_id_to_est_eqtl_effects[gene_id], ordered_gtex_tissues)

			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 10:
				continue

			# Load in data for gene
			# eQTL
			eqtl_effects, eqtl_variant_alleles, eqtl_effect_ses, eqtl_mafs = load_in_snp_gene_eqtl_data(ordered_cis_variants, gene_id_to_est_eqtl_effects[gene_id], ordered_gtex_tissues)

			# Borzoi
			borzoi_effects, borzoi_variant_alleles = load_in_borzoi_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])

			# Load in LD
			cis_genotype_indices = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				geno_allele_string = geno_alleles[0] + '_' + geno_alleles[1]
				geno_allele_alt_string = geno_alleles[1] + '_' + geno_alleles[0]
				
				# Also flip signs of eqtls to match LD
				if np.isnan(borzoi_effects[var_index,0]) == False:
					if borzoi_variant_alleles[var_index] == geno_allele_alt_string:
						borzoi_effects[var_index,:] = -1.0*borzoi_effects[var_index,:]
						#print('flip!')
					else:
						print('no flip!')
				for tissue_index, tissue_name in enumerate(ordered_gtex_tissues):
					if np.isnan(eqtl_effects[var_index, tissue_index]) == False:
						if eqtl_variant_alleles[var_index, tissue_index] == geno_allele_alt_string:
							eqtl_effects[var_index, tissue_index] = -1.0*eqtl_effects[var_index, tissue_index]
							#print('flip2! ' + str(var_index))
						else:
							print('no flip!')

			
			# Extract genotype
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = G[cis_genotype_indices,:].compute()




			row_means = np.nanmean(geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(geno_mat))
			geno_mat[nan_rows, nan_cols] = row_means[nan_rows]

			LD = np.corrcoef(geno_mat)

			# Save files
			np.save(output_stem + '_' + gene_id + '_eqtl_effects.npy', eqtl_effects)
			np.save(output_stem + '_' + gene_id + '_eqtl_effect_ses.npy', eqtl_effect_ses)
			np.save(output_stem + '_' + gene_id + '_eqtl_mafs.npy', eqtl_mafs)
			np.save(output_stem + '_' + gene_id + '_borzoi_effects.npy', borzoi_effects)
			np.save(output_stem + '_' + gene_id + '_LD.npy', LD)

			# Print to global output file
			t.write(gene_id + '\t')
			t.write(output_stem + '_' + gene_id + '_eqtl_effects.npy' + '\t')
			t.write(output_stem + '_' + gene_id + '_eqtl_effect_ses.npy' + '\t')
			t.write(output_stem + '_' + gene_id + '_eqtl_mafs.npy' + '\t')
			t.write(output_stem + '_' + gene_id + '_borzoi_effects.npy' + '\t')
			t.write(output_stem + '_' + gene_id + '_LD.npy' + '\n')

	return t




#############################
# Command line args
#############################
borzoi_gtex_independent_target_names_file = sys.argv[1]
gtex_v10_pc_genes_gtf = sys.argv[2]
eqtl_sumstats_dir = sys.argv[3]
borzoi_results_dir = sys.argv[4]
genotype_stem = sys.argv[5]
bayes_input_data_stem = sys.argv[6] # OUTPUT STEM


# Load in target file
target_df = np.loadtxt(borzoi_gtex_independent_target_names_file,dtype=str,delimiter='\t')
ordered_gtex_tissues = target_df[1:,-1]
borzoi_target_indices = target_df[1:,1].astype(int)

# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)



# Global output file 
global_output_file = bayes_input_data_stem + '_per_gene_summary.txt'
t = open(global_output_file,'w')
t.write('gene_id\teqtl_effect_file\teqtl_effect_se_file\teqtl_maf_file\tborzoi_effect_file\tLD_file\n')



for chrom_num in range(1,23):


	###**#*#**#*#*# UPDATE CHROMS HERE
	# Extract eqtl data
	eqtl_data_obj = load_in_eqtl_data(ordered_gtex_tissues, pc_genes, eqtl_sumstats_dir, chrom_num)
	#eqtl_pickle_file = bayes_input_data_stem + '_eqtl_data.pkl'
	#save_pickle_object(eqtl_pickle_file, eqtl_data_obj)
	#eqtl_data_obj = load_pickle_object(eqtl_pickle_file)


	# Load in Borzoi data
	borzoi_data_obj = load_in_borzoi_data(eqtl_data_obj, borzoi_results_dir, borzoi_target_indices)
	#borzoi_pickle_file = bayes_input_data_stem + '_borzoi_data.pkl'
	#save_pickle_object(borzoi_pickle_file, borzoi_data_obj)
	#borzoi_data_obj = load_pickle_object(borzoi_pickle_file)


	# Get all organized together /save to output
	t = generate_bayesian_ld_corr_input_data(borzoi_data_obj, eqtl_data_obj, genotype_stem, bayes_input_data_stem, ordered_gtex_tissues, chrom_num, t)

t.close()



