import numpy as np
import os
import sys
import pdb





def extract_trait_and_tissue_from_matched_run_name(run_name):
	known_tissues = [
		'Skin_Sun_Exposed_Lower_leg',
		'Cells_Cultured_fibroblasts',
		'Whole_Blood',
		'Artery_Aorta',
		'Muscle_Skeletal',
		'Pituitary',
		'Spleen',
		'Liver',
		'Lung']
	run_prefix = run_name.split('_GTEX-')[0]
	for tissue_name in known_tissues:
		tissue_suffix = '_' + tissue_name
		if run_prefix.endswith(tissue_suffix):
			return run_prefix[:-len(tissue_suffix)], tissue_name
	print('Could not parse trait/tissue from matched run name: ' + run_name)
	return run_name, 'NA'


def run_matched_tissue_fsr_analysis(twas_dir, matched_runs, matched_tissue_fsr_summary_file):
	twas_zeds = []
	borzoi_zeds = []
	borzoi_probs = []
	hybrid_probs = []
	cis_h2s = []
	borzoi_fsrs = []
	t = open(matched_tissue_fsr_summary_file,'w')
	t.write('trait_name\ttissue\tgene\ttwas_zed\tborzoi_zed\tborzoi_fsr\tborzoi_prob\tcis_h2\n')
	for run_name in matched_runs:
		filer = twas_dir + 'twas_summary_' + run_name + '_borzoi_effect_sizes_ldscore_grid_squared_False.txt'
		trait_name, tissue_name = extract_trait_and_tissue_from_matched_run_name(run_name)

		run_gene_count = 0
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			twas_z = float(data[8])
			borzoi_twas_z = float(data[11])
			borzoi_prob = float(data[15])
			cis_h2_p = float(data[19])
			borzoi_fsr = float(data[17])
			hybrid_prob = float(data[14])
			gene_id = data[0]
			
			twas_zeds.append(twas_z)
			borzoi_zeds.append(borzoi_twas_z)
			borzoi_probs.append(borzoi_prob)
			cis_h2s.append(cis_h2_p)
			borzoi_fsrs.append(borzoi_fsr)
			hybrid_probs.append(hybrid_prob)
			run_gene_count = run_gene_count + 1
			t.write(trait_name + '\t' + tissue_name + '\t' + gene_id + '\t' + str(twas_z) + '\t' + str(borzoi_twas_z) + '\t' + str(borzoi_fsr) + '\t' + str(borzoi_prob) + '\t' + str(cis_h2_p) + '\n')
		f.close()

	twas_zeds = np.asarray(twas_zeds)
	borzoi_zeds = np.asarray(borzoi_zeds)
	borzoi_probs = np.asarray(borzoi_probs)
	cis_h2s = np.asarray(cis_h2s)
	borzoi_fsrs = np.asarray(borzoi_fsrs)
	hybrid_probs = np.asarray(hybrid_probs)

	hybrid_probs[np.isnan(hybrid_probs)] = 0.5
	borzoi_probs[np.isnan(borzoi_probs)] = 0.5

	t.close()
	return

def benjamini_hochberg_qvalues(pvalues):
	pvalues = np.asarray(pvalues)
	qvalues = np.ones(len(pvalues))
	valid_indices = np.where(np.isfinite(pvalues))[0]
	if len(valid_indices) == 0:
		return qvalues

	valid_pvalues = pvalues[valid_indices]
	order = np.argsort(valid_pvalues)
	ordered_pvalues = valid_pvalues[order]
	n_tests = len(ordered_pvalues)
	ordered_qvalues = ordered_pvalues*n_tests/(np.arange(n_tests) + 1.0)
	ordered_qvalues = np.minimum.accumulate(ordered_qvalues[::-1])[::-1]
	ordered_qvalues[ordered_qvalues > 1.0] = 1.0

	qvalues[valid_indices[order]] = ordered_qvalues
	return qvalues


def run_matched_tissue_borzoi_fsr_filtered_fdr_analysis(twas_dir, matched_runs, output_file):
	t = open(output_file,'w')
	t.write('trait_name\ttissue\ttwas_result_file\tn_genes\tn_borzoi_fsr_lt_0.1\tn_twas_fdr_lt_0.05\n')

	for run_name in matched_runs:
		twas_result_file = twas_dir + 'twas_summary_' + run_name + '_borzoi_effect_sizes_ldscore_grid_squared_False.txt'
		trait_name, tissue_name = extract_trait_and_tissue_from_matched_run_name(run_name)
		twas_pvalues = []
		n_genes = 0
		n_borzoi_fsr_filtered = 0

		f = open(twas_result_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			n_genes = n_genes + 1
			borzoi_fsr = float(data[17])
			if borzoi_fsr >= .1:
				continue
			n_borzoi_fsr_filtered = n_borzoi_fsr_filtered + 1
			twas_pvalues.append(float(data[10]))
		f.close()

		twas_qvalues = benjamini_hochberg_qvalues(twas_pvalues)
		n_twas_fdr_hits = np.sum(twas_qvalues < .05)
		t.write(trait_name + '\t' + tissue_name + '\t' + twas_result_file + '\t' + str(n_genes) + '\t' + str(n_borzoi_fsr_filtered) + '\t' + str(n_twas_fdr_hits) + '\n')

	t.close()
	print(output_file)


def run_matched_tissue_downsample_hit_summary(twas_dir, matched_runs, training_sample_sizes, output_file):
	anno_method = 'borzoi_effect_sizes'
	distribution = 'ldscore_grid_squared'
	standardize_geno = 'False'
	z_score_columns = [
		('standard_twas_z', 8),
		('hybrid_twas_z', 5)]
	abs_z_thresholds = np.asarray([3.0, 5.0, 10.0])

	t = open(output_file,'w')
	t.write('training_sample_size\tz_score_type\tabs_z_threshold\tcis_snp_h2_pvalue_threshold\tn_expected_trait_runs\tn_observed_trait_runs\tn_missing_trait_runs\tmean_n_tested_genes_per_trait\tse_n_tested_genes_per_trait\tmean_n_hit_genes_per_trait\tse_n_hit_genes_per_trait\n')

	for training_sample_size in training_sample_sizes:
		n_tested_gene_counts = []
		hit_counts = {}
		for z_score_type, column_index in z_score_columns:
			hit_counts[z_score_type] = {}
			for abs_z_threshold in abs_z_thresholds:
				hit_counts[z_score_type][abs_z_threshold] = []

		n_missing_trait_runs = 0
		for run_name in matched_runs:
			twas_result_file = twas_dir + 'twas_summary_' + run_name + '_' + anno_method + '_' + distribution + '_' + standardize_geno + '_downsample_' + str(training_sample_size) + '.txt'
			if os.path.exists(twas_result_file) == False:
				n_missing_trait_runs = n_missing_trait_runs + 1
				continue

			run_hit_counts = {}
			for z_score_type, column_index in z_score_columns:
				run_hit_counts[z_score_type] = {}
				for abs_z_threshold in abs_z_thresholds:
					run_hit_counts[z_score_type][abs_z_threshold] = 0

			run_n_tested_genes = 0
			f = open(twas_result_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				cis_snp_h2_pvalue = float(data[19])
				if cis_snp_h2_pvalue >= .05:
					continue
				run_n_tested_genes = run_n_tested_genes + 1
				for z_score_type, column_index in z_score_columns:
					z_score = float(data[column_index])
					for abs_z_threshold in abs_z_thresholds:
						if np.abs(z_score) >= abs_z_threshold:
							run_hit_counts[z_score_type][abs_z_threshold] = run_hit_counts[z_score_type][abs_z_threshold] + 1
			f.close()
			n_tested_gene_counts.append(run_n_tested_genes)

			for z_score_type, column_index in z_score_columns:
				for abs_z_threshold in abs_z_thresholds:
					hit_counts[z_score_type][abs_z_threshold].append(run_hit_counts[z_score_type][abs_z_threshold])

		for z_score_type, column_index in z_score_columns:
			for abs_z_threshold in abs_z_thresholds:
				counts = np.asarray(hit_counts[z_score_type][abs_z_threshold])
				tested_counts = np.asarray(n_tested_gene_counts)
				n_observed_trait_runs = len(counts)
				if n_observed_trait_runs == 0:
					mean_n_tested_genes = 'NA'
					se_n_tested_genes = 'NA'
					mean_n_hit_genes = 'NA'
					se_n_hit_genes = 'NA'
				elif n_observed_trait_runs == 1:
					mean_n_tested_genes = str(np.mean(tested_counts))
					se_n_tested_genes = '0.0'
					mean_n_hit_genes = str(np.mean(counts))
					se_n_hit_genes = '0.0'
				else:
					mean_n_tested_genes = str(np.mean(tested_counts))
					se_n_tested_genes = str(np.std(tested_counts, ddof=1)/np.sqrt(n_observed_trait_runs))
					mean_n_hit_genes = str(np.mean(counts))
					se_n_hit_genes = str(np.std(counts, ddof=1)/np.sqrt(n_observed_trait_runs))
				t.write(str(training_sample_size) + '\t' + z_score_type + '\t' + str(abs_z_threshold) + '\t0.05\t' + str(len(matched_runs)) + '\t' + str(n_observed_trait_runs) + '\t' + str(n_missing_trait_runs) + '\t' + mean_n_tested_genes + '\t' + se_n_tested_genes + '\t' + mean_n_hit_genes + '\t' + se_n_hit_genes + '\n')

	t.close()
	print(output_file)
	return


def run_matched_tissue_downsample_fdr_summary(twas_dir, matched_runs, training_sample_sizes, output_file):
	anno_method = 'borzoi_effect_sizes'
	distribution = 'ldscore_grid_squared'
	standardize_geno = 'False'
	pvalue_columns = [
		('standard_twas', 10),
		('hybrid_twas', 7)]

	t = open(output_file,'w')
	t.write('training_sample_size\ttwas_type\tcis_snp_h2_pvalue_threshold\tfdr_threshold\tn_expected_trait_runs\tn_observed_trait_runs\tn_missing_trait_runs\tmean_n_tested_genes_per_trait\tse_n_tested_genes_per_trait\tmean_n_fdr_significant_genes_per_trait\tse_n_fdr_significant_genes_per_trait\n')

	for training_sample_size in training_sample_sizes:
		n_tested_gene_counts = []
		fdr_hit_counts = {}
		for twas_type, column_index in pvalue_columns:
			fdr_hit_counts[twas_type] = []

		n_missing_trait_runs = 0
		for run_name in matched_runs:
			twas_result_file = twas_dir + 'twas_summary_' + run_name + '_' + anno_method + '_' + distribution + '_' + standardize_geno + '_downsample_' + str(training_sample_size) + '.txt'
			if os.path.exists(twas_result_file) == False:
				n_missing_trait_runs = n_missing_trait_runs + 1
				continue

			run_pvalues = {}
			for twas_type, column_index in pvalue_columns:
				run_pvalues[twas_type] = []

			f = open(twas_result_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				cis_snp_h2_pvalue = float(data[19])
				if cis_snp_h2_pvalue >= .05:
					continue
				for twas_type, column_index in pvalue_columns:
					run_pvalues[twas_type].append(float(data[column_index]))
			f.close()

			n_tested_gene_counts.append(len(run_pvalues[pvalue_columns[0][0]]))
			for twas_type, column_index in pvalue_columns:
				qvalues = benjamini_hochberg_qvalues(run_pvalues[twas_type])
				fdr_hit_counts[twas_type].append(np.sum(qvalues < .05))

		for twas_type, column_index in pvalue_columns:
			counts = np.asarray(fdr_hit_counts[twas_type])
			tested_counts = np.asarray(n_tested_gene_counts)
			n_observed_trait_runs = len(counts)
			if n_observed_trait_runs == 0:
				mean_n_tested_genes = 'NA'
				se_n_tested_genes = 'NA'
				mean_n_fdr_hits = 'NA'
				se_n_fdr_hits = 'NA'
			elif n_observed_trait_runs == 1:
				mean_n_tested_genes = str(np.mean(tested_counts))
				se_n_tested_genes = '0.0'
				mean_n_fdr_hits = str(np.mean(counts))
				se_n_fdr_hits = '0.0'
			else:
				mean_n_tested_genes = str(np.mean(tested_counts))
				se_n_tested_genes = str(np.std(tested_counts, ddof=1)/np.sqrt(n_observed_trait_runs))
				mean_n_fdr_hits = str(np.mean(counts))
				se_n_fdr_hits = str(np.std(counts, ddof=1)/np.sqrt(n_observed_trait_runs))
			t.write(str(training_sample_size) + '\t' + twas_type + '\t0.05\t0.05\t' + str(len(matched_runs)) + '\t' + str(n_observed_trait_runs) + '\t' + str(n_missing_trait_runs) + '\t' + mean_n_tested_genes + '\t' + se_n_tested_genes + '\t' + mean_n_fdr_hits + '\t' + se_n_fdr_hits + '\n')

	t.close()
	print(output_file)
	return


def count_twas_result_genes(twas_result_file):
	if os.path.exists(twas_result_file) == False:
		return 'NA'
	line_count = 0
	f = open(twas_result_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_count = line_count + 1
	f.close()
	return str(line_count)


def run_non_gtex_twas_file_loop(twas_dir, non_gtex_runs, non_gtex_twas_file_summary_file):
	anno_method = 'borzoi_effect_sizes'
	distribution = 'ldscore_grid_squared'
	standardize_geno = 'False'
	ref_tissue = 'Muscle_Skeletal'


	for run_info in non_gtex_runs:
		trait_name = run_info[0]
		non_gtex_tissue = run_info[1]
		non_gtex_sample = run_info[2]
		target_gtex_tissue = run_info[3]
		target_gtex_sample = run_info[4]

		non_gtex_result_file = twas_dir + 'twas_summary_' + trait_name + '_' + non_gtex_tissue + '_' + target_gtex_tissue + '_' + ref_tissue + '_' + anno_method + '_' + distribution + '_' + standardize_geno + '.txt'
		matched_proxy_result_file = twas_dir + 'twas_summary_' + trait_name + '_matched_proxy_for_' + non_gtex_tissue + '_' + target_gtex_tissue + '_' + ref_tissue + '_' + anno_method + '_' + distribution + '_' + standardize_geno + '.txt'

		for analysis_type, twas_result_file in (('non_gtex', non_gtex_result_file), ('matched_proxy', matched_proxy_result_file)):
			file_exists = os.path.exists(twas_result_file)
			n_genes = count_twas_result_genes(twas_result_file)


			if analysis_type != 'non_gtex':
				continue

			namer = twas_result_file.split('/twas_summary_')[1].split('_Muscle')[0]

			f = open(twas_result_file)
			head_count = 0
			for line in f:
				line =line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				twas_z = float(data[8])
				borzoi_twas_z = float(data[11])
				borzoi_prob = float(data[15])
				cis_h2_p = float(data[19])
				borzoi_fsr = float(data[17])
				hybrid_prob = float(data[14])

				gene_name = data[0]


				#if np.isnan(borzoi_prob) == False and np.min([borzoi_prob, 1.0 - borzoi_prob]) < .1:
				if borzoi_fsr < .1 and np.abs(borzoi_twas_z) > 3:
					print(namer + '\t' + gene_name + '\t' + str(twas_z) + '\t' + str(borzoi_twas_z) + '\t' + str(borzoi_prob) + '\t' + str(borzoi_fsr) + '\t' + str(cis_h2_p))

			f.close()



	return

def load_non_gtex_runs():
	non_gtex_run_lines = [
		'bmd_HEEL_TSCOREz\tosteoblast\tENCFF953CRN\tMuscle_Skeletal\tGTEX-13QJ3-0726-SM-5SI68.1',
		'bmd_HEEL_TSCOREz\tosteocyte\tENCFF092QIW\tMuscle_Skeletal\tGTEX-13QJ3-0726-SM-5SI68.1',
		'body_HEIGHTz\tknee_articular_chondrocyte\tENCFF314QIG\tMuscle_Skeletal\tGTEX-13QJ3-0726-SM-5SI68.1',
		'body_HEIGHTz\tchondrocyte\tENCFF958WEW\tMuscle_Skeletal\tGTEX-13QJ3-0726-SM-5SI68.1',
		'body_HEIGHTz\tfetal_skeletal_muscle\tENCFF285OHV\tMuscle_Skeletal\tGTEX-13QJ3-0726-SM-5SI68.1',
		'body_HEIGHTz\tskeletal_muscle_satellite_cell\tENCFF804YRV\tMuscle_Skeletal\tGTEX-13QJ3-0726-SM-5SI68.1',
		'disease_ASTHMA_DIAGNOSED\tbronchial_epithelial_cell\tENCFF153YEN\tLung\tGTEX-1399S-1726-SM-5L3DI.1',
		'disease_ASTHMA_DIAGNOSED\tbronchial_smooth_muscle_cell\tENCFF383NTW\tLung\tGTEX-1399S-1726-SM-5L3DI.1',
		'lung_FEV1FVCzSMOKE\talveolar_epithelial_cell\tENCFF734VZN\tLung\tGTEX-1399S-1726-SM-5L3DI.1',
		'lung_FVCzSMOKE\tairway_epithelial_cell\tENCFF595CZA\tLung\tGTEX-1399S-1726-SM-5L3DI.1',
		'disease_RESPIRATORY_ENT\tbronchus_fibroblast\tENCFF383RPD\tLung\tGTEX-1399S-1726-SM-5L3DI.1',
		'lung_FEV1FVCzSMOKE\tfetal_lung\tENCFF892OBT\tLung\tGTEX-1399S-1726-SM-5L3DI.1',
		'disease_T2D\tpancreatic_beta_cell\tENCFF995AUL\tPancreas\tGTEX-11I78-0626-SM-5A5LZ.1',
		'biochemistry_HbA1c\tpancreatic_beta_cell\tENCFF995AUL\tPancreas\tGTEX-11I78-0626-SM-5A5LZ.1',
		'biochemistry_Glucose\tendocrine_pancreas_progenitor_cell\tENCFF225CYB\tPancreas\tGTEX-11I78-0626-SM-5A5LZ.1',
		'disease_T2D\tendocrine_pancreas_progenitor_cell\tENCFF225CYB\tPancreas\tGTEX-11I78-0626-SM-5A5LZ.1',
		'bp_SYSTOLICadjMEDz\taortic_smooth_muscle_cell\tENCFF281BWX\tArtery_Aorta\tGTEX-1JK1U-0426-SM-CYPSP.1',
		'bp_DIASTOLICadjMEDz\taortic_smooth_muscle_cell\tENCFF281BWX\tArtery_Aorta\tGTEX-1JK1U-0426-SM-CYPSP.1',
		'disease_HYPERTENSION_DIAGNOSED\taortic_smooth_muscle_cell\tENCFF281BWX\tArtery_Aorta\tGTEX-1JK1U-0426-SM-CYPSP.1',
		'disease_CARDIOVASCULAR\tcoronary_artery_smooth_muscle_cell\tENCFF537AIY\tArtery_Aorta\tGTEX-1JK1U-0426-SM-CYPSP.1',
		'disease_CARDIOVASCULAR\tcoronary_artery_endothelial_cell\tENCFF226UWU\tArtery_Aorta\tGTEX-1JK1U-0426-SM-CYPSP.1',
		'bp_SYSTOLICadjMEDz\tfetal_metanephros\tENCFF367PUX\tKidney_Cortex\tGTEX-13112-2126-SM-5GCO4.1',
		'disease_AID_ALL\tregulatory_t_cell\tENCFF772RIQ\tWhole_Blood\tGTEX-1LB8K-0005-SM-DIPED.1',
		'disease_PSORIASIS\tactivated_cd4_t_cell\tENCFF362PNM\tWhole_Blood\tGTEX-1LB8K-0005-SM-DIPED.1',
		'disease_ALLERGY_ECZEMA_DIAGNOSED\tkeratinocyte\tENCFF703WGS\tSkin_Sun_Exposed_Lower_leg\tGTEX-13U4I-0126-SM-5LU38.1',
		'repro_MENARCHE_AGE\tfetal_diencephalon\tENCFF355IUO\tPituitary\tGTEX-12WSC-3126-SM-5GCNB.1',
		'cov_EDU_COLLEGE\tfetal_frontal_cortex\tENCFF217HQN\tBrain_Cortex\tGTEX-1H3O1-1726-SM-9WYSR.1',
		'mental_NEUROTICISM\tganglionic_eminence_neurosphere\tENCFF674RUW\tBrain_Cortex\tGTEX-1H3O1-1726-SM-9WYSR.1',
		'pigment_TANNING\tadult_skin_melanocyte\tENCFF230ZSG\tSkin_Sun_Exposed_Lower_leg\tGTEX-13U4I-0126-SM-5LU38.1',
		'body_BALDING1\thair_follicle_dermal_papilla_cell\tENCFF798NEI\tSkin_Sun_Exposed_Lower_leg\tGTEX-13U4I-0126-SM-5LU38.1',
		'bp_SYSTOLICadjMEDz\tglomerular_endothelial_cell\tENCFF863JIL\tKidney_Cortex\tGTEX-13112-2126-SM-5GCO4.1',
		'disease_HYPERTENSION_DIAGNOSED\tmesangial_cell\tENCFF346QDJ\tKidney_Cortex\tGTEX-13112-2126-SM-5GCO4.1',
		'bp_DIASTOLICadjMEDz\tproximal_tubule_epithelial_cell\tENCFF767MLU\tKidney_Cortex\tGTEX-13112-2126-SM-5GCO4.1',
		'repro_NumberChildrenEverBorn_Pooled\tplacental_epithelial_cell\tENCFF591NKN\tUterus\tGTEX-13FTX-1026-SM-5J2O5.1',
		'repro_NumberChildrenEverBorn_Pooled\tplacental_pericyte\tENCFF779FMR\tUterus\tGTEX-13FTX-1026-SM-5J2O5.1',
		'repro_NumberChildrenEverBorn_Pooled\tvillous_mesenchyme_fibroblast\tENCFF649CUP\tUterus\tGTEX-13FTX-1026-SM-5J2O5.1',
		'cov_EDU_COLLEGE\tneural_progenitor_cell\tENCFF403DKN\tBrain_Cortex\tGTEX-1H3O1-1726-SM-9WYSR.1',
		'mental_NEUROTICISM\tcortex_neurosphere\tENCFF615IJV\tBrain_Cortex\tGTEX-1H3O1-1726-SM-9WYSR.1',
		'other_MORNINGPERSON\tfetal_diencephalon\tENCFF355IUO\tPituitary\tGTEX-12WSC-3126-SM-5GCNB.1',
		'body_HEIGHTz\tneural_crest_cell\tENCFF024MEB\tBrain_Cortex\tGTEX-1H3O1-1726-SM-9WYSR.1',
		'disease_ALLERGY_ECZEMA_DIAGNOSED\tdermal_fibroblast\tENCFF343FWG\tSkin_Sun_Exposed_Lower_leg\tGTEX-13U4I-0126-SM-5LU38.1',
		'body_BALDING1\thair_follicular_keratinocyte\tENCFF544MJM\tSkin_Sun_Exposed_Lower_leg\tGTEX-13U4I-0126-SM-5LU38.1',
		'disease_CARDIOVASCULAR\taortic_adventitial_fibroblast\tENCFF256ZZS\tArtery_Aorta\tGTEX-1JK1U-0426-SM-CYPSP.1',
		'disease_CARDIOVASCULAR\tcardiac_ventricle_fibroblast\tENCFF287LRZ\tHeart_Left_Ventricle\tGTEX-18465-0926-SM-731AY.1',
		'disease_AID_ALL\tactivated_t_cell_il2_cd3_cd28\tENCFF122XJD\tWhole_Blood\tGTEX-1LB8K-0005-SM-DIPED.1',
		'disease_PSORIASIS\tactivated_cd4_memory_t_cell\tENCFF330EUN\tWhole_Blood\tGTEX-1LB8K-0005-SM-DIPED.1']

	non_gtex_runs = []
	for line in non_gtex_run_lines:
		non_gtex_runs.append(line.split('\t'))
	return np.asarray(non_gtex_runs)









######################
# Command line args
######################
twas_dir = sys.argv[1]


#########
# 0. loop through non-GTEx TWAS files and matched proxy TWAS files
non_gtex_runs = load_non_gtex_runs()
non_gtex_twas_file_summary_file = twas_dir + 'non_gtex_twas_file_summary.txt'
run_non_gtex_twas_file_loop(twas_dir, non_gtex_runs, non_gtex_twas_file_summary_file)


#########
# 1. run analysis where we have matched tissues and look at FSR
matched_runs = []
matched_runs.append('disease_AID_ALL_Spleen_GTEX-14PJ4-0526-SM-6871G.1')
matched_runs.append('blood_MONOCYTE_COUNT_Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1')
matched_runs.append('blood_MEAN_CORPUSCULAR_HEMOGLOBIN_Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1')
matched_runs.append('blood_MEAN_PLATELET_VOL_Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1')
matched_runs.append('blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT_Cells_Cultured_fibroblasts_GTEX-139TS-0008-SM-62LDG.1')
matched_runs.append('disease_ALLERGY_ECZEMA_DIAGNOSED_Skin_Sun_Exposed_Lower_leg_GTEX-13U4I-0126-SM-5LU38.1')
matched_runs.append('biochemistry_VitaminD_Skin_Sun_Exposed_Lower_leg_GTEX-13U4I-0126-SM-5LU38.1')
matched_runs.append('biochemistry_Cholesterol_Liver_GTEX-11EQ9-0526-SM-5A5JZ.1')
matched_runs.append('repro_MENARCHE_AGE_Pituitary_GTEX-12WSC-3126-SM-5GCNB.1')
matched_runs.append('bp_DIASTOLICadjMEDz_Artery_Aorta_GTEX-1JK1U-0426-SM-CYPSP.1')
matched_runs.append('lung_FVCzSMOKE_Lung_GTEX-1399S-1726-SM-5L3DI.1')
matched_runs.append('lung_FEV1FVCzSMOKE_Lung_GTEX-1399S-1726-SM-5L3DI.1')
matched_runs.append('body_HEIGHTz_Cells_Cultured_fibroblasts_GTEX-139TS-0008-SM-62LDG.1')
matched_runs.append('bmd_HEEL_TSCOREz_Cells_Cultured_fibroblasts_GTEX-139TS-0008-SM-62LDG.1')

matched_runs = np.asarray(matched_runs)

matched_tissue_fsr_summary_file = twas_dir + 'matched_tissue_fsr_summary.txt'
run_matched_tissue_fsr_analysis(twas_dir, matched_runs, matched_tissue_fsr_summary_file)

matched_tissue_borzoi_fsr_filtered_fdr_summary_file = twas_dir + 'matched_tissue_borzoi_fsr_filtered_fdr_summary.txt'
run_matched_tissue_borzoi_fsr_filtered_fdr_analysis(twas_dir, matched_runs, matched_tissue_borzoi_fsr_filtered_fdr_summary_file)

matched_tissue_downsample_hit_summary_file = twas_dir + 'matched_tissue_downsample_hit_summary.txt'
training_sample_sizes = np.asarray([25, 50, 100, 150, 200, 250])
run_matched_tissue_downsample_hit_summary(twas_dir, matched_runs, training_sample_sizes, matched_tissue_downsample_hit_summary_file)

matched_tissue_downsample_fdr_summary_file = twas_dir + 'matched_tissue_downsample_fdr_summary.txt'
run_matched_tissue_downsample_fdr_summary(twas_dir, matched_runs, training_sample_sizes, matched_tissue_downsample_fdr_summary_file)
