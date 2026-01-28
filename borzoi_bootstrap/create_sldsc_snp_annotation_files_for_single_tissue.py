import numpy as np
import os
import sys
import pdb
import gzip



def create_mapping_from_borzoi_sig_to_variant_id(borzoi_variant_effect_file):
	f = open(borzoi_variant_effect_file)
	dicti = {}
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue


		bs_vals = np.asarray(data[5].split(';')).astype(float)

		borzoi_effect = np.mean(bs_vals)
		np_p_plus = (1 + np.sum(bs_vals >= 0))/(1+len(bs_vals))
		np_p_minus = (1 + np.sum(bs_vals <=0))/(1+len(bs_vals))
		borzoi_pvalue = 2.0*np.min((np_p_plus, np_p_minus))



		variant_id = data[2]

		if variant_id not in dicti:
			dicti[variant_id] = borzoi_pvalue
		else:
			dicti[variant_id] = np.min([borzoi_pvalue, dicti[variant_id]])
		counter = counter +1

	f.close()

	return dicti

def create_mapping_from_rsid_to_variant_id(genotype_dir):
	dicti = {}
	for chrom_num in range(1,23):
		chrom_bim_file = genotype_dir + str(chrom_num) + '.bim'
		f = open(chrom_bim_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			rsid = data[1]
			variant_id1 = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5] + '_b38'
			variant_id2 = 'chr' + data[0] + '_' + data[3] + '_' + data[5] + '_' + data[4] + '_b38'

			if rsid in dicti:
				print('assumption eororr')
				pdb.set_trace()

			dicti[rsid] = (variant_id1, variant_id2)


		f.close()
	return dicti



################
# Command line args
tissue_name = sys.argv[1]
sldsc_processed_anno_dir = sys.argv[2]
borzoi_res_file = sys.argv[3]
ldsc_snp_annotation_dir = sys.argv[4]
genotype_1000G_plink_stem = sys.argv[5]



# First extract mapping from variants to borzoi significnace (in at least one gene)
variant_id_to_borzoi_sig = create_mapping_from_borzoi_sig_to_variant_id(borzoi_res_file)

# Now create mapping from rsid to variant id
rsid_to_variantid = create_mapping_from_rsid_to_variant_id(genotype_1000G_plink_stem)



# Loop through chromosomes
for chrom_num in range(1,23):
	print(chrom_num)
	chrom_string = str(chrom_num)

	# Open output file handles
	t = {}
	t['borzoi_sig0'] = gzip.open(sldsc_processed_anno_dir + tissue_name + '_borzoi_sig0.' + chrom_string + '.annot.gz','wt')
	t['borzoi_sig1'] = gzip.open(sldsc_processed_anno_dir + tissue_name + '_borzoi_sig1.' + chrom_string + '.annot.gz','wt')
	t['borzoi_sig2'] = gzip.open(sldsc_processed_anno_dir + tissue_name + '_borzoi_sig2.' + chrom_string + '.annot.gz','wt')
	t['borzoi_sig3'] = gzip.open(sldsc_processed_anno_dir + tissue_name + '_borzoi_sig3.' + chrom_string + '.annot.gz','wt')



	ref_file = ldsc_snp_annotation_dir + 'baselineLD.' + chrom_string + '.annot.gz'
	f = gzip.open(ref_file,'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t['borzoi_sig0'].write('\t'.join(np.asarray(data)[:4]) + '\t' + tissue_name + '_borzoi_sig0' + '\n')
			t['borzoi_sig1'].write('\t'.join(np.asarray(data)[:4]) + '\t' + tissue_name + '_borzoi_sig1' + '\n')
			t['borzoi_sig2'].write('\t'.join(np.asarray(data)[:4]) + '\t' + tissue_name + '_borzoi_sig2' + '\n')
			t['borzoi_sig3'].write('\t'.join(np.asarray(data)[:4]) + '\t' + tissue_name + '_borzoi_sig3' + '\n')
			continue

		rsid = data[2]

		# check if variant is tested
		tested_variant = False
		sig_val = None
		if rsid in rsid_to_variantid:
			variant_ids = rsid_to_variantid[rsid]
			if variant_ids[0] in variant_id_to_borzoi_sig:
				sig_val = variant_id_to_borzoi_sig[variant_ids[0]]
				tested_variant = True
			elif variant_ids[1] in variant_id_to_borzoi_sig:
				sig_val = variant_id_to_borzoi_sig[variant_ids[1]]
				tested_variant = True
		
		row_lead = np.asarray(data)[:4]

		if tested_variant:
			if sig_val <= 0.05:
				t['borzoi_sig0'].write('\t'.join(row_lead) + '\t' + '1.0' + '\n')
				t['borzoi_sig1'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig2'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig3'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
			elif sig_val > 0.05 and sig_val <= 0.2:
				t['borzoi_sig0'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig1'].write('\t'.join(row_lead) + '\t' + '1.0' + '\n')
				t['borzoi_sig2'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig3'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
			elif sig_val > 0.2 and sig_val <= 0.5:
				t['borzoi_sig0'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig1'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig2'].write('\t'.join(row_lead) + '\t' + '1.0' + '\n')
				t['borzoi_sig3'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
			elif sig_val > 0.5:
				t['borzoi_sig0'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig1'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig2'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
				t['borzoi_sig3'].write('\t'.join(row_lead) + '\t' + '1.0' + '\n')
			else:
				print('assumpitoneororneonre')
				pdb.set_trace()
		else:
			t['borzoi_sig0'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
			t['borzoi_sig1'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
			t['borzoi_sig2'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')
			t['borzoi_sig3'].write('\t'.join(row_lead) + '\t' + '0.0' + '\n')



	f.close()
	t['borzoi_sig0'].close()
	t['borzoi_sig1'].close()
	t['borzoi_sig2'].close()
	t['borzoi_sig3'].close()


