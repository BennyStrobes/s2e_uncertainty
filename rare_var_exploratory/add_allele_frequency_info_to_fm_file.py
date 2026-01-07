import numpy as np
import os 
import sys
import pdb





def create_mapping_from_variant_id_to_af(hg19_1000G_genotype_dir):
	dicti = {}

	for chrom_num in range(1,23):
		# First create mapping from rsid to variant
		rsid_to_variant_id = {}
		bim_file = hg19_1000G_genotype_dir + '1000G.EUR.QC.' + str(chrom_num) + '.bim'
		f = open(bim_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			rsid = data[1]
			snp_id = data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5]
			snp_id_alt= data[0] + '_' + data[3] + '_' + data[5] + '_' + data[4]
			if rsid in rsid_to_variant_id:
				print('asssumpiotneornoer')
				pdb.set_trace()
			rsid_to_variant_id[rsid] = (snp_id, snp_id_alt)
		f.close()

		# Now add mapping from variant to af
		af_file = hg19_1000G_genotype_dir + '1000G.EUR.QC.' + str(chrom_num) + '.frq'
		f = open(af_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if len(data) != 6:
				print('assumption eroror')
				pdb.set_trace()
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[1]
			af = data[4]
			if float(af) > .5:
				print('asumpitoneornor')
				pdb.set_trace()
			(snp_id, snp_id_alt) = rsid_to_variant_id[rsid]

			dicti[snp_id] = af
			dicti[snp_id_alt] = af
		f.close()


	return dicti








#########################
# Command line args
#########################
fm_eqtl_sumstats_file = sys.argv[1]
fm_and_af_eqtl_sumstats_file = sys.argv[2]
hg19_1000G_genotype_dir = sys.argv[3]


# First create mapping from variant id to allele frequency
variant_id_to_af = create_mapping_from_variant_id_to_af(hg19_1000G_genotype_dir)


f = open(fm_eqtl_sumstats_file)
t = open(fm_and_af_eqtl_sumstats_file,'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'maf\n')
		continue
	variant_id = data[0] + '_' + data[1] + '_' + data[3] + '_' + data[4]

	if variant_id in variant_id_to_af:
		maf = variant_id_to_af[variant_id]
	else:
		maf = 'NA'
	t.write(line + '\t' + maf + '\n')
f.close()
t.close()
print(fm_and_af_eqtl_sumstats_file)
