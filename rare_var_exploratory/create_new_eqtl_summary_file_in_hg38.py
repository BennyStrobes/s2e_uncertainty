import numpy as np
import os
import sys
import pdb
from pyfaidx import Fasta
import pysam







######################
# Command line args
######################
fm_and_af_eqtl_sumstats_file = sys.argv[1]
liftover_result = sys.argv[2]
fm_and_af_eqtl_sumstats_hg38_file = sys.argv[3]
sardinia_hg38_vcf = sys.argv[4]
fasta_file = sys.argv[5]

# First create mapping from hg19 to hg38
hg19_to_hg38 = {}
f = open(liftover_result)
for line in f:
	line = line.rstrip()
	data = line.split()
	hg19_var_info = data[3].split(':')
	hg19_chrom = hg19_var_info[0]
	hg19_pos = str(int(hg19_var_info[1].split('-')[0]) - 1)
	hg38_chrom = data[0]
	hg38_pos = data[1]
	if hg19_chrom != hg38_chrom:
		hg19_to_hg38[hg19_chrom + ':'+ hg19_pos] = (False, 'na', 'na')
	else:
		hg19_to_hg38[hg19_chrom + ':'+ hg19_pos] = (True, hg38_chrom, hg38_pos)
f.close()

#fasta = Fasta(fasta_file)
genome_open = pysam.Fastafile(fasta_file)
# convert files
f = open(fm_and_af_eqtl_sumstats_file)
t = open(fm_and_af_eqtl_sumstats_hg38_file,'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	hg19_chrom = 'chr' + data[0]
	hg19_pos = data[1]
	hg19_var = hg19_chrom + ':' + hg19_pos

	hg38_info = hg19_to_hg38[hg19_var]
	if hg38_info[0] == False:
		continue

	#ref = fasta[hg38_info[1]][int(hg38_info[2])+1]

	ref = genome_open.fetch(hg38_info[1], int(hg38_info[2])-1, int(hg38_info[2])).upper()

	if data[3] == ref:
		var_id = 'chr' + data[0] + '_' + hg38_info[2] + '_' +data[3] + '_' + data[4] + '_b38'
		t.write(data[0] + '\t' + hg38_info[2] + '\t' + var_id + '\t' + '\t'.join(data[3:]) + '\n')
	elif data[4] == ref:
		var_id = 'chr' + data[0] + '_' + hg38_info[2] + '_' +data[4] + '_' + data[3] + '_b38'
		t.write(data[0] + '\t' + hg38_info[2] + '\t' + var_id + '\t' + data[4] + '\t' + data[3] + '\t' + data[5] + '\t' + data[6] + '\t' + data[7] + '\t' + str(-float(data[8])) + '\t' + data[9] + '\t' + data[10] + '\t' + data[11] + '\n')
	else:
		print('ref allelel assumption oeroorr')
		continue
f.close()
t.close()
genome_open.close()


f = open(fm_and_af_eqtl_sumstats_hg38_file)
t = open(sardinia_hg38_vcf,'w')
used = {}
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	var_id = data[2]
	if var_id in used:
		continue
	used[var_id] = 1
	var_info = var_id.split('_')

	t.write(var_info[0] + '\t' + var_info[1] + '\t' + var_id + '\t' + var_info[2] + '\t' + var_info[3] + '\n')

f.close()
t.close()


print(sardinia_hg38_vcf)





