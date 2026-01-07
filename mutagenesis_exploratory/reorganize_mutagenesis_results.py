import numpy as np
import os
import sys
import pdb











#####################
# Command line args
#####################
raw_mutagenesis_file = sys.argv[1]
reorganized_mutagenesis_file = sys.argv[2]
mutagenesis_variant_vcf = sys.argv[3]

sig_threshold = .001

valid_chromosomes = {}
for chrom_num in range(1,23):
	valid_chromosomes[str(chrom_num)] = 1

gene_to_ensembl = {
	"BCL11A": "ENSG00000119866",
	"FOXE1": "ENSG00000178919",
	"GP1BA": "ENSG00000185245",
	"HBB": "ENSG00000244734",
	"HBG1": "ENSG00000213934",
	"HNF4A": "ENSG00000140551",
	"IRF4": "ENSG00000137265",
	"IRF6": "ENSG00000185559",
	"LDLR": "ENSG00000130164",
	"MSMB": "ENSG00000096696",
	"MYC": "ENSG00000136997",
	"PKLR": "ENSG00000104055",
	"RET": "ENSG00000165731",
	"SORT1": "ENSG00000134323",
	"TCF7L2": "ENSG00000148737",
	"TERT": "ENSG00000164362",
	"ZFAND3": "ENSG00000125266",
	"LDLR.2": "ENSG00000130164",
	"MYCrs11986220": "ENSG00000136997",
	"MYCrs6983267": "ENSG00000136997"
}


f = open(raw_mutagenesis_file)
t = open(reorganized_mutagenesis_file,'w')
used_genes = {}
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'variant_id' + '\t' + 'ensamble_id' + '\n')
		continue
	line_chrom = data[0]
	if line_chrom not in valid_chromosomes:
		continue
	alt_allele = data[3]
	if alt_allele == '-':
		continue
	used_genes[data[-1]] = 1
	pvalue = float(data[-2])
	if pvalue > sig_threshold:
		continue
	if data[-1] not in gene_to_ensembl:
		continue
	variant_id = 'chr' + line_chrom + '_' + data[1] + '_' + data[2] + '_' + data[3] + '_b38'
	ensamble_id = gene_to_ensembl[data[-1]]
	t.write(line + '\t' + variant_id + '\t' + ensamble_id + '\n')

f.close()
t.close()

head_count = 0
f = open(reorganized_mutagenesis_file)
t = open(mutagenesis_variant_vcf,'w')
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	variant_id = data[10]
	var_info = variant_id.split('_')
	t.write(var_info[0] + '\t' + var_info[1] + '\t' + variant_id + '\t' + var_info[2] + '\t' + var_info[3] + '\n')
f.close()
t.close()