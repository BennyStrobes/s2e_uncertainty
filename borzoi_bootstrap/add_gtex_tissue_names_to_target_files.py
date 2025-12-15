import numpy as np
import os
import sys
import pdb
import re



def to_underscore(s: str) -> str:
	s = s.replace(")", "")
	# replace " (" OR any run of spaces and hyphens with "_"
	return re.sub(r"(?:\s*\(|[ -]+)", "_", s)



################################
# Command line args
gtex_target_file = sys.argv[1]
full_gtex_target_file = sys.argv[2]
gtex_sample_attributes_file = sys.argv[3]


#####################
# First create mapping from gtex sample name to tissue name
#####################
gtex_sample_to_tissue_name = {}
gtex_tissues = {}
f = open(gtex_sample_attributes_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gtex_sample_id = data[0]
	tissue_type = to_underscore(data[6])
	if tissue_type == 'Cells_EBV_transformed_lymphocytes':
		tissue_type = 'Cells_EBV-transformed_lymphocytes'
	elif tissue_type == 'Brain_Spinal_cord_cervical_c_1':
		tissue_type = 'Brain_Spinal_cord_cervical_c-1'
	gtex_sample_to_tissue_name[gtex_sample_id] = tissue_type
	gtex_tissues[tissue_type] = 1
f.close()


head_count = 0
f = open(gtex_target_file)
t = open(full_gtex_target_file,'w')
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write('\t'.join(data) + '\t' + 'tissue\n')
		continue
	short_gtex_sample_name = data[1].split('.')[0]
	if short_gtex_sample_name not in gtex_sample_to_tissue_name:
		print('assumption eroror')
		pdb.set_trace()
	tissue_name = gtex_sample_to_tissue_name[short_gtex_sample_name]
	t.write('\t'.join(data) + '\t' + tissue_name + '\n')
f.close()
t.close()

print(full_gtex_target_file)



