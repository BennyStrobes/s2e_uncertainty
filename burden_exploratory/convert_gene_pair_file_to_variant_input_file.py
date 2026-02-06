import numpy as np
import os
import sys
import pdb










######################
# Commmand line args
######################
vg_file = sys.argv[1] # input file
variant_vcf_file = sys.argv[2] # output file


f = open(vg_file)
t = open(variant_vcf_file,'w')
head_count = 0
used = {}
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count== 0:
		head_count = head_count + 1
		continue
	variant_id = data[1]
	# dont repeat variants
	if variant_id in used:
		continue
	used[variant_id] = 1
	t.write(data[0] + '\t' + data[2] + '\t' + data[1] + '\t' + data[3] + '\t' + data[4] + '\n')

	# Quick error check
	if data[1].split('_')[2] != data[3]:
		print('assumption eornror')
		pdb.set_trace()
	if data[1].split('_')[3] != data[4]:
		print('assumption eornror')
		pdb.set_trace()
f.close()
t.close()