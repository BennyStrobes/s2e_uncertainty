import numpy as np
import os
import sys
import pdb






######################
# Command line args
######################
fm_summary_file = sys.argv[1]
hg19_variant_pos_file = sys.argv[2]


t = open(hg19_variant_pos_file,'w')
f = open(fm_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	var_id = 'chr' + data[0] + ' ' + data[1] + ' ' +  str(int(data[1]) + 1)
	t.write(var_id + '\n')

t.close()
f.close()
