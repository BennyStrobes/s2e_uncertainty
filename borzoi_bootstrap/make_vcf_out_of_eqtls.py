import numpy as np
import os
import sys
import pdb



input_file = sys.argv[1]
output_file = sys.argv[2]



f = open(input_file)
t = open(output_file,'w')
t.write('##fileformat=VCFv4.2\n')
used = {}
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		header = np.copy(np.asarray(data))
		continue

	if data[3] not in used:
		t.write(data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5] + '\n')

		used[data[3]] = 1


f.close()
t.close()