import numpy as np
import os
import sys
import pdb





def concatenate_results(file_stem, borzoi_threshs):
	output_file = file_stem + '.txt'
	t = open(output_file,'w')

	itera = 0

	for borzoi_thresh in borzoi_threshs:
		for chrom_num in range(1,23):
			filer = file_stem + '_chr' + str(chrom_num) + '_thresh_' + borzoi_thresh + '.txt'
			if chrom_num == 4:
				continue
			f = open(filer)
			head_count = 0
			for line in f:
				line = line.rstrip()
				if head_count == 0:
					head_count = head_count + 1
					if itera == 0:
						t.write(line + '\t' 'borzoi_thresh' + '\n')
					continue
				t.write(line + '\t' + borzoi_thresh + '\n')
			itera = itera + 1
			f.close()
	t.close()

	return




#########################
# COmmand line args
########################
file_stem = sys.argv[1]


borzoi_threshs = ['0.05', '0.075', '0.1', '0.15']
concatenate_results(file_stem, borzoi_threshs)
output_file = file_stem + '.txt'



aa = np.loadtxt(output_file,dtype=str,delimiter='\t')
bb = aa[1:,:]
thresholds = bb[:,-1]
noncoding_effect = bb[:,-3].astype(float)
coding_effect = bb[:,5].astype(float)

for borzoi_thresh in borzoi_threshs:

	indices = thresholds == borzoi_thresh

	n_agree = np.sum(np.sign(coding_effect[indices]) != np.sign(noncoding_effect[indices]))
	total = len(coding_effect[indices])

	print(borzoi_thresh + '\t' + str(n_agree/total) + '\t' + str(n_agree) + '\t' + str(total))











