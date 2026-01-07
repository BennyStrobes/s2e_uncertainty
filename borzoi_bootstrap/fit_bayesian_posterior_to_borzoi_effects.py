import numpy as np
import os
import sys
import pdb
from gibbs_ash import model



def load_in_bootstrapped_effect_sizes_plus_standard_errors(organized_borzoi_bs_effecs_file):
	f = open(organized_borzoi_bs_effecs_file)
	effect_sizes = []
	ses = []
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		bs_vals = np.asarray(data[5].split(';')).astype(float)
		mean = np.mean(bs_vals)
		se = np.std(bs_vals, ddof=1)

		effect_sizes.append(mean)
		ses.append(se)
		counter = counter + 1
		if counter == 100000:
			break

	f.close()

	return np.asarray(effect_sizes), np.asarray(ses)








#####################
# Command line args
#####################
organized_borzoi_bs_effecs_file = sys.argv[1]
bayes_fit_output_root = sys.argv[2]


# Load in vector of effect sizes and standard errors
effect_sizes, ses = load_in_bootstrapped_effect_sizes_plus_standard_errors(organized_borzoi_bs_effecs_file)

mm = model(max_iter=500, burn_in_iter=400)
mm.fit(effect_sizes, ses)

pdb.set_trace()





