import numpy as np
import os
import sys
import pdb











#####################
# Command line args
#####################
sim_results_dir = sys.argv[1]
organized_sim_results_dir = sys.argv[2]


sample_sizes = ['50000', '100000']
n_detecteds = ['1', '2', '3', '4', '5']
sample_sizes = ['50000']
n_detecteds = ['0', '1', '2', '3', '4', '5']

t1e_output_file = organized_sim_results_dir + 'organized_t1e_sim_results.txt'
power_output_file = organized_sim_results_dir + 'organized_power_sim_results.txt'

t_t1e = open(t1e_output_file,'w')
t_power = open(power_output_file,'w')

t_t1e.write('method\tsample_size\tn_detected\tt1e\tt1e_lb\tt1e_ub\n')
t_power.write('method\tsample_size\tn_detected\tpower\tpower_lb\tpower_ub\n')


for sample_size in sample_sizes:
	for n_detected in n_detecteds:
		pvals = []
		labels = []
		methods = []
		for sim_iter in range(1,101):
			sim_results_file = sim_results_dir + 'sim_' + str(sim_iter) + '_' + sample_size + '_' + n_detected + '_sim_summary.txt'
			f = open(sim_results_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				methods.append(data[1])
				if data[2] == '0.0':
					labels.append(0.0)
				else:
					labels.append(1.0)
				pvals.append(float(data[-1]))
			f.close()
		labels = np.asarray(labels)
		pvals = np.asarray(pvals)
		methods = np.asarray(methods)

		# T1E
		indices = (labels == 0.0) & (methods == 'burden')
		t1e = np.sum(pvals[indices] < 0.05)/len(pvals[indices])
		t1e_se = np.sqrt(t1e*(1.0-t1e)/len(pvals[indices]))
		t1e_lb = t1e - (1.96*t1e_se)
		t1e_ub = t1e + (1.96*t1e_se)
		t_t1e.write('burden\t' + sample_size + '\t' + n_detected + '\t' + str(t1e) + '\t' + str(t1e_lb) + '\t' + str(t1e_ub) + '\n')

		indices = (labels == 0.0) & (methods == 'dRVAT')
		t1e = np.sum(pvals[indices] < 0.05)/len(pvals[indices])
		t1e_se = np.sqrt(t1e*(1.0-t1e)/len(pvals[indices]))
		t1e_lb = t1e - (1.96*t1e_se)
		t1e_ub = t1e + (1.96*t1e_se)
		t_t1e.write('dRVAT\t' + sample_size + '\t' + n_detected + '\t' + str(t1e) + '\t' + str(t1e_lb) + '\t' + str(t1e_ub) + '\n')

		# power
		indices = (labels == 1.0) & (methods == 'burden')
		t1e = np.sum(pvals[indices] < 0.05)/len(pvals[indices])
		t1e_se = np.sqrt(t1e*(1.0-t1e)/len(pvals[indices]))
		t1e_lb = t1e - (1.96*t1e_se)
		t1e_ub = t1e + (1.96*t1e_se)
		t_power.write('burden\t' + sample_size + '\t' + n_detected + '\t' + str(t1e) + '\t' + str(t1e_lb) + '\t' + str(t1e_ub) + '\n')

		indices = (labels == 1.0) & (methods == 'dRVAT')
		t1e = np.sum(pvals[indices] < 0.05)/len(pvals[indices])
		t1e_se = np.sqrt(t1e*(1.0-t1e)/len(pvals[indices]))
		t1e_lb = t1e - (1.96*t1e_se)
		t1e_ub = t1e + (1.96*t1e_se)
		t_power.write('dRVAT\t' + sample_size + '\t' + n_detected + '\t' + str(t1e) + '\t' + str(t1e_lb) + '\t' + str(t1e_ub) + '\n')


t_t1e.close()
t_power.close()

