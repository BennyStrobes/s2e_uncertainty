import numpy as np
import sys
import pdb
import scipy.stats
import os



def get_same_sign_indices(borzoi_preds):
	nrows = borzoi_preds.shape[0]
	arr = []
	for row_iter in range(nrows):
		sorted_row_preds = np.sort(borzoi_preds[row_iter, :])
		miny = sorted_row_preds[0]
		maxy = sorted_row_preds[-1]
		kk = np.sum(np.sign(sorted_row_preds) == 1.0)
		nn = len(sorted_row_preds)
		result = scipy.stats.binomtest(kk, nn, p=0.5, alternative='two-sided')
		if result.pvalue < .01:
			print(result.pvalue)
			arr.append(True)
		else:
			arr.append(False)

	return np.asarray(arr)


def compute_sign_test_pvalue(borzoi_preds):
	kk = np.sum(np.sign(borzoi_preds) == 1.0)
	nn = len(borzoi_preds)
	result = scipy.stats.binomtest(kk, nn, p=0.5, alternative='two-sided')
	return result.pvalue


######################
# Command line args
######################
bs_iter_start_inclusive = int(sys.argv[1])
bs_iter_end_inclusive = int(sys.argv[2])
results_dir = sys.argv[3]
organized_bs_eqtl_output = sys.argv[4]


# Open input files
tissues = []
snps = []
genes = []
borzoi_preds = []
for bs_iter in range(bs_iter_start_inclusive, bs_iter_end_inclusive + 1):

	results_file = results_dir + 'bs' + str(bs_iter) + '_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt'
	raw_data = np.loadtxt(results_file, dtype=str, delimiter='\t')

	borzoi_preds.append(raw_data[1:,-1].astype(float))
	tissues.append(raw_data[1:,0])
	snps.append(raw_data[1:,3])
	genes.append(raw_data[1:,7])

tissues = np.transpose(np.asarray(tissues))
snps = np.transpose(np.asarray(snps))
genes = np.transpose(np.asarray(genes))
borzoi_preds = np.transpose(np.asarray(borzoi_preds))

ref_bs_iter = bs_iter_start_inclusive + 0
orig_results_file = results_dir + 'bs' + str(bs_iter) + '_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt'
f = open(orig_results_file)
t = open(organized_bs_eqtl_output,'w')
head_count = 0
line_counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write('\t'.join(data[:-1]) + '\t' + 'mean_borzoi_log_sed\tsdev_borzoi_log_sed\tbs_borzoi_log_sed\n')
		continue

	if data[0] != tissues[line_counter,0] or data[3] != snps[line_counter,0] or data[7] != genes[line_counter,0]:
		print('assumption erororo')
		pdb.set_trace()

	if len(np.unique(tissues[line_counter,:])) != 1 or len(np.unique(snps[line_counter,:])) != 1 or len(np.unique(genes[line_counter,:])) != 1:
		print('asssumption eroror')
		pdb.set_trace()

	t.write('\t'.join(data[:-1]) + '\t')

	preds = np.copy(borzoi_preds[line_counter,:])
	#preds = np.random.choice(borzoi_preds[line_counter,:], 10, replace=False)

	meany = np.mean(preds)
	sdev = np.std(preds)
	#sign_test_pvalue = compute_sign_test_pvalue(preds)
	vals = ';'.join(preds.astype(str))

	t.write(str(meany) + '\t' + str(sdev) + '\t' + vals + '\n')

	line_counter = line_counter + 1
f.close()
t.close()

aa = np.loadtxt(organized_bs_eqtl_output,dtype=str,delimiter='\t')

susie = aa[1:,-5].astype(float)
borzoi_mean = aa[1:,-3].astype(float)
borzoi_sd = aa[1:,-2].astype(float)
zeds = borzoi_mean/borzoi_sd

print(scipy.stats.pearsonr(susie,borzoi_mean))


sig_indices = np.abs(zeds) > 1.0
print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))
sig_indices = np.abs(zeds) > 2.0
print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))
sig_indices = np.abs(zeds) > 3.0
#print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))
#print(np.sum(sig_indices))



'''
aa = np.loadtxt(organized_bs_eqtl_output,dtype=str,delimiter='\t')

susie = aa[1:,-5].astype(float)
borzoi_mean = aa[1:,-3].astype(float)
borzoi_sd = aa[1:,-2].astype(float)
zeds = borzoi_mean/borzoi_sd


print(scipy.stats.pearsonr(susie,borzoi_mean))


sig_indices = np.abs(zeds) > 1.0
print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))
sig_indices = np.abs(zeds) > 2.0
print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))
sig_indices = np.abs(zeds) > 3.0
print(scipy.stats.pearsonr(susie[sig_indices],borzoi_mean[sig_indices]))

same_sign_indices = get_same_sign_indices(borzoi_preds)
print(scipy.stats.pearsonr(susie[same_sign_indices],borzoi_mean[same_sign_indices]))
'''