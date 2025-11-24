import json
import os
import time
import warnings
import sys
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import pysam
import pyfaidx
import pybedtools
import csv
import tensorflow as tf

from baskerville import seqnn
from baskerville import gene as bgene
from baskerville import dna

from borzoi_helpers import *

import pdb
'''
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
'''


def get_test_range(x, num_jobs, T):
	"""
	Determine the range of test indices assigned to job x (0-based).

	Args:
		x (int): Job index (0-based, from 0 to num_jobs - 1).
		num_jobs (int): Total number of parallel jobs.
		T (int): Total number of tests.

	Returns:
		(int, int): Tuple of (first_test_index, last_test_index), both 0-based.
	"""
	assert 0 <= x < num_jobs, "Job index x must be between 0 and num_jobs - 1"

	tests_per_job = T // num_jobs
	extras = T % num_jobs

	if x < extras:
		start = x * (tests_per_job + 1)
		end = start + tests_per_job
	else:
		start = extras * (tests_per_job + 1) + (x - extras) * tests_per_job
		end = start + tests_per_job - 1

	return start, end



def extract_total_number_of_tests(processed_fm_eqtl_output_file):
	f = open(processed_fm_eqtl_output_file)
	counter = 0
	for line in f:
		counter = counter + 1
	f.close()
	return counter - 1

###################################
processed_fm_eqtl_output_file = sys.argv[1]
borzoi_downloads_dir = sys.argv[2]
model_train_dir = sys.argv[3]
borzoi_eqtl_output_file = sys.argv[4]
sim_iter = sys.argv[5]
total_sims = sys.argv[6]
###################################

#################################
# For parallelization purposes
#################################
n_tests = extract_total_number_of_tests(processed_fm_eqtl_output_file)
parallel_lb, parallel_ub = get_test_range(int(sim_iter), int(total_sims), n_tests) # [parallel_ub, parallel_lb]



#################################
# Load in stuff
#################################
params_file = model_train_dir + 'micro_models/f0c0/params.json'
targets_file = model_train_dir + 'micro_models/f0c0/data0/targets.txt' #Subset of targets_human.txt

pyfaidx.Faidx(borzoi_downloads_dir +'hg38/assembly/ucsc/hg38.fa')

#Model configuration
seq_len = 393216
rc = True         #Average across reverse-complement prediction
#Read model parameters
with open(params_file) as params_open :
	
	params = json.load(params_open)
	
	params_model = params['model']
	params_train = params['train']

#Remove cropping
#params_model['trunk'][-2]['cropping'] = 0

#Read targets

targets_df = pd.read_csv(targets_file, index_col=0, sep='\t')
target_index = targets_df.index

#Create local index of strand_pair (relative to sliced targets)
if rc :
	strand_pair = targets_df.strand_pair
	
	target_slice_dict = {ix : i for i, ix in enumerate(target_index.values.tolist())}
	slice_pair = np.array([
		target_slice_dict[ix] if ix in target_slice_dict else ix for ix in strand_pair.values.tolist()
	], dtype='int32')

#Load genome fasta and gene annotations

#Initialize fasta sequence extractor
fasta_open = pysam.Fastafile(borzoi_downloads_dir + 'hg38/assembly/ucsc/hg38.fa')

#Load gene/exon annotation
gtf_file = borzoi_downloads_dir + 'hg38/genes/gencode41/gencode41_basic_nort_protein.gtf'

transcriptome = bgene.Transcriptome(gtf_file)

#Get gene span bedtool
bedt_span = transcriptome.bedtool_span()

#Load APA atlas
apa_df = pd.read_csv(borzoi_downloads_dir + 'hg38/genes/polyadb/polyadb_human_v3.csv.gz', sep='\t', compression='gzip')
apa_df = apa_df[['pas_id', 'gene', 'chrom', 'position_hg38', 'strand', 'site_num', 'num_sites', 'site_type', 'pas_type', 'total_count']]

apa_df.loc[apa_df['pas_type'] == 'NoPAS', 'pas_type'] = 'No_CSE'

#Only consider 3' UTR sites
apa_df_utr = apa_df.query("site_type == '3\\' most exon'").copy().reset_index(drop=True)

#Or intronic sites
apa_df_intron = apa_df.query("site_type == 'Intron' and pas_type != 'No_CSE'").copy().reset_index(drop=True)

#Load TSS atlas
tss_df = pd.read_csv(borzoi_downloads_dir + 'hg38/genes/gencode41/gencode41_basic_tss2.bed', sep='\t', names=['chrom', 'position_hg38', 'end', 'tss_id', 'feat1', 'strand'])
tss_df['gene'] = tss_df['tss_id'].apply(lambda x: x.split("/")[1] if "/" in x else x)


#Tracks
track_index = np.nonzero((targets_df['description'] == 'RNA:blood').values)[0].tolist()





#################################
# Run main analysis
#################################

f = open(processed_fm_eqtl_output_file)
t = open(borzoi_eqtl_output_file,'w')
head_count = 0
test_counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'borzoi_mean\n')
		continue
	if test_counter < parallel_lb or test_counter > parallel_ub:
		test_counter = test_counter + 1
		continue
	print(test_counter)
	test_counter = test_counter + 1

	booler = True
	#if booler:
	try:
		chrom = data[0]
		search_gene = data[6].split('.')[0]
		center_pos = int(data[2].split('_')[1])
		poses = [center_pos]
		alts = data[4]
		ref_allele = data[3]


		start = center_pos - seq_len // 2
		end = center_pos + seq_len // 2

		#Get exon bin range
		gene_keys = [gene_key for gene_key in transcriptome.genes.keys() if search_gene in gene_key]

		if len(gene_keys) == 0:
			print('skipped gene due to it not being in the transcriptome')
			continue

		gene = transcriptome.genes[gene_keys[0]]
		gene_strand = gene.strand

		if chrom is None or start is None or end is None :
			chrom = gene.chrom
			g_start, g_end = gene.span()
			mid = (g_start + g_end) // 2
			start = mid - seq_len // 2
			end = mid + seq_len // 2



		# (~6 minutes on CPU w 1 replicate; ~2 minutes on GPU)
		sequence_one_hot_wt = process_sequence(fasta_open, chrom, start, end, seq_len=seq_len)

		#Induce mutation(s)
		sequence_one_hot_mut = np.copy(sequence_one_hot_wt)

		for pos, alt in zip(poses, alts) :
			alt_ix = -1
			if alt == 'A' :
				alt_ix = 0
			elif alt == 'C' :
				alt_ix = 1
			elif alt == 'G' :
				alt_ix = 2
			elif alt == 'T' :
				alt_ix = 3

			sequence_one_hot_mut[pos-start-1] = 0.
			sequence_one_hot_mut[pos-start-1, alt_ix] = 1.

		# Quick error check
		if ref_allele == 'A':
			if sequence_one_hot_wt[pos-start-1][0] != 1:
				print('assumption eroror')
				pdb.set_trace()
		elif ref_allele == 'C':
			if sequence_one_hot_wt[pos-start-1][1] != 1:
				print('assumption eroror')
				pdb.set_trace()
		elif ref_allele == 'G':
			if sequence_one_hot_wt[pos-start-1][2] != 1:
				print('assumption eroror')
				pdb.set_trace()
		elif ref_allele == 'T':
			if sequence_one_hot_wt[pos-start-1][3] != 1:
				print('assumption eroror')
				pdb.set_trace()


		models = []	
		model_file = model_train_dir + 'micro_models/f0c0/train/model_best.h5'

		seqnn_model = seqnn.SeqNN(params_model)
		seqnn_model.restore(model_file, 0)
		seqnn_model.build_slice(target_index)
		if rc :
			seqnn_model.strand_pair.append(slice_pair)
		seqnn_model.build_ensemble(rc, [0])
	
		models.append(seqnn_model)


		#Determine output sequence start
		seq_out_start = start + seqnn_model.model_strides[0]*seqnn_model.target_crops[0]
		seq_out_len = seqnn_model.model_strides[0]*seqnn_model.target_lengths[0]

		#Determine output positions of gene exons
		gene_slice = gene.output_slice(seq_out_start, seq_out_len, seqnn_model.model_strides[0], False, old_version=True)



		#Make predictions
		import time
		t1 = time.time()
		y_wt = predict_tracks(models, sequence_one_hot_wt)
		y_mut = predict_tracks(models, sequence_one_hot_mut)
		t2 = time.time()
		print(t2-t1)
		track_scale = 0.01
		track_transform = 3./4.
		soft_clip = 384.
		untransform_old = True

		#Plot coverage
		qtl_effect_size = extract_eqtl_effect(
			y_wt,
			track_index,
			track_scale,
			track_transform,
			soft_clip,
			y_2_in=y_mut,
			gene_slice=gene_slice,
			untransform_old=untransform_old,
		)


		t.write(line + '\t' + str(qtl_effect_size) + '\n')

	except:
		print('missed test ' + str(test_counter))

	t.flush()
f.close()
t.close()