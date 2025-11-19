import numpy as np 
import pandas as pd
import SAGEnet.data 
import SAGEnet.tools
from SAGEnet.models import pSAGEnet
from torch.utils.data import DataLoader
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping, LearningRateMonitor
import pytorch_lightning as pl
from pytorch_lightning.loggers import WandbLogger
import glob
import os
import pysam
import pdb
import sys
from scipy.stats import pearsonr


def get_best_model_file_name(model_save_dir):
	best_file_name = 'NANA'
	counter = 0
	for file_name in os.listdir(model_save_dir):
		if file_name.startswith("epoch"):
			best_file_name = file_name
			counter = counter +1

	if counter != 1:
		print('assumptinoe rororor')
		pdb.set_trace()

	if best_file_name == 'NANA':
		print('assumption eorororo')
		pdb.set_trace()

	return best_file_name


#################
# Command line args
#################
processed_fm_eqtl_output_file = sys.argv[1]
hg38_path = sys.argv[2]
expr_path = sys.argv[3]
tss_data_path = sys.argv[4]
protein_coding_genes_file = sys.argv[5]
model_stem = sys.argv[6]
output_file = sys.argv[7]
input_len=40000


# Load in Tss data file
gene_meta_info = pd.read_csv(tss_data_path, sep="\t",index_col='region_id')
gene_meta_info['chr'] = gene_meta_info['chr_hg38'].str.replace('chr', '', regex=False)
gene_meta_info['tss'] = pd.to_numeric(gene_meta_info['tss_hg38'], errors='coerce').astype('Int64')
# Ben added
gene_meta_info['pos'] = pd.to_numeric(gene_meta_info['tss_hg38'], errors='coerce').astype('Int64')
gene_meta_info = gene_meta_info[gene_meta_info['chr'].isin(np.arange(1,23).astype(str))]


# Restrict genes to protein-coding genes that are present in the expression data and in the metadata
protein_coding_gene_list = np.loadtxt(protein_coding_genes_file,delimiter=',',dtype=str)
use_gene_list = np.intersect1d(protein_coding_gene_list,  gene_meta_info.index)
valid_genes = {}
for gene in use_gene_list:
	valid_genes[gene] = 1

test_gene_meta = gene_meta_info.loc[use_gene_list]
test_gene_meta['region_id'] = test_gene_meta.index

# Create mapping from gene to tss
gene_to_tss = {}
ordered_genes = np.asarray(test_gene_meta['region_id'])
ordered_tss = np.asarray(test_gene_meta['tss'])
for ii, gene_name in enumerate(ordered_genes):
	gene_to_tss[gene_name] = ordered_tss[ii]


# Create variant info data frame
var_info_gene_names = []
var_info_chr = []
var_info_pos = []
var_info_ref = []
var_info_alt = []
var_eqtl_beta = []
f = open(processed_fm_eqtl_output_file)
skipped = 0
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	line_gene_id = data[6].split('.')[0]
	if line_gene_id not in valid_genes:
		continue
	gene_tss = gene_to_tss[line_gene_id]

	line_chrom_num = data[0].split('hr')[1]
	line_pos = data[1]
	line_alt_allele = data[4]
	line_ref_allele = data[3]

	if np.abs(gene_tss-int(line_pos)) > (input_len//2):
		skipped = skipped + 1
		continue

	var_info_gene_names.append(line_gene_id)
	var_info_chr.append(line_chrom_num)
	var_info_pos.append(line_pos)
	var_info_ref.append(line_ref_allele)
	var_info_alt.append(line_alt_allele)
	var_eqtl_beta.append(float(data[-2]))
f.close()
var_eqtl_beta = np.asarray(var_eqtl_beta)
print(skipped)
variant_info = pd.DataFrame({'region_id':var_info_gene_names, 'chr':var_info_chr, 'pos':np.asarray(var_info_pos).astype(int), 'ref':var_info_ref, 'alt':var_info_alt})


variant_dataset = SAGEnet.data.VariantDataset(metadata=test_gene_meta, hg38_file_path=hg38_path,variant_info=variant_info, input_len=input_len, single_seq=True, insert_variants=True)
variant_dataloader = DataLoader(variant_dataset, shuffle=False)

ref_dataset = SAGEnet.data.VariantDataset(metadata=test_gene_meta, hg38_file_path=hg38_path,variant_info=variant_info, input_len=input_len, single_seq=True, insert_variants=False)
ref_dataloader = DataLoader(ref_dataset, shuffle=False)



###############
# DEBUGGING
################
variant_seqs = []
for info in variant_dataloader:
	variant_seqs.append(np.asarray(info[0][0,:,:]))

ref_seqs = []
for info in ref_dataloader:
	ref_seqs.append(np.asarray(info[0][0,:,:]))

for ii in range(len(variant_seqs)):
	if np.sum((variant_seqs[ii] - ref_seqs[ii]) != 0) != 2:
		print('assumption eroror')
		pdb.set_trace()




# Load in pre-trained model
#UPDATE
model = SAGEnet.models.rSAGEnet.load_from_checkpoint(checkpoint_path=model_stem + get_best_model_file_name(model_stem))

# Predict
trainer = pl.Trainer(accelerator="gpu", devices=1, strategy="auto", precision=16,logger=False)
variant_preds = trainer.predict(model, variant_dataloader)
variant_preds=np.concatenate(variant_preds)

# Predict
trainer = pl.Trainer(accelerator="gpu", devices=1, strategy="auto", precision=16,logger=False)
ref_preds = trainer.predict(model, ref_dataloader)
ref_preds=np.concatenate(ref_preds)
deltas = variant_preds - ref_preds

print(pearsonr(deltas, var_eqtl_beta))

