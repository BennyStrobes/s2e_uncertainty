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
hg19_path = sys.argv[1]
expr_path = sys.argv[2]
tss_data_path = sys.argv[3]
protein_coding_genes_file = sys.argv[4]
model_save_dir = sys.argv[5]
model_evaluation_dir = sys.argv[6]
output_stem = sys.argv[7]

model_save_dir = model_save_dir + output_stem + '/'


# Load in Tss data file
gene_meta_info = pd.read_csv(tss_data_path, sep="\t",index_col='region_id')
gene_meta_info['chr'] = gene_meta_info['chr_hg38'].str.replace('chr', '', regex=False)
gene_meta_info['tss'] = pd.to_numeric(gene_meta_info['tss_hg38'], errors='coerce').astype('Int64')
# Ben added
gene_meta_info['pos'] = pd.to_numeric(gene_meta_info['tss_hg38'], errors='coerce').astype('Int64')
gene_meta_info = gene_meta_info[gene_meta_info['chr'].isin(np.arange(1,23).astype(str))]


# Load the processed expression data
#orig_expr_df = pd.read_csv(expr_path, skiprows=2, sep='\t')
orig_expr_df = pd.read_csv(expr_path, sep='\t')



# Randomly select 80% of individuals as training individuals
all_expr_data_individauls=orig_expr_df.columns[2:]
n_train_individuals = int(.8*len(all_expr_data_individauls))
np.random.seed(12)
shuffled_indices = np.random.permutation(len(all_expr_data_individauls))
shuffled_individs = np.array(all_expr_data_individauls)[shuffled_indices]  # convert to numpy array for indexing
train_individuals = shuffled_individs[:n_train_individuals]

# Put the expression data into the format required by ReferenceGenomeDataset (indexed by gene IDs, column values are sample names):
expr_df = orig_expr_df[all_expr_data_individauls]
expr_df.index=orig_expr_df['Name'].str.split('.').str[0]

# Filter out zeros??
# expr_df[(expr_df != 0.0).any(axis=1)] # Corr=.31

# Convert to log scale
#expr_df = np.log(expr_df + 1)

# Restrict genes to protein-coding genes that are present in the expression data and in the metadata
protein_coding_gene_list = np.loadtxt(protein_coding_genes_file,delimiter=',',dtype=str)
use_gene_list = np.intersect1d(np.intersect1d(protein_coding_gene_list, expr_df.index), gene_meta_info.index)

# Split by chromosome to get train, validaiton, and test gene sets
train_genes, val_genes, test_genes = SAGEnet.tools.get_train_val_test_genes(use_gene_list,tss_data_path=tss_data_path, use_enformer_gene_assignments=False)


# Select train and validation gene meta information:
test_gene_meta = gene_meta_info.loc[test_genes]
input_len=40000
test_dataset = SAGEnet.data.ReferenceGenomeDataset(metadata=test_gene_meta, hg38_file_path=hg19_path, y_data=expr_df, input_len=input_len,majority_seq=False,train_subs=train_individuals)
test_dataloader = DataLoader(test_dataset, shuffle=False)


# Load in pre-trained model
#UPDATE

ckpt_path = model_save_dir + get_best_model_file_name(model_save_dir)
model = SAGEnet.models.rSAGEnet.load_from_checkpoint(checkpoint_path=ckpt_path)

# Predict
trainer = pl.Trainer(accelerator="gpu", devices=1, strategy="auto", precision=16,logger=False)
preds = trainer.predict(model, test_dataloader)
preds=np.concatenate(preds)

obs = []
for batch_idx, batch_data in enumerate(test_dataloader):
	obs.append(np.asarray(batch_data[1])[0])
obs = np.asarray(obs)

print(np.corrcoef(obs, preds))

