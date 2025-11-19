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






#################
# Command line args
#################
hg38_path = sys.argv[1]
expr_path = sys.argv[2]
tss_data_path = sys.argv[3]
protein_coding_genes_file = sys.argv[4]
model_save_dir = sys.argv[5]
output_stem = sys.argv[6]


model_save_dir = model_save_dir + output_stem + '/'

# make output dir if doesn't exist
os.makedirs(model_save_dir, exist_ok=True)

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
#expr_df[(expr_df != 0.0).any(axis=1)]

# Convert to log scale
#expr_df = np.log(expr_df + 1)

# Restrict genes to protein-coding genes that are present in the expression data and in the metadata
protein_coding_gene_list = np.loadtxt(protein_coding_genes_file,delimiter=',',dtype=str)
use_gene_list = np.intersect1d(np.intersect1d(protein_coding_gene_list, expr_df.index), gene_meta_info.index)

# Split by chromosome to get train, validaiton, and test gene sets
train_genes, val_genes, test_genes = SAGEnet.tools.get_train_val_test_genes(use_gene_list,tss_data_path=tss_data_path, use_enformer_gene_assignments=False)

# Select train and validation gene meta information:
train_gene_meta = gene_meta_info.loc[train_genes]
val_gene_meta = gene_meta_info.loc[val_genes]


# Ben added: not sure why needed
# Make sure on valid chromosomes
#train_gene_meta = train_gene_meta[train_gene_meta['chr'].isin(np.arange(1,23).astype(str))]
#val_gene_meta = val_gene_meta[val_gene_meta['chr'].isin(np.arange(1,23).astype(str))]

# Initialize datasets and dataloaders
input_len=40000
train_dataset = SAGEnet.data.ReferenceGenomeDataset(metadata=train_gene_meta, hg38_file_path=hg38_path, y_data=expr_df, input_len=input_len,majority_seq=False,train_subs=train_individuals)
val_dataset = SAGEnet.data.ReferenceGenomeDataset(metadata=val_gene_meta, hg38_file_path=hg38_path, y_data=expr_df, input_len=input_len,majority_seq=False,train_subs=train_individuals)
train_dataloader = DataLoader(train_dataset,  shuffle=True)
val_dataloader = DataLoader(val_dataset, shuffle=False)


# Initialize an r-SAGE-net model
my_model = SAGEnet.models.rSAGEnet(input_length=input_len)


# Set up for model training
device=1 # which GPU 
max_epochs=50
wandb_logger = WandbLogger(project='test_project_name', name='test_job_name', id='test_job_name', resume="allow", anonymous="must") # change these based on your logging preferences 

es = EarlyStopping(monitor="val_pearson", patience=10,mode='max')
checkpoint_callback = ModelCheckpoint(dirpath=model_save_dir, monitor="val_pearson", save_top_k=1, mode="max", save_last=True, every_n_epochs=1)
lr_monitor = LearningRateMonitor(logging_interval='epoch')
callbacks=[es,checkpoint_callback,lr_monitor]

# ORIG
'''
trainer = pl.Trainer(
accelerator="gpu", 
devices=[int(device)] if device else 1, 
num_nodes=1, 
strategy="ddp" if not device else 'auto', 
callbacks=callbacks, 
max_epochs=max_epochs, 
benchmark=False, 
profiler='simple', 
gradient_clip_val=1, 
logger=wandb_logger, 
log_every_n_steps=10)
'''
trainer = pl.Trainer(
    accelerator="gpu",
    devices=1,                 # use one visible GPU (auto-picks id 0)
    strategy="auto",
    callbacks=callbacks,
    max_epochs=max_epochs,
    benchmark=False,
    profiler='simple',
    gradient_clip_val=1,
    logger=wandb_logger,
    log_every_n_steps=10,
)


# TRAIN
trainer.fit(my_model, train_dataloader, val_dataloaders=val_dataloader)
