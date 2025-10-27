##################
# Input data
##################
# Downloaded from (on 10/20/25):
# curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz  
# gunzip hg19.fa.gz
fasta_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/GEUVADIS/hg19.fa"
# Downloaded from (on 10/20/25):
# https://github.com/ni-lab/finetuning-enformer/tree/main/process_geuvadis_data/log_tpm/corrected_log_tpm.annot.csv.gz
expr_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/GEUVADIS/corrected_log_tpm.annot.csv.gz"
# Downloaded from https://github.com/mostafavilabuw/SAGEnet/tree/main/input_data on 10/20/25
tss_data_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/gene-ids-and-positions.tsv"
# Downloaded from https://github.com/mostafavilabuw/SAGEnet/tree/main/input_data on 10/20/25
protein_coding_gene_list="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/protein_coding_genes.csv"

##################
# Output data
##################
# Directory to keep track of model fits
model_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/model_training/"
# directory containing evaluation results
model_evaluation_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/model_evaluation/"



##################
# Load in environment to run SAGE
##################
if false; then
conda activate SAGEnet
fi



##################
# Run analysis
##################
# Train model based on rsage notebook
if false; then
sbatch train_rsage_based_on_notebook.sh $fasta_file $expr_file $tss_data_file $protein_coding_gene_list $model_training_dir
fi

# Evaluate model based on rsage notebook
if false; then
sh evaluate_rsage_based_on_notebook.sh $fasta_file $expr_file $tss_data_file $protein_coding_gene_list $model_training_dir $model_evaluation_dir
fi