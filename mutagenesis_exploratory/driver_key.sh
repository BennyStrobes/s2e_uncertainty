
#######################
# Input data
#######################

# Borzoi boootstrapped model training dir
model_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_bootstrap/model_train/"

# Fasta file
hg38_fasta_file="/home/ch271704/tools/borzoi/examples/hg38/assembly/ucsc/hg38.fa"

# GTex target file used in training
full_gtex_target_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/gtex_targets_full.txt"

raw_mutagenesis_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/mutagenesis_exploratory/input_data/GRCh38_ALL.tsv"


#######################
# Output data
#######################
# Output root directory
root_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/mutagenesis_exploratory/"

# Directory containing processed eqtl sumstats
processed_mutagenesis_dir=${root_dir}"processed_mutagenesis_results/"

# Directory containing borzoi predicted effects
borzoi_eqtl_effects_dir=${root_dir}"borzoi_predicted_effects/"


# Visualize results
visualization_dir=${root_dir}"visualization/"

#######################
# Run analysis
#######################
reorganized_mutagenesis_file=${processed_mutagenesis_dir}"processed_mutagenesis_results.txt"
mutagenesis_variant_vcf=${processed_mutagenesis_dir}"processed_mutagenesis_variants.vcf"
if false; then
sh reorganize_mutagnesis_results.sh $raw_mutagenesis_file $reorganized_mutagenesis_file $mutagenesis_variant_vcf
fi

if false; then
for bs_iter in {1..100}; do
	output_file=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_borzoi_pred_eqtl_effects_borzoi_sed_results.txt"
	sbatch fast_borzoi_sed.sh ${output_file} ${mutagenesis_variant_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi




# Organize results across parallel runs
if false; then
sh organize_mutagenesis_borzoi_sed_results.sh ${borzoi_eqtl_effects_dir} $full_gtex_target_file $borzoi_eqtl_effects_dir $reorganized_mutagenesis_file
fi

##############################################################
# Visualize fine-mapped eQTL analysis
source ~/.bashrc
conda activate borzoi
python visualize_borzoi_estimated_effects.py $borzoi_eqtl_effects_dir"fm_organized_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects.txt" $visualization_dir



