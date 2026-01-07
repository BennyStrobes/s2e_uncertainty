

#######################
# Input data
#######################

# Sardinia eqtl summary stat file
sardinia_raw_sumstat_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rare_var_exploratory/input_data/eQTLs.Leuko.full_summary_stats.tsv.gz"

# hg19 genotype dir
hg19_1000G_genotype_dir="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg19/"


# Borzoi boootstrapped model training dir
model_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_bootstrap/model_train/"

# Fasta file
hg38_fasta_file="/home/ch271704/tools/borzoi/examples/hg38/assembly/ucsc/hg38.fa"

full_gtex_target_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/gtex_targets_full.txt"


#######################
# Output data
#######################
# Output root directory
root_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rare_var_exploratory/"

# Directory containing processed eqtl sumstats
processed_eqtl_sumstats_dir=${root_dir}"processed_eqtl_sumstats/"

# Directory containing borzoi predicted effects
borzoi_eqtl_effects_dir=${root_dir}"borzoi_predicted_effects/"




#######################
# Run analysis
#######################

# Fine-mapping eqtls using single causal variant fine-mapping
pip_threshold="0.9"
fm_eqtl_sumstats_file=${processed_eqtl_sumstats_dir}"sardinia_scv_fm_eqtls_pip_"${pip_threshold}".txt"
if false; then
sh fine_map_eqtl_sumstats.sh $sardinia_raw_sumstat_file $fm_eqtl_sumstats_file $pip_threshold
fi

# Add AF to fine-mapped QTLs files
fm_and_af_eqtl_sumstats_file=${processed_eqtl_sumstats_dir}"sardinia_scv_fm_eqtls_pip_"${pip_threshold}_w_af".txt"
if false; then
sh add_allele_frequency_info_to_fm_file.sh $fm_eqtl_sumstats_file $fm_and_af_eqtl_sumstats_file $hg19_1000G_genotype_dir
fi
# NOte: a bit manual as a filled in the MAF NAs via webportal gnomad


# Liftover
hg19_variant_id_file=${processed_eqtl_sumstats_dir}"hg19_variant_ids.txt"
if false; then
source ~/.bashrc
conda activate borzoi
python create_variant_id_file.py $fm_and_af_eqtl_sumstats_file $hg19_variant_id_file
fi
# Run Liftover on ucsc genome browser to get new variant positions
liftover_result=${processed_eqtl_sumstats_dir}"hglft_genome_13bf3b_407810.bed"

fm_and_af_eqtl_sumstats_hg38_file=${processed_eqtl_sumstats_dir}"sardinia_scv_fm_eqtls_pip_"${pip_threshold}"_w_af_hg38.txt"
sardinia_hg38_vcf=${processed_eqtl_sumstats_dir}"sardinia_scv_fm_eqtls_pip_"${pip_threshold}"_hg38.vcf"
if false; then
source ~/.bashrc
conda activate borzoi
python create_new_eqtl_summary_file_in_hg38.py $fm_and_af_eqtl_sumstats_file $liftover_result $fm_and_af_eqtl_sumstats_hg38_file $sardinia_hg38_vcf $hg38_fasta_file
fi

if false; then
for bs_iter in {1..100}; do
	output_file=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results.txt"
	sbatch fast_borzoi_sed.sh ${output_file} ${sardinia_hg38_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi

# Organize results across parallel runs
sh organize_sardinia_borzoi_sed_results_fm.sh ${borzoi_eqtl_effects_dir} $full_gtex_target_file $borzoi_eqtl_effects_dir $fm_and_af_eqtl_sumstats_hg38_file




