#################
# Input data
#################


# Directory containing results of borzoi runs
borzoi_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/"

# Directory containing genotype data
plink2_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/process_gtex_genotype_data/processed_genotype/"

# Directory containing eQTL summary statistics
eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/eqtl_sumstats/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"


#################
# Output directories
#################
# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/eqtl_borzoi_integration_processed_data/"

LD_output_dir=${output_root}"LD/"
eqtl_ss_output_dir=${output_root}"marginal_eqtl_ss/"
borzoi_output_dir=${output_root}"borzoi/"





#################
# Run analysis
#################
if false; then
for chrom_num in {1..22}; do
	sbatch integrate_eqtl_borzoi_ld_data.sh ${borzoi_results_dir} ${plink2_genotype_data_dir} ${eqtl_sumstats_dir} ${borzoi_output_dir} ${LD_output_dir} ${eqtl_ss_output_dir} $chrom_num $gtex_v10_pc_genes_gtf
done
fi








