########################################################
# This will be a bit messy because half on o2, and half on e3
# o2 because thats where individual level UKBB data is
# e3 because thats where we can run borzoi
########################################################

###########
#* CONSIDER MATDCHING SAMPLES TO UKBB BURDEN TEST SAMPLES


#####################
# Input data (o2)
#####################
burden_test_data_dir="/n/groups/price/UKBiobank/sumstats/bolt_337K_unrelStringentBrit_burden_WES/"


ukbb_pheno_file_v2="/n/groups/price/UKBiobank/app19808mosaic/bloodQC/ukb4777.blood_v2.covars.tab"
ukbb_pheno_file_v2="/n/groups/price/UKBiobank/app10438assoc/phenotype_steven/UKB_v3.061518.tab"
ukbb_pheno_file_v3="/n/groups/price/UKBiobank/app10438assoc/ukb4777.processed_and_post.plinkPCs.tab.gz"
ukbb_pheno_file_v4="/n/groups/price/UKBiobank/app10438assoc/phenotype_steven/UKB_biochemistry.051619.tab"




genotype_dir="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"
genotype_sample_names="/n/groups/price/UKBiobank/download_500K/ukb1404_imp_chr1_v3_withdrawn3.sample"


wb_unrelated_samples_file="/n/groups/price/UKBiobank/sampleQC/samples_337K.txt"

trait_names_file="/n/groups/price/ben/non_coding_burden/input_data/trait_names.txt"

gene_annotation_summary_file="/n/groups/price/ben/gene_annotation_files/gencode.v37lift37.genes.tsv"

high_info_snp_dir="/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/snp_info/"

#####################
# Input data (e3)
#####################

borzoi_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/"

liftover_chain="/lab-share/CHIP-Strober-e2/Public/ben/tools/hg19ToHg38.over.chain.gz"


#####################
# Output data (o2)
#####################
o2_output_root="/n/groups/price/ben/non_coding_burden/"

sig_burden_genes_dir=${o2_output_root}"sig_burden_genes/"

variant_gene_pairs_dir=${o2_output_root}"variant_gene_pairs/"

o2_s2e_results_dir=${o2_output_root}"s2e_results/"

non_coding_burden_results_dir=${o2_output_root}"non_coding_burden/"

visualization_dir=${o2_output_root}"visualization/"



#####################
# Output data (e3)
#####################
e3_output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/non_coding_burden_exploratory/"

# input data resulting from o2
e3_input_data_dir=${e3_output_root}"e3_input_data_dir/"

# Directory containing vcf files to test borzoi
vcf_dir=${e3_output_root}"vcfs/"

# Directory containing results of borzoi run
borzoi_pred_dir=${e3_output_root}"borzoi_preds/"

# Organized borzoi prediction results
organized_borzoi_pred_dir=${e3_output_root}"organized_borzoi_pred_results/"


#####################
# Run code
#####################


#########################
# o2 chunk
#########################
if false; then
sh extract_significant_burden_genes.sh $trait_names_file $burden_test_data_dir $sig_burden_genes_dir $gene_annotation_summary_file
fi
burden_genes_summary_file=${sig_burden_genes_dir}"sig_burden_genes_PTV_0.01.txt"

if false; then
for chrom_num in {1..22}; do
	sbatch extract_variant_gene_pairs_to_test.sh $burden_genes_summary_file $genotype_dir $variant_gene_pairs_dir $chrom_num
done
fi

if false; then
source ~/.bash_profile
python organize_variant_gene_pairs_across_chrom_runs.py $variant_gene_pairs_dir
fi
#########################
# e3 chunk
#########################
# File created from o2
vg_pairs_to_test_file=${e3_input_data_dir}"burden_variant_gene_pairs_dist_25000.txt"

# Convert variant gene pairs to variant file to input to borzoi
variant_vcf_file=${vcf_dir}"burdent_variants.vcf"
if false; then
source ~/.bashrc
conda activate borzoi
python convert_gene_pair_file_to_variant_input_file.py $vg_pairs_to_test_file $variant_vcf_file
fi

variant_vcf_hg38_file=${vcf_dir}"burdent_variants_hg38.vcf"
if false; then
sh liftover_variants.sh ${vcf_dir}"burdent_variants.vcf" $variant_vcf_hg38_file $liftover_chain $vcf_dir
fi

# Split into 4 equal sized chunks
if false; then
split -n l/4 "$variant_vcf_hg38_file" "${variant_vcf_hg38_file}.part."
fi


# Run borzoi
if false; then
for model_num in {0..3}; do
	split_number="aa"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number

	split_number="ab"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number

	split_number="ac"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number

	split_number="ad"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number
done
fi

if false; then
model_num="1"
	split_number="ab"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number
	split_number="ac"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number
	split_number="ad"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number
model_num="2"
	split_number="aa"
	sbatch borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_"${split_number} ${vcf_dir}"burdent_variants_hg38.vcf.part."${split_number} $borzoi_training_dir $model_num $split_number
fi

if false; then
organized_pred_results_file=${organized_borzoi_pred_dir}"organized_variant_gene_pair_pred_results.txt"
sh organize_borzoi_sed_results.sh $borzoi_pred_dir"model_" $organized_pred_results_file $vg_pairs_to_test_file
fi



#########################
# Back to o2
#########################
# Copied over from e3
organized_pred_results_file_o2=${o2_s2e_results_dir}"organized_variant_gene_pair_pred_results.txt"


if false; then
for chrom_num in {1..22}; do
	burden_results_stem=${non_coding_burden_results_dir}"non_coding_burden_tests_chr"${chrom_num}
	sbatch run_non_coding_burden_tests.sh $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $burden_results_stem
done
fi

if false; then
for chrom_num in {1..22}; do
	burden_results_stem=${non_coding_burden_results_dir}"non_coding_burden_tests_well_imputed_chr"${chrom_num}
	sbatch run_non_coding_burden_tests_well_imputed.sh $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $burden_results_stem $high_info_snp_dir
done
fi

if false; then
source ~/.bash_profile
python concatenate_non_coding_burden_tests.py ${non_coding_burden_results_dir}"non_coding_burden_tests"
fi


if false; then
source ~/.bash_profile
Rscript visualize_non_coding_burden_results.R ${non_coding_burden_results_dir}"non_coding_burden_tests.txt" $visualization_dir
fi



