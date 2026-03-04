#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





ukbb_pheno_file_v2="${1}"
ukbb_pheno_file_v3="${2}"
ukbb_pheno_file_v4="${3}"
organized_pred_results_file_o2="${4}"
burden_genes_summary_file="${5}"
genotype_dir="${6}"
trait_names_file="${7}"
chrom_num="${8}"
genotype_sample_names="${9}"
wb_unrelated_samples_file="${10}"
burden_results_stem="${11}"

source ~/.bash_profile


borzoi_effect_thresh="0.1"
burden_output_file=$burden_results_stem"_thresh_"$borzoi_effect_thresh".txt"
python run_non_coding_burden_tests.py $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $borzoi_effect_thresh $burden_output_file

borzoi_effect_thresh="0.075"
burden_output_file=$burden_results_stem"_thresh_"$borzoi_effect_thresh".txt"
python run_non_coding_burden_tests.py $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $borzoi_effect_thresh $burden_output_file

borzoi_effect_thresh="0.05"
burden_output_file=$burden_results_stem"_thresh_"$borzoi_effect_thresh".txt"
python run_non_coding_burden_tests.py $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $borzoi_effect_thresh $burden_output_file

borzoi_effect_thresh="0.15"
burden_output_file=$burden_results_stem"_thresh_"$borzoi_effect_thresh".txt"
python run_non_coding_burden_tests.py $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $borzoi_effect_thresh $burden_output_file

borzoi_effect_thresh="0.2"
burden_output_file=$burden_results_stem"_thresh_"$borzoi_effect_thresh".txt"
python run_non_coding_burden_tests.py $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $borzoi_effect_thresh $burden_output_file

borzoi_effect_thresh="0.25"
burden_output_file=$burden_results_stem"_thresh_"$borzoi_effect_thresh".txt"
python run_non_coding_burden_tests.py $ukbb_pheno_file_v2 $ukbb_pheno_file_v3 $ukbb_pheno_file_v4 $organized_pred_results_file_o2 $burden_genes_summary_file $genotype_dir $trait_names_file $chrom_num $genotype_sample_names $wb_unrelated_samples_file $borzoi_effect_thresh $burden_output_file
