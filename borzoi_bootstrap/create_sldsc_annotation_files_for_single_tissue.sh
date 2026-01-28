#!/bin/bash
#SBATCH -t 0-50:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=12GB  


tissue_name="${1}"
sldsc_processed_anno_dir="${2}"
borzoi_res_file="${3}"
ldsc_snp_annotation_dir="${4}"
genotype_1000G_plink_stem="${5}"
ldsc_code_dir="${6}"
ldsc_weights_dir="${7}"


if false; then
source ~/.bashrc
conda activate borzoi
python create_sldsc_snp_annotation_files_for_single_tissue.py $tissue_name $sldsc_processed_anno_dir $borzoi_res_file $ldsc_snp_annotation_dir $genotype_1000G_plink_stem
fi


source ~/.bashrc
conda activate ldsc
for chrom_num in {1..22}; do
	anno_name=$tissue_name"_borzoi_sig0"
	python ${ldsc_code_dir}ldsc.py --l2 --bfile ${genotype_1000G_plink_stem}${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num}".annot.gz" --out ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num} --print-snps ${ldsc_weights_dir}"hm3_noMHC."${chrom_num}".rsid"
done

for chrom_num in {1..22}; do
	anno_name=$tissue_name"_borzoi_sig1"
	python ${ldsc_code_dir}ldsc.py --l2 --bfile ${genotype_1000G_plink_stem}${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num}".annot.gz" --out ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num} --print-snps ${ldsc_weights_dir}"hm3_noMHC."${chrom_num}".rsid"
done

for chrom_num in {1..22}; do
	anno_name=$tissue_name"_borzoi_sig2"
	python ${ldsc_code_dir}ldsc.py --l2 --bfile ${genotype_1000G_plink_stem}${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num}".annot.gz" --out ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num} --print-snps ${ldsc_weights_dir}"hm3_noMHC."${chrom_num}".rsid"
done

for chrom_num in {1..22}; do
	anno_name=$tissue_name"_borzoi_sig3"
	python ${ldsc_code_dir}ldsc.py --l2 --bfile ${genotype_1000G_plink_stem}${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num}".annot.gz" --out ${sldsc_processed_anno_dir}${anno_name}"."${chrom_num} --print-snps ${ldsc_weights_dir}"hm3_noMHC."${chrom_num}".rsid"
done