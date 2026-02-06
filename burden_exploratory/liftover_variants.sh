#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-8:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu-pe                           # Partition to run in
#SBATCH --mem=10GB  



hg19_vcf=${1}
hg38_vcf=${2}
liftover_chain=${3}
vcf_dir=${4}

source ~/.bashrc
conda activate plink_env

# 1) TSV -> BED (0-based start, 1-based end)
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' $hg19_vcf > ${vcf_dir}variants.b37.bed



# 3) Lift
CrossMap bed $liftover_chain ${vcf_dir}variants.b37.bed ${vcf_dir}variants.hg38.bed
# unmapped will be written as variants.b37.bed.unmap

# 4) BED -> hg38 positions (POS is the 3rd column, i.e. end)
awk 'BEGIN{OFS="\t"} {print $1, $3, $4}' ${vcf_dir}variants.hg38.bed > ${vcf_dir}lifted.hg38.pos.tsv



awk 'BEGIN{OFS="\t"} {print $3,$0}' "$hg19_vcf" | sort -k1,1 > "${vcf_dir}orig.byid.tsv"
sort -k3,3 "${vcf_dir}lifted.hg38.pos.tsv" > "${vcf_dir}lifted.sorted.tsv"

join -t $'\t' -1 1 -2 3 "${vcf_dir}orig.byid.tsv" "${vcf_dir}lifted.sorted.tsv" \
| awk 'BEGIN{OFS="\t"} {print $7, $8, $1, $5, $6}' > "$hg38_vcf"