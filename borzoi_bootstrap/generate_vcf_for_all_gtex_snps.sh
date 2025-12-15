#!/bin/bash
#SBATCH -t 0-15:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB 




gtex_sumstat_dir="${1}"
output_vcf_file="${2}"

source ~/.bashrc
conda activate borzoi

python generate_vcf_for_all_gtex_snps.py $gtex_sumstat_dir $output_vcf_file



file=$output_vcf_file
lines=$(wc -l < "$file")
base=$((lines / 4))
rem=$((lines % 4))

awk -v base="$base" -v rem="$rem" '
{
  # first "rem" chunks get one extra line
  target = (chunk < rem) ? base + 1 : base
  print > sprintf("%s.part.%02d", FILENAME, chunk)
  count++
  if (count == target) {
    chunk++
    count = 0
  }
}
' "$file"