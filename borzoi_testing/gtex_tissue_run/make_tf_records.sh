#!/bin/bash
#SBATCH -t 0-12:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=50GB                         # Memory total in MiB (for all cores)
#SBATCH -c 8              # <-- request 8 CPU cores




gtex_target_file="${1}"
tfrecords_dir="${2}"
baskerville_code_dir="${3}"

source ~/.bashrc
conda activate borzoi


####################
# THIS WAS  TAKEN FROM https://github.com/calico/borzoi/blob/main/tutorials/latest/make_data/Makefile
FASTA_HUMAN="$BORZOI_HG38/assembly/gnomad/hg38.ml.fa"
GAPS_HUMAN="$BORZOI_HG38/assembly/ucsc/hg38_gaps.bed"
UMAP_HUMAN="$BORZOI_HG38/mappability/umap_k36_t10_l32.bed"
BLACK_HUMAN="$BORZOI_HG38/blacklist/blacklist_hg38_all.bed"

FASTA_MOUSE="$BORZOI_MM10/assembly/ucsc/mm10.ml.fa"
GAPS_MOUSE="$BORZOI_MM10/assembly/ucsc/mm10_gaps.bed"
UMAP_MOUSE="$BORZOI_MM10/mappability/umap_k36_t10_l32.bed"
BLACK_MOUSE="$BORZOI_MM10/blacklist/blacklist_mm10_all.bed"

ALIGN="$BORZOI_HG38/align/hg38.mm10.syn.net.gz"

OUT=${tfrecords_dir}

# mini borzoi configuration
LENGTH=393216
TSTRIDE=131087   # 393216/3 - 15
CROP=0
WIDTH=32
FOLDS=8

AOPTS="--break 2097152 -c $CROP --nf 524288 --no 393216 -l $LENGTH --stride $TSTRIDE -f $FOLDS --umap_t 0.5 -w $WIDTH"
DOPTS="-c $CROP -d 2 -f $FOLDS -l $LENGTH -p 64 -r 16 --umap_clip 0.5 -w $WIDTH"
###########################################

if false; then
cat "$UMAP_HUMAN" "$BLACK_HUMAN" \
  | awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' \
  | bedtools sort -i - \
  | bedtools merge -i - \
  > ${tfrecords_dir}umap_human.bed
fi

if false; then
cat "$UMAP_MOUSE" "$BLACK_MOUSE" \
  | awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' \
  | bedtools sort -i - \
  | bedtools merge -i - \
  > ${tfrecords_dir}umap_mouse.bed
fi



if false; then
python ${baskerville_code_dir}hound_data_align.py \
  -a hg38,mm10 \
  -g "$GAPS_HUMAN","$GAPS_MOUSE" \
  -u ${tfrecords_dir}"umap_human.bed",${tfrecords_dir}"umap_mouse.bed" \
  $AOPTS \
  -o "$OUT" \
  "$ALIGN" \
  "$FASTA_HUMAN","$FASTA_MOUSE"
fi






python ${baskerville_code_dir}hound_data.py \
  --restart \
  $DOPTS \
  -b "$BLACK_HUMAN" \
  -o "$OUT/hg38" \
  "$FASTA_HUMAN" \
  -u {tfrecords_dir}umap_human.bed \
  ${gtex_target_file} \
  --local


