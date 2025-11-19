

#################
# Input data
#################
gtex_target_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/qc/gtex_targets.txt"

baskerville_code_dir="/home/ch271704/tools/baskerville/src/baskerville/scripts/"



#################
# Output directories
#################
# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_run/"

# Directory containing processed tfrecords
tfrecords_dir=${output_root}"tfrecords/"




#################
# Run analysis
#################

# Make TF-records 
# This is basically the Makefile in the borzoi tutorial (https://github.com/calico/borzoi/blob/main/tutorials/latest/make_data/Makefile)
# But in a shell script because i'm not fancy
if false; then
sbatch make_tf_records.sh ${gtex_target_file} ${tfrecords_dir} ${baskerville_code_dir}
fi