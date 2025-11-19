

################
# Input data
################
borzoi_micro_json_file="/home/ch271704/tools/borzoi/tutorials/latest/train_model/params_micro.json"

training_input_data_dir="/home/ch271704/tools/borzoi/tutorials/latest/make_data/data/hg38"


################
# Output data
################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_testing/tutorial/"
model_training_dir=${output_root}"model_training/"





if false; then
source ~/.bashrc
conda activate borzoi
fi

if false; then
python westminster_train_folds_custom.py -e borzoi -f 2 -c 1 -o ${model_training_dir}micro_models ${borzoi_micro_json_file} ${training_input_data_dir}
fi

if false; then
sbatch practice_sub.sh
fi