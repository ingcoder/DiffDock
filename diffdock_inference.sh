#!/bin/bash
# Shell script to run the inference code using Python

# Get command-line arguments
PROTEIN_LIGAND_CSV="$1"
OUT_DIR="$2"
INFERENCE_STEPS="$3"
SAMPLES_PER_COMPLEX="$4"
BATCH_SIZE="$5"
ACTUAL_STEPS="$6"

# Activate conda env
source /home/ec2-user/miniconda3/bin/activate diffdock

# Run the Python script
# Print the variables
echo "PROTEIN_LIGAND_CSV: $PROTEIN_LIGAND_CSV"
echo "OUT_DIR: $OUT_DIR"
echo "INFERENCE_STEPS: $INFERENCE_STEPS"
echo "SAMPLES_PER_COMPLEX: $SAMPLES_PER_COMPLEX"
echo $PATH
echo $(pwd)
# Get the directory of the script itself
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change to that directory
cd "$DIR"

python /home/ec2-user/DiffDock/inference.py --protein_ligand_csv "$PROTEIN_LIGAND_CSV" --out_dir "$OUT_DIR" --inference_steps "$INFERENCE_STEPS" --samples_per_complex "$SAMPLES_PER_COMPLEX" --batch_size "$BATCH_SIZE" --actual_steps "$ACTUAL_STEPS" --no_final_step_noise
