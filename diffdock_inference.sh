#!/bin/bash
# Shell script to run the inference code using Python

# Exit immediately if a command exits with a non-zero status.
set -e

# Function to print info messages
log_info() {
    echo -e "[INFO]- $1"
}

# Function to print error messages
log_error() {
    echo -e "[ERROR] - $1" >&2
    if [ -n "$2" ]; then
        echo -e "[ERROR DETAILS] - $2" >&2
    fi
    exit 1
}

# Get command-line arguments
PROTEIN_LIGAND_CSV="$1"
RESULTS_DIR="$2"
INFERENCE_STEPS="$3"
SAMPLES_PER_COMPLEX="$4"
BATCH_SIZE="$5"
ACTUAL_STEPS="$6"
HOME_DIR="/home/ec2-user"
LOG_FILE="/var/log/my_log_file.log"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # Get the directory of the script itself


# Print the variables
log_info "===== Start diffdock_inference.sh ====="
log_info "PROTEIN_LIGAND_CSV: $PROTEIN_LIGAND_CSV"
log_info "RESULTS_DIR: $RESULTS_DIR"
log_info "INFERENCE_STEPS: $INFERENCE_STEPS"
log_info "SAMPLES_PER_COMPLEX: $SAMPLES_PER_COMPLEX"
log_info "PATH: $PATH"
log_info "Current directory: $(pwd)"
log_info "HOME_DIR: $HOME_DIR"
log_info "Shell script directory: $SCRIPT_DIR"

# Activate conda env
log_info "Activating conda environment"

# source /home/ec2-user/miniconda3/bin/activate diffdock >> "$LOG_FILE" 2>&1 || log_error "Failed to activate conda environment"
if ! source "$HOME_DIR/miniconda3/bin/activate" diffdock >> "$LOG_FILE" 2>&1; then
    log_error "Failed to source Conda activate script"
fi

log_info "Conda environment successfully activated"
log_info "Conda Path: $PATH"

# Change to that directory ~/DiffDock
cd "$SCRIPT_DIR"

# Run the Python script
log_info "Running Python script for DiffDock Inference..."
python "$SCRIPT_DIR/inference.py" --protein_ligand_csv "$PROTEIN_LIGAND_CSV" --out_dir "$RESULTS_DIR" --inference_steps "$INFERENCE_STEPS" --samples_per_complex "$SAMPLES_PER_COMPLEX" --batch_size "$BATCH_SIZE" --actual_steps "$ACTUAL_STEPS" --no_final_step_noise
log_info "===== End diffdock_inference.sh ====="