#!/bin/bash
# Shell script to run python postpress code with Gninina. 

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

DIFFDOCK_RESULTS_DIR="$1"
OUT_DIR="$2"
PROTEIN_LIGAND_CSV="$3"
FILENAME="$4"
HOME_DIR="/home/ec2-user"
LOG_FILE="/home/ec2-user/DiffDock/my_log_file.log"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # Get the directory where this .sh script is located
CURRENT_USER=$(whoami)
S3_BUCKET="$5"

# Print the parameters
log_info "===== Start diffdock_postprocess.sh ====="
log_info "DIFFDOCK_RESULTS_DIR: $DIFFDOCK_RESULTS_DIR"
log_info "OUT_DIR: $OUT_DIR"
log_info "PROTEIN_LIGAND_CSV: $PROTEIN_LIGAND_CSV"
log_info "PATH: $PATH"
log_info "HOME_DIR: $HOME_DIR"
log_info "Current Directory: $(pwd)"
log_info "Shell script directory: $SCRIPT_DIR"
log_info "CURRENT_USER: $CURRENT_USER"

# Activate conda env 
log_info "Activating conda environment"

# source /home/ec2-user/miniconda3/bin/activate diffdock >> "$LOG_FILE" 2>&1 || log_error "Failed to activate conda environment"
if ! source "$HOME_DIR/miniconda3/bin/activate" diffdock >> "$LOG_FILE" 2>&1; then
    log_error "Failed to source Conda activate script"
fi

log_info "Conda environment successfully activated"
log_info "Conda Path: $PATH"

# Change to that directory
cd "$SCRIPT_DIR" || log_error "Failed to change directory to $PROJECT_DIR"

# Run the Python script
log_info "Running Python script for Gnina Docking"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # Get the directory where this .sh script is located
python "$SCRIPT_DIR/diffdock_postprocess.py" "$DIFFDOCK_RESULTS_DIR" "$OUT_DIR" "$PROTEIN_LIGAND_CSV" "$FILENAME">> "$LOG_FILE" 2>&1 || log_error "Python script execution failed"


# Check if the Python script executed successfully
if [ $? -eq 0 ]; then
    log_info "Python script executed successfully. Start uploading results to S3..."

    # Modify the below path and bucket/folder as needed
    if aws s3 cp "$SCRIPT_DIR/$DIFFDOCK_RESULTS_DIR" "$S3_BUCKET" --recursive >> "$LOG_FILE" 2>&1; then
        log_info "Successfully uploaded results to S3."
    else
        log_error "Failed to upload results to S3."
    fi
else
    ERROR_MESSAGE=$(tail -n 1 "$LOG_FILE")
    log_error "Python script execution failed" "$ERROR_MESSAGE"
    # log_error "Python script execution failed, not uploading to S3."
fi
log_info "===== End diffdock_postprocess.sh ====="
