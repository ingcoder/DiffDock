#!/bin/bash
# Shell script to run the inference code using Python

# Get command-line arguments
PROTEIN_LIGAND_CSV="$1"
OUT_DIR="$2"

# # Activate conda env
source /home/ec2-user/miniconda3/bin/activate diffdock

# # Run the Python script
# # Print the parameters
echo "PROTEIN_LIGAND_CSV: $PROTEIN_LIGAND_CSV"
echo "OUT_DIR: $OUT_DIR"
# echo "INFERENCE_STEPS: $INFERENCE_STEPS"
# echo "SAMPLES_PER_COMPLEX: $SAMPLES_PER_COMPLEX"
echo $PATH
echo $(pwd)

# # Get the directory where this .sh script is located 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# # Change to that directory
cd "$DIR"
echo $(pwd)

python /home/ec2-user/DiffDock/diffdock_postprocess.py "$PROTEIN_LIGAND_CSV" "$OUT_DIR"

# Check if the Python script executed successfully
if [ $? -eq 0 ]; then
    echo "Python script executed successfully, uploading results to S3."

    # Modify the below path and bucket/folder as needed
    aws s3 cp "$PROTEIN_LIGAND_CSV/$OUT_DIR" s3://diffdock/ --recursive

    if [ $? -eq 0 ]; then
        echo "Successfully uploaded results to S3."
    else
        echo "Failed to upload results to S3."
    fi
else
    echo "Python script execution failed, not uploading to S3."
fi
