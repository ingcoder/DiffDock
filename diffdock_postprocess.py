import os
import re
import subprocess
import pandas as pd
from glob import glob
from shlex import quote
from datetime import datetime
from tqdm.auto import tqdm
from rdkit import Chem

def run_command(command):
    # process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    # stdout, _ = process.communicate()
    # return stdout.decode('utf-8')

    try:
        # Run the command and capture stdout and stderr
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Wait for the command to complete and get the output
        stdout, stderr = process.communicate()

        # Check if the command was successful
        if process.returncode != 0:
            # The command failed, raise an exception with stderr
            raise subprocess.CalledProcessError(process.returncode, command, stderr)

        # Decode and return the stdout
        return stdout.decode('utf-8')

    # except subprocess.CalledProcessError as e:
    #     # Handle subprocess errors (e.g., command failed, invalid command, etc.)
    #     error_message = e.stderr.decode('utf-8')
    #     print(f"Command failed with error: {error_message}")
    #     return None

    except Exception as e:
        # Handle other exceptions
        print(f"An error occurred: {str(e)}")
        return None

def main(results_directory, download_results=False):
    print("inside main")
    os.chdir(results_directory)
    results_dirs = glob("./complex*")

    results_pdb_df = pd.read_csv('/home/ec2-user/DiffDock/data/protein_ligand_example_csv.csv')
    print("results_pdb_df", results_pdb_df)
    results_pdb_file = [os.path.join('/home/ec2-user/DiffDock', file_path) for file_path in results_pdb_df['protein_path']][0]
    print('results pdb file:', results_pdb_file)

    # Remove HETATM records from the PDB file for further processing
    results_pdb_file_no_hetatms = f"{results_pdb_file}_nohet.pdb"
    cmd = f"grep -v '^HETATM' {results_pdb_file} > {results_pdb_file_no_hetatms}"
    run_command(cmd)

    rows = []

    for results_dir in tqdm(results_dirs, desc="runs"):
    #     # Extract the PDB filename and molecule's SMILES representation from the directory name
    #     # Replace. Instead extract from $PROTEIN_LIGAND_CSV"
    #     # results_pdb_file = "/tmp/pdb/" + re.findall("tmp-pdb-(.+\.pdb)", results_dir)[0] # extract from csv
    #     # results_smiles = re.findall("pdb_+(.+)", results_dir)[0] # create from sdf file

        # results_smiles = None
        # def sdf_to_smiles(sdf_filename):
        #     supplier = Chem.SDMolSupplier(sdf_filename)
        #     for mol in supplier:
        #         if mol:  # Ensure the molecule is not None
        #             return Chem.MolToSmiles(mol)
        #     return None

        # sdf_filename = 'path_to_your_sdf_file.sdf'
        # smiles = sdf_to_smiles(sdf_filename)
        # print(smiles)

        
        # Collect all SDF files corresponding to molecule conformations with confidence values in their filenames
        results_sdfs = [os.path.join(results_dir, f) for f in os.listdir(results_dir) if "confidence" in f and f.endswith(".sdf")]
       
       
        # Print the current directory
        current_directory = os.getcwd()
        print(f"Current directory: {current_directory}")
        print(len(results_sdfs))

    #     # Remove HETATM records from the PDB file for further processing
    #     results_pdb_file_no_hetatms = f"{results_pdb_file}_nohet.pdb"
    #     cmd = f"grep -v '^HETATM' {results_pdb_file} > {results_pdb_file_no_hetatms}"
    #     run_command(cmd)

        for results_sdf in tqdm(results_sdfs[:2], leave=False, desc="files"):
            # Extract confidence value associated with the conformation from the SDF filename
            confidence = re.findall("confidence([\-\.\d]+)\.sdf", results_sdf)[0]
            print('Confidence', confidence)

            # Get the initial predicted binding affinity of the molecule using gnina
            cmd = f"/home/ec2-user/gnina --score_only -r '{results_pdb_file_no_hetatms}' -l '{results_sdf}'"
            scored_stdout = run_command(cmd)
            print(scored_stdout)
    #         scored_affinity = re.findall("Affinity:\s*([\-\.\d+]+)", scored_stdout)[0]
    #         print('Scored affinity', scored_affinity)

    #         # Optimize the molecule's pose locally and get the new predicted binding affinity using gnina
    #         cmd = f"/home/ec2-user/gnina --local_only --minimize -r '{results_pdb_file_no_hetatms}' -l '{results_sdf}' --autobox_ligand '{results_sdf}' --autobox_add 2"
    #         minimized_stdout = run_command(cmd)
    #         minimized_affinity = re.findall("Affinity:\s*([\-\.\d+]+)", minimized_stdout)[0]

    #         # Store the results for this molecule conformation in the rows list
    #         rows.append((results_pdb_file.split('/')[-1], results_smiles, float(confidence), float(scored_affinity), float(minimized_affinity), results_sdf))
    #     print(rows)

    # # # Convert results to a Pandas DataFrame and save to a TSV file
    # df_results = pd.DataFrame(rows, columns=["pdb_file", "smiles", "diffdock_confidence", "gnina_scored_affinity", "gnina_minimized_affinity", "sdf_file"])
    # df_results_tsv = "df_diffdock_results.tsv"
    # df_results.to_csv(df_results_tsv, sep='\t', index=None)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process the results.')
    parser.add_argument('results_directory', type=str, help='Path to the results directory.')
    parser.add_argument('--download', dest='download_results', action='store_true', help='Flag to initiate download of results.')
    args = parser.parse_args()
    print('main', args.results_directory)
   
    
    main(args.results_directory)