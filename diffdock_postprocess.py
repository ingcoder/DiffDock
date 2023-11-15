import os
import re
import subprocess
import pandas as pd
from glob import glob
from tqdm.auto import tqdm
from rdkit import Chem
import argparse


def run_command(command):
    """Executes a shell command and returns the output."""
    try:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command, stderr)
        return stdout.decode('utf-8')
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


def sdf_to_smile(sdf_filename):
    """Converts an SDF file to a SMILES string."""
    supplier = Chem.SDMolSupplier(sdf_filename)
    mol = next(supplier, None)
    return Chem.MolToSmiles(mol) if mol else None


def remove_hetatms(pdb_file):
    """Removes HETATM records from a PDB file."""
    pdb_file_no_hetatms = f"{pdb_file}_nohet.pdb"
    cmd = f"grep -v '^HETATM' {pdb_file} > {pdb_file_no_hetatms}"
    run_command(cmd)
    return pdb_file_no_hetatms


def process_results_directory(results_directory, results_pdb_file_no_hetatms):
    """Processes each result directory and returns a DataFrame with the results."""
    df_rows = []
    results_dirs = glob(os.path.join(results_directory, "complex*"))

    for results_dir in tqdm(results_dirs, desc="runs"):
        results_sdfs = [os.path.join(results_dir, f) for f in os.listdir(results_dir) if "confidence" in f and f.endswith(".sdf")]

        for results_sdf in tqdm(results_sdfs[:2], leave=False, desc="files"):
            process_sdf_file(results_sdf, results_pdb_file_no_hetatms, df_rows)

    return pd.DataFrame(df_rows, columns=["pdb_file", "smile", "diffdock_confidence", "gnina_scored_affinity", "gnina_minimized_affinity", "sdf_file"])


def process_sdf_file(results_sdf, results_pdb_file_no_hetatms, df_rows):
    """Processes an SDF file and appends the results to the provided list."""
    smile = sdf_to_smile(results_sdf)
    confidence = re.findall("confidence([\-\.\d]+)\.sdf", results_sdf)[0]
    scored_affinity = get_affinity(results_pdb_file_no_hetatms, results_sdf, score_only=True)
    minimized_affinity = get_affinity(results_pdb_file_no_hetatms, results_sdf, local_only=True)
    df_rows.append((os.path.basename(results_sdf), smile, float(confidence), scored_affinity, minimized_affinity, results_sdf))


def get_affinity(pdb_file, sdf_file, score_only=False, local_only=False):
    """Gets the affinity using gnina."""
    cmd = f"/home/ec2-user/gnina --{'score_only' if score_only else 'local_only --minimize'} -r '{pdb_file}' -l '{sdf_file}'"
    if local_only:
        cmd += f" --autobox_ligand '{sdf_file}' --autobox_add 2"
    stdout = run_command(cmd)
    return float(re.findall("Affinity:\s*([\-\.\d+]+)", stdout)[0])


def save_results(df_results, output_directory):
    """Saves the results DataFrame to a TSV file."""
    tsv_file = os.path.join(output_directory, "df_diffdock_results.tsv")
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    df_results.to_csv(tsv_file, sep='\t', index=None)
    return tsv_file


def main(results_directory, output_directory):
    # Load PDB file path
    results_pdb_df = pd.read_csv('/home/ec2-user/DiffDock/data/protein_ligand_example_csv.csv')
    results_pdb_file = os.path.join('/home/ec2-user/DiffDock', results_pdb_df['protein_path'][0])

    # Remove HETATM records
    results_pdb_file_no_hetatms = remove_hetatms(results_pdb_file)

    # Process results directory
    df_results = process_results_directory(results_directory, results_pdb_file_no_hetatms)

    # Save results
    tsv_file = save_results(df_results, output_directory)
    print(f"Results saved to {tsv_file}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process the results.')
    parser.add_argument('results_dir', type=str, help='Path to the results directory.')
    parser.add_argument('output_dir', type=str, help='Path to the gnina output directory.')
    args = parser.parse_args()
    print('main results directory', args.results_dir)
    print('output directory', args.output_dir)
    main(args.results_dir, args.output_dir)