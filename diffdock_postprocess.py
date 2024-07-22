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
        print(f"An shell command error occurred: {str(e)}")
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


def process_results_directory(results_directory, results_pdb_file_no_hetatms, pdb_filepath):
    """Processes each result directory and returns a DataFrame with the results."""
    df_rows = []
    results_dirs = glob(os.path.join(results_directory, "complex_0"))
    print('DiffDock folders containing sdf predictions:', results_dirs)

    for results_dir in tqdm(results_dirs, desc="runs"):
        results_sdfs = [os.path.join(results_dir, f) for f in os.listdir(results_dir) if "confidence" in f and f.endswith(".sdf")]
        print('DiffDock folder:', results_sdfs)
        
        for results_sdf in tqdm(results_sdfs, leave=False, desc="files"):
            filename = results_sdf.split("/")[-1]
            print("SDF filename:", filename)
            process_sdf_file(results_sdf, results_pdb_file_no_hetatms, df_rows, pdb_filepath)

    return pd.DataFrame(df_rows, columns=["sdf_file", "smile", "diffdock_confidence", "gnina_scored_affinity", "gnina_minimized_affinity", "pdb_file"])


def process_sdf_file(results_sdf, results_pdb_file_no_hetatms, df_rows, pdb_filepath):
    """Processes an SDF file and appends the results to the provided list."""
    smile = sdf_to_smile(results_sdf)	
    confidence = re.findall("confidence([\-\.\d]+)\.sdf", results_sdf)[0]
    scored_affinity = get_affinity(results_pdb_file_no_hetatms, results_sdf, score_only=True)
    minimized_affinity = get_affinity(results_pdb_file_no_hetatms, results_sdf, local_only=True)
    sdf_filename = results_sdf.split("/")[-1]
    pdb_filename = pdb_filepath.split("/")[-1]
    print('PDB Filepath', pdb_filepath)
    df_rows.append((sdf_filename, smile, float(confidence), scored_affinity, minimized_affinity, pdb_filename))


def get_affinity(pdb_file, sdf_file, score_only=False, local_only=False):
    """Gets the affinity using gnina."""
    cmd = f"/home/ec2-user/gnina --{'score_only' if score_only else 'local_only --minimize'} -r '{pdb_file}' -l '{sdf_file}'"
    if local_only:
        cmd += f" --autobox_ligand '{sdf_file}' --autobox_add 2"
    stdout = run_command(cmd)
    return float(re.findall("Affinity:\s*([\-\.\d+]+)", stdout)[0])


def save_results(df_results, out_dir, results_dir, filename):
    """Saves the results DataFrame to a TSV file."""
    # tsv_file = os.path.join(output_directory, "df_diffdock_results.tsv")

    # Construct the absolute output directory path
    output_directory = os.path.abspath(os.path.join(results_dir, out_dir))
    print("TSV filename", filename)
    print("Output directory", output_directory)

    tsv_file = os.path.join(output_directory, filename)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    print("Results dataframe", df_results)
    df_results.to_csv(tsv_file, sep='\t', index=None)
    return tsv_file

def main(results_directory, output_directory, protein_ligand_csv, filename):
    # Load PDB file path
    project_dir = os.getcwd()
    
    csv_path = os.path.join(project_dir, protein_ligand_csv)
    protein_ligand_df = pd.read_csv(csv_path)
    pdb_filepath = os.path.join(project_dir, protein_ligand_df['protein_path'][0])
    results_directory = os.path.join(project_dir, results_directory)

    print('main function:')
    print('project_dir', project_dir)
    print('pdb_filepathh', pdb_filepath)
    print('csv_path', csv_path)
    print('output_directory', output_directory)
    print('=====')

    # Remove HETATM records
    results_pdb_file_no_hetatms = remove_hetatms(pdb_filepath)

    # Process results directory
    df_results = process_results_directory(results_directory, results_pdb_file_no_hetatms, pdb_filepath)

    # Save results
    tsv_file = save_results(df_results, output_directory, results_directory, filename)
    print(f"Results saved to {tsv_file}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process the results.')
    parser.add_argument('results_dir', type=str, help='Path to the results directory.')
    parser.add_argument('output_dir', type=str, help='Path to the gnina output directory.')
    parser.add_argument('protein_ligand_csv', type=str, help='Name of protein liagand batch csv.')
    parser.add_argument('filename', type=str, help='Name of summary tsv file.')
    args = parser.parse_args()
    print('DiffDock inference results directory:', args.results_dir)
    print('Gnina results output directory:', args.output_dir)
    print('Name of protein liagand batch csv:', args.protein_ligand_csv)
    print('filename:', args.filename)
    main(args.results_dir, args.output_dir, args.protein_ligand_csv, args.filename)
    
