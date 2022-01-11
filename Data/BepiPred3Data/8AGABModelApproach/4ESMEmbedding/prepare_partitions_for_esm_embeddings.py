### IMPORTS ###
import pickle
from pathlib import Path
import sys
import subprocess
import argparse

#Expect a directory with pickle files containing the antibody antigen interactions.

### STATIC PATHS AND VARIABLES ###
ROOT_DIR = Path.cwd()

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Takes pickle formated files, single chain ab-ag complexes, converts to probable light and heavy chain. And joins with ag-ab light and heavy chain separated files.")
parser.add_argument("-i", required=True, action="store", dest="in_dir", type=Path, help="Input directory containing antibody-antigen interactions. Should be pickle formatted files")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target directory. Will write results to this directory.")

args = parser.parse_args()
in_dir = args.in_dir
target_dir = args.target_dir
pickle_files = in_dir.glob("*.pickle*")

if not target_dir.is_dir():
    target_dir.mkdir(exist_ok=False, parents=True)

### FUNCTIONS ###
def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)


def get_light_ag_and_heavy_chains(ag_ab_interactions):
    
    heavy_chains = list()
    light_chains = list()
    antigen_chains = list()

    for ag_ab_interaction in ag_ab_interactions:
        light_chain_id = ag_ab_interaction[0]
        light_chain = ag_ab_interaction[1]
        heavy_chain_id = ag_ab_interaction[2]
        heavy_chain = ag_ab_interaction[3]
        ag_acc = ag_ab_interaction[4]
        ag_chain = ag_ab_interaction[5]
        unique_ag_ab_complex_id = ag_ab_interaction[6]


        light_chains.append( (unique_ag_ab_complex_id, light_chain.upper()) )
        heavy_chains.append( (unique_ag_ab_complex_id, heavy_chain.upper()) )
        antigen_chains.append( (unique_ag_ab_complex_id, ag_chain.upper()) )

    return light_chains, heavy_chains, antigen_chains

### MAIN ###

for pickle_file in pickle_files:

    partition_target_dir = str(pickle_file.name).split(".pickle")[0]
    partition_target_dir = target_dir / partition_target_dir
    
    if not partition_target_dir.is_dir():
        partition_target_dir.mkdir(exist_ok=False, parents=True)
        for d_name in ("ESM_heavy_chains", "ESM_light_chains", "ESM_antigens"):
            esm_dir = partition_target_dir / d_name
            esm_dir.mkdir(exist_ok=False, parents=True)

    ag_ab_interactions = pickle.load(open(pickle_file, "rb"))
    light_chains, heavy_chains, antigen_chains = get_light_ag_and_heavy_chains(ag_ab_interactions)
    data_to_fasta_format(light_chains, partition_target_dir / "light_chains.fasta")
    data_to_fasta_format(heavy_chains, partition_target_dir / "heavy_chains.fasta")
    data_to_fasta_format(antigen_chains, partition_target_dir / "antigens.fasta")
