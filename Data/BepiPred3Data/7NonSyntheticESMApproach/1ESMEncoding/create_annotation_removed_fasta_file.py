### IMPORTS ###

import argparse
from pathlib import Path
import sys

parser = argparse.ArgumentParser("Takes fasta file with antigens where epitopes are annoated capitilized letters as input. Removes annotations and capitilizes everythings. Saves to new fasta file", )
parser.add_argument("-i", required=True, action="store", dest="infile", type=Path, help="Input file. Should be a fasta formatted file with antigens where epitopes are capitilized.")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target file.")

args = parser.parse_args()
infile = args.infile
target_dir = args.target_dir

if not infile.is_file():
    sys.exit(f"The input file was not a valid file.\nInput file given:{infile.name})")

if not target_dir.parent.is_dir():
    target_dir.parent.mkdir(exist_ok=False, parents=True)

#infile = open(infile, "r")
#readfile = infile.readlines()
#infile.close()

### FUNCTIONS ###
def read_accs_and_sequences_from_fasta(infile):
    """
    Input: readfile: Fasta file. 
    Outputs: pdb_accs_and_sequences: List of tuples. Containing accs and sequences, e.g. [(acc, aTHNtem..)..()]. 
    """
    pdb_accs = list()
    sequences = list()
    seq = ""
    read_pdb_acc = False
    
    if not infile.is_file():
        print(f"The input file was invalid. Invalid file was {infile}")
        
    infile = open(infile, "r")
    readfile = infile.readlines()
    infile.close()

    for line in readfile:
        line = line.strip()
        if line.startswith(">"):
            pdb_acc = line.split(">")[1]
            if read_pdb_acc:
                pdb_accs.append(pdb_acc)
                sequences.append(seq)
                #reset sequence string
                seq = ""
            #catch first pdb accesion. First pdb acc is unique
            else:
                pdb_accs.append(pdb_acc)
        else:
            seq += line
            read_pdb_acc = True

    #get last sequence
    sequences.append(seq)
    pdb_accs_and_sequences = tuple( zip(pdb_accs, sequences) )
    return pdb_accs_and_sequences

def create_annotation_removed_pdb_accs(readfile):
    """
    Inputs: readfile: Valid fasta file path containing antigens
                      where eptiope are annotated with capitlized letters.         
    Outputs: not_annotated_pdb_accs: List of tuples, containing accesions and non-annotated seuqences e.g. [(acc, ATAGAH..)..()]  
    """
    
    not_annotated_pdb_accs_and_sequences = list()
    pdb_accs_and_sequences = read_accs_and_sequences_from_fasta(readfile)
    
    #convert all sequences to uppercase
    for pdb_acc_and_sequence in pdb_accs_and_sequences:
        pdb_acc = pdb_acc_and_sequence[0]
        pdb_sequence = pdb_acc_and_sequence[1].upper()
        not_annotated_pdb_accs_and_sequences.append( (pdb_acc, pdb_sequence) ) 
    
    return not_annotated_pdb_accs_and_sequences

def data_to_fasta_format(pdb_accs, outfile_path):
    """
    Inputs: pdb_accs: List of tuples, consisting of accesions and sequences [(accs, ATGRE..)..]
    Outputs: Fasta file at outfile_path
    """
    outfile_path = Path(outfile_path)
#    try:
#        outfile_path.mkdir(parents=True, exist_ok=False)
#    except FileExistsError:
#        print("Save directory was already there. Saving data there.")
#    else:
#        print("Save directory not found. Made new one. ")
        
    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

not_annotated_pdb_accs_and_sequences = create_annotation_removed_pdb_accs(infile)
data_to_fasta_format(not_annotated_pdb_accs_and_sequences, target_dir)
