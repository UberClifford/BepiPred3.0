### IMPORTS ###

import sys 
import argparse
from pathlib import Path

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Takes fasta file with antigens where epitopes are annoated capitilized letters as input. Removes identical sequences.\
        Separates sequences with unknown AA characters and peptides. Combines sequence epitope annotations for sequneces which are identical except for their epitope annotatios.", )
parser.add_argument("-i", required=True, action="store", dest="infile", type=Path, help="Input file. Should be a fasta formatted file with antigens where epitopes are annoated capitilized letters.")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target directory. Will write results to this directory.")
parser.add_argument("-pl", required=False, action="store", dest="peptide_len", type=int, default = 40, help="Peptide length. Entries with sequence length lower than this, will be regarded as a peptide and separated from main results.")
parser.add_argument("-sn", required=False, action="store_true", dest="short_name", default = False, help="Use short name for combined pdb accs (just the first accession in combined accs)")


args = parser.parse_args()
infile = args.infile
target_dir = args.target_dir
peptide_len = args.peptide_len
short_name = args.short_name

if not infile.is_file():
    sys.exit(f"The input file was not a valid file.\nInput file given:{infile.name})")
if not target_dir.is_dir():
    target_dir.mkdir(exist_ok=False, parents=True)

infile = open(infile, "r")
readfile = infile.readlines()
infile.close()

def read_accs_and_sequences_from_fasta(readfile):
    pdb_accs = list()
    sequences = list()
    seq = ""
    read_pdb_acc = False
    
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

def detect_sequence_with_unknowns(pdb_seq):
    lowercase_pdb_seq = pdb_seq.lower()
    seq_with_unknowns = None
    
    if "x" in lowercase_pdb_seq:
        seq_with_unknowns = pdb_seq
    else:
        pass

    return seq_with_unknowns, lowercase_pdb_seq

def detect_peptide(pdb_seq, allowed_protein_len = 40):
    seq_len = len(pdb_seq)
    peptide = None
    if seq_len < allowed_protein_len:
        peptide = pdb_seq
    else:
        pass
    return peptide

def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

def combine_epitope_annotations(similar_accs, lowercase_pdb_seq, use_short_name = False):

    #get sequence length.
    seq_len = len(lowercase_pdb_seq)
    combined_seq = ""
    acc_check_list = list()

    for i in range(seq_len):
        #if there is an epitope annoated in any of the similar accesions
        similar_accs_residue = [similar_acc[1][i] for similar_acc in similar_accs if similar_acc[1][i] == lowercase_pdb_seq[i].upper()]

        if similar_accs_residue:
            combined_seq += lowercase_pdb_seq[i].upper()
        else:
            combined_seq += lowercase_pdb_seq[i]

    pdb_acc_list = [similar_acc[0] for similar_acc in similar_accs]
    if use_short_name:
        combined_pdb_acc = (similar_accs[0][0], combined_seq)
    else:
        combined_acc_names = "|".join(pdb_acc_list)
        combined_pdb_acc = (combined_acc_names, combined_seq)

    return combined_pdb_acc, pdb_acc_list

def check_for_epitopes(pdb_seq):
    """
    Checking that there is at least one epitope annotation.
    """
    pos = False
    for res in pdb_seq:
        pos = res.isupper()
        if pos:
            break
    return pos

def combine_epitope_annotations_and_find_unknowns(pdb_accs, short_name = False, allowed_protein_len = 40):

    """
    Inputs: pdb_accs. List of tuples, [(pdb_acc, seq)...]
    Outputs: combined_pdb_accs. List of tuples, [(pdb_acc, seq)...] Accesions which only differ in epitope annotations (so have different lower and uppercase annotation).
                                are combined. 
                 pdb_accs_with_unknowns. List of tuples. [(pdb_acc, seq)...]. Accesions which contained unknown characters. Marked with x/X.   
    """
    combined_pdb_accs = list()
    pdb_accs_with_unknowns = list()
    peptides = list()
    too_big_sequences = list()
    complete_negatives = list()
    pdb_acc_check_list = list()

    for pdb_acc in pdb_accs:
        pdb_seq = pdb_acc[1]
        seq_with_unknowns, lowercase_pdb_seq = detect_sequence_with_unknowns( pdb_seq )
        peptide = detect_peptide(pdb_seq, allowed_protein_len)
        epitope_exists = check_for_epitopes(pdb_seq)
        too_big_seq = check_maximum_length(pdb_seq)

        #append sequneces with unknown characters
        if seq_with_unknowns != None:
            pdb_accs_with_unknowns.append(pdb_acc)
        #append peptides
        if peptide != None:
            peptides.append(pdb_acc)
        #append complete negatives
        if not epitope_exists:
            complete_negatives.append(pdb_acc)
        #append sequences that are too big for ESM encoding (above 1024 sequence length)
        if too_big_seq != None:
            too_big_sequences.append(pdb_acc)

        if seq_with_unknowns == None and peptide == None and epitope_exists and too_big_seq == None and pdb_acc[0] not in pdb_acc_check_list:
            #finding accesions which identical, except for epitope annotation
            similar_accs = [pdb_acc for pdb_acc in pdb_accs if pdb_acc[1].lower() == lowercase_pdb_seq]
            combined_pdb_acc, pdb_acc_list = combine_epitope_annotations(similar_accs, lowercase_pdb_seq, short_name)
            combined_pdb_accs.append(combined_pdb_acc)
            pdb_acc_check_list.extend(pdb_acc_list)

    return combined_pdb_accs, pdb_accs_with_unknowns, peptides, complete_negatives, too_big_sequences

def check_maximum_length(pdb_seq, maximum_length=1024):
    """
    The ESM can only encode proteins of a maximum length of 1024.
    """

    seq_len = len(pdb_seq)
    too_big_seq = None

    if seq_len > maximum_length:
        too_big_seq = pdb_seq
    else:
        pass

    return too_big_seq

pdb_accs_and_sequences = read_accs_and_sequences_from_fasta(readfile)
combined_pdb_accs, pdb_accs_with_unknowns, peptides, complete_negatives, too_big_sequences = combine_epitope_annotations_and_find_unknowns(pdb_accs_and_sequences, short_name, peptide_len)

#some statistics
print(f"Number of pdb accs in input file were: {len(pdb_accs_and_sequences)}")
print(f"After removing sequences with unknown characters, proteins shorter than {peptide_len}, complete negative, sequences longer than 1024, and combining epitope annotations.\
 There were {len(combined_pdb_accs)} left.") 
print(f"Number of seqeunces with unknown characters: {len(pdb_accs_with_unknowns)}")
print(f"Number of seqeunces which length of {peptide_len} and below: {len(peptides)}")
print(f"Number of seqeuences with length above 1024: {len(too_big_sequences)}")
print(f"Number of seqeunces with no epitopes at all: {len(complete_negatives)}")

#saving results
data_to_fasta_format(combined_pdb_accs, target_dir /  "combined_unique_antigen_epitopes.fasta")
data_to_fasta_format(pdb_accs_with_unknowns, target_dir / "antigen_epitopes_with_unknown_chars.fasta")
data_to_fasta_format(peptides, target_dir / "peptide_antigen_epitopes.fasta")
data_to_fasta_format(too_big_sequences, target_dir / "long_antigen_epitopes.fasta")
data_to_fasta_format(complete_negatives, target_dir / "complete_negatives_antigen_epitopes.fasta")
