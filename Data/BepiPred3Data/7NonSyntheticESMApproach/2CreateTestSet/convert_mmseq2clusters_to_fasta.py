### IMPORTS ###
import sys 
import argparse
from pathlib import Path
from Bio import pairwise2

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Takes fasta file with antigens where epitopes are annoated capitilized letters as input. Removes identical sequences.\
        Separates sequences with unknown AA characters and peptides. Combines sequence epitope annotations for sequneces which are identical except for their epitope annotatios.", )
parser.add_argument("-i", required=True, action="store", dest="infile", type=Path, help="Input file. Should be a fasta formatted file with antigens where epitopes are annoated capitilized letters.")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target directory. Will write results to this directory.")

args = parser.parse_args()
infile = args.infile
target_dir = args.target_dir

if not infile.is_file():
    sys.exit(f"The input file was not a valid file.\nInput file given:{infile.name})")
if not target_dir.is_dir():
    target_dir.mkdir(exist_ok=False, parents=True)

infile = open(infile, "r")
readfile = infile.readlines()
infile.close()

def get_cluster_accs_from_mmseq_output(readfile):
    
    cluster_accs = list()
    cluster_sequences = list()
    all_cluster_acc_and_seqs = list()
    cluster_sequence = str()
    past_first_cluster = False
    read_seq = False
    i = 0
    total_acc_count = 0

    for line in readfile:
        line = line.strip()
        if line.startswith(">"):
            total_acc_count += 1
            pdb_acc = line.split(">")[1]
            
            #at new cluster
            if cluster_accs and pdb_acc == cluster_accs[-1]:
                total_acc_count -= 1

                if past_first_cluster:
                    cluster_accs = cluster_accs[:-1]
                    #run pairwise alignments and get cluster representative
                    cluster_accs_and_seqs = list( zip(cluster_accs, cluster_sequences) )
                    all_cluster_acc_and_seqs.extend(cluster_accs_and_seqs)
                    #reset cluster
                    cluster_accs = list()
                    cluster_sequences = list()
                    cluster_accs.append(pdb_acc)

                else:
                    past_first_cluster = True

            else:
                cluster_accs.append(pdb_acc)
                if past_first_cluster:
                    cluster_sequences.append(cluster_sequence)
                    cluster_sequence = str()
        
        #add to cluster sequence string
        else:
            cluster_sequence += line

    #run pairwise alignment on last sequence.
    #if last cluster only has one sequence
    if len(cluster_accs) == 1:
        cluster_accs_and_seqs = [(cluster_accs[0], cluster_sequence) ]
    #if last cluster has multiple sequences
    else:
    	#append last sequence 
        cluster_sequences.append(cluster_sequence)
        cluster_accs_and_seqs = list( zip(cluster_accs, cluster_sequences) )

    all_cluster_acc_and_seqs.extend(cluster_accs_and_seqs)

    return all_cluster_acc_and_seqs, total_acc_count

def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

all_cluster_acc_and_seq, total_acc_count = get_cluster_accs_from_mmseq_output(readfile)
data_to_fasta_format(all_cluster_acc_and_seq, target_dir /  "bepipred3_test_set_removed.fasta")
