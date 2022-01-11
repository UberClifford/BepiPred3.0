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
    all_cluster_acc_and_seqs = dict()
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
                    all_cluster_acc_and_seqs[cluster_accs[0]] = list( zip(cluster_accs, cluster_sequences) )
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
        all_cluster_acc_and_seqs[cluster_accs[0]] = [ (cluster_accs[0], cluster_sequence) ]
    #if last cluster has multiple sequences
    else:
        #append last sequence 
        cluster_sequences.append(cluster_sequence)
        all_cluster_acc_and_seqs[cluster_accs[0]] = list( zip(cluster_accs, cluster_sequences) )

    return all_cluster_acc_and_seqs, total_acc_count

def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

def get_sorted_cluster_sizes(all_cluster_acc_and_seq):
    cluster_sizes = list()
    for cluster_rep in all_cluster_acc_and_seq.keys():
        cluster_sizes.append( (cluster_rep, len(all_cluster_acc_and_seq[cluster_rep])) )
    cluster_sizes = sorted(cluster_sizes, key=lambda pair: pair[1], reverse=True)
    
    return cluster_sizes

def create_5_partitions(all_cluster_acc_and_seq, cluster_sizes, target_dir):
    
    i = 1
    fold1 = str()
    fold2 = str()
    fold3 = str()
    fold4 = str()
    fold5 = str()

    for cluster in cluster_sizes:
        #cluster = (cluster_rep, cluster_size)
        cluster_rep = cluster[0]
        
        cluster_acc_and_seqs = all_cluster_acc_and_seq[cluster_rep]
        output = str()
        for cluster_acc_and_seq in cluster_acc_and_seqs:
            #distinguish cluster representatives using mmseq2 format
            acc_name = cluster_acc_and_seq[0]
            seq = cluster_acc_and_seq[1]
            if acc_name == cluster_rep:
                output += f">{acc_name}\n>{acc_name}\n{seq}\n"
            else:
                output += f">{acc_name}\n{seq}\n"

        if i == 1:
            fold1 += output
        elif i == 2:
            fold2 += output
        elif i == 3:
            fold3 += output
        elif i == 4:
            fold4 += output
        elif i == 5:
            fold5 += output

        i += 1
        if i == 6:
            i = 1

    fold1 = fold1[:-1]
    fold2 = fold2[:-1]
    fold3 = fold3[:-1]
    fold4 = fold4[:-1]
    fold5 = fold5[:-1]

    #write to fasta files
    with open(target_dir / "Fold1.fasta", "w") as outfile:
        outfile.write(fold1)
    with open(target_dir / "Fold2.fasta", "w") as outfile:
        outfile.write(fold2)
    with open(target_dir / "Fold3.fasta", "w") as outfile:
        outfile.write(fold3)
    with open(target_dir / "Fold4.fasta", "w") as outfile:
        outfile.write(fold4)
    with open(target_dir / "Fold5.fasta", "w") as outfile:
        outfile.write(fold5)

all_cluster_acc_and_seq, total_acc_count = get_cluster_accs_from_mmseq_output(readfile)
cluster_sizes = get_sorted_cluster_sizes(all_cluster_acc_and_seq)
create_5_partitions(all_cluster_acc_and_seq, cluster_sizes, target_dir)
