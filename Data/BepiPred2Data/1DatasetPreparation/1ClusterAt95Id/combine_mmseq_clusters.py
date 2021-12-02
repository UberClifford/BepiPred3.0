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

def pairwise_alignment_on_cluster(cluster_sequences):
    """
    Inputs: List of tuples: [(cluter_acc, cluster_seq)...]
    Output: Aligned sequence
    """

    #sort list from longest to shortest sequence is first.
    cluster_size = len(cluster_sequences)
    if cluster_size == 1:
        cluster_representative = cluster_sequences[0]
    else:
        sorted_cluster_sequences = sorted(cluster_sequences, key=lambda tup: len(tup[1]), reverse=True)
        template_acc = sorted_cluster_sequences[0][0]
        template_seq = sorted_cluster_sequences[0][1]

        for i in range(1, cluster_size):
            clus_seq = sorted_cluster_sequences[i][1]
            #pairwise alignment on non-annotated sequences
            alignment = pairwise2.align.globalxx(template_seq.upper(), clus_seq.upper() )[0]
            template_alignment = alignment.seqA
            paired_seq_alignment = alignment.seqB
            seq_len = len(template_alignment)
            
            #put annotations back into alignment
            annotated_template_alignment = str()
            annotated_paired_seq_alignment = str()
            j = 0
            k = 0
            for i in range(seq_len):
                
                if template_alignment[i] == "-":
                    annotated_template_alignment += "-"
                else:
                    annotated_template_alignment += template_seq[j]
                    j += 1
                
                if paired_seq_alignment[i] == "-":
                    annotated_paired_seq_alignment += "-"
                else:
                    annotated_paired_seq_alignment += clus_seq[k]
                    k += 1

            new_template_seq = str()

            for i in range(seq_len):
                template_res = annotated_template_alignment[i]
                paired_seq_res = annotated_paired_seq_alignment[i]
 
                #if both sequences agree, (A/a = A/a) 
                if template_res == paired_seq_res:
                    new_template_seq += template_res
                #if gap at one sequence and residue in other sequence, (A/a , -) or (- , A/a) 
                elif template_res != "-" and paired_seq_res == "-":
                    new_template_seq += template_res
                elif template_res == "-" and paired_seq_res != "-":
                    new_template_seq += paired_seq_res
                #if epitope annotation in paired sequence and not current template
                elif template_res == paired_seq_res.lower():
                    new_template_seq += template_res.upper()
                #if epitope annotation in template sequence and not in paired sequence
                elif template_res.lower() == paired_seq_res:
                    new_template_seq += template_res
                #misalignment
                elif paired_seq_res.lower() != template_res.lower():
                    print(f"Misalignment: {template_res, paired_seq_res}")
                else:
                    print("something else occured")
                    print(template_res, paired_seq_res)

            #updated template sequence
            template_seq = new_template_seq

        cluster_representative = (template_acc, template_seq)

    return cluster_representative

def get_cluster_accs_from_mmseq_output(readfile):
    
    cluster_accs = list()
    cluster_sequences = list()
    cluster_representatives = list()
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
                    cluster_rep = pairwise_alignment_on_cluster(cluster_accs_and_seqs)
                    cluster_representatives.append(cluster_rep)
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
        print(cluster_accs_and_seqs)
    #if last cluster has multiple sequences
    else:
        cluster_accs_and_seqs = list( zip(cluster_accs, cluster_sequences) )
    cluster_rep = pairwise_alignment_on_cluster(cluster_accs_and_seqs)
    cluster_representatives.append(cluster_rep)

    return cluster_representatives, total_acc_count

def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

cluster_representatives, total_num_accs = get_cluster_accs_from_mmseq_output(readfile)
#some statistics
print(f"Total number of pdb accs in clusters: {total_num_accs}")
print(f"Total number of sequences after combining clusters: {len(cluster_representatives)}")
#saving results
data_to_fasta_format(cluster_representatives, target_dir /  "clustered_combined_antigens.fasta")
