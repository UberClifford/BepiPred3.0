### IMPORTS ###
from pathlib import Path
import numpy as np
import sys
import torch
import pickle
import matplotlib.pyplot as plt
import numpy as np
import random
import argparse
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline

### STATICS PATHS ###
ROOT_DIR = Path.cwd()
ESM_EMBEDDINGS = ROOT_DIR / "../1ESMEncoding/ESMEncodedProteins/"
PARTITION_PATH = ROOT_DIR / "../4DivideProteinsInto5Partitions/Results/"
TEST_DATA_PATH = ROOT_DIR / "../2CreateTestSet" / "bepipred3_test.fasta"
CROSS_VALIDATION_SAVE_PATH = ROOT_DIR / "../6CrossValidationNonSynthESM"

ESM_EMBEDDINGS = ESM_EMBEDDINGS.resolve()
PARTITION_PATH = PARTITION_PATH.resolve()
TEST_DATA_PATH = TEST_DATA_PATH.resolve()
CROSS_VALIDATION_SAVE_PATH = CROSS_VALIDATION_SAVE_PATH.resolve()

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

def get_cluster_accs_from_mmseq_output(infile):

    """
    readfile:mmseq formatted file
    Outputs: Dict {cluster_rep:[(cluster_rep, seq)..(cluster_acc, seq)..],
                   cluster_rep:[(cluster_rep, seq)..(cluster_acc, seq)..]...}
    """
    
    cluster_accs = list()
    cluster_sequences = list()
    all_cluster_acc_and_seqs = dict()
    cluster_sequence = str()
    past_first_cluster = False
    read_seq = False
    i = 0
    total_acc_count = 0

    if not infile.is_file():
        print(f"The input file was invalid. Invalid file was {infile}")

    infile = open(infile, "r")
    readfile = infile.readlines()
    infile.close()

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
        all_cluster_acc_and_seqs[cluster_accs[0]] = list( zip(cluster_accs, cluster_sequences) )

    return all_cluster_acc_and_seqs, total_acc_count

def get_epitope_labels_as_binary_array_from_sequence(sequence):
    """
    Input: sequence: String. Amino acid sequence where epitope are marked with capitilized letters. 
    Outputs: binary_epitope_annotated_array: Np array. Binary array where 0 = Non-epitope and 1 = Epitope. 
    """
    
    seq_len = len(sequence)
    binary_epitope_annotated_array = np.zeros(seq_len, dtype= np.uint8)
    for i in range (seq_len):
        #if it is an epitope
        if sequence[i].isupper():
            binary_epitope_annotated_array[i] = 1
        #if it is not an epitope
        elif sequence[i].islower():
            binary_epitope_annotated_array[i] = 0
        #if something else, it's weird.
        else:
            print(f"Something weird in sequence detected: {sequence[i]}")
            break
    
    return binary_epitope_annotated_array

def get_labels_dict(cluster_acc_and_seqs, esm_encoded_data_path):
    """
    Output: Dict {cluster_acc: (binary_annotations, string_seq_anno, esm representation)...}
    """

    label_dict = dict()
    for cluster_acc_and_seq in cluster_acc_and_seqs:

        cluster_acc = cluster_acc_and_seq[0]
        cluster_seq = cluster_acc_and_seq[1]

        #load esm embedding
        esm_encoded_acc_path = esm_encoded_data_path /  f"{cluster_acc}.pt"
        esm_encoded_acc = torch.load(esm_encoded_acc_path)
        esm_representations = esm_encoded_acc["representations"][33]

        binary_seq = get_epitope_labels_as_binary_array_from_sequence(cluster_seq)
        label_dict[cluster_acc] = (binary_seq, cluster_seq, esm_representations)

    return label_dict

def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

def put_annotation_back_into_msa(msa_accs_and_seqs, label_dict):
    
    num_seqs = len(msa_accs_and_seqs)
    msa_len = len(msa_accs_and_seqs[0][1])
    annotated_msa_accs_and_seqs = np.zeros((num_seqs, msa_len+1), dtype=object)

#   #put annotations back into muscle multiple sequence alignments
#   #also pad esm embeddings and binary labels at gaps. 
#   #if gap, binary label = 2 and esm_embedding torch.zeros(1280)
    for j in range(num_seqs):
        i = 0
        annotated_seq = str()
        msa_acc = msa_accs_and_seqs[j][0]
        msa_seq = msa_accs_and_seqs[j][1]
        binary_annotations = label_dict[msa_acc][0]
        annotations = label_dict[msa_acc][1]
        seq_len = len(annotations)
        esm_embeddings = label_dict[msa_acc][2]
        padded_esm_embeddings = torch.zeros((msa_len, 1280))
        padded_binary_annotations = np.zeros(msa_len)

        for k in range(msa_len):
            if msa_seq[k] == "-":
                annotated_seq += "-"
                padded_esm_embeddings[k] = torch.zeros(1280)
                padded_binary_annotations[k] = 2

            #if residue is part of an epitope
            elif msa_seq[k] == annotations[i]:
                annotated_seq += annotations[i]
                padded_esm_embeddings[k] = esm_embeddings[i]
                padded_binary_annotations[k] = binary_annotations[i]
                i += 1
            #if residue is not part of an epitope
            elif msa_seq[k].lower() == annotations[i]:
                annotated_seq += annotations[i]
                padded_esm_embeddings[k] = esm_embeddings[i]
                padded_binary_annotations[k] = binary_annotations[i]
                i += 1
            else:
                print("Something is wrong!")


        #update dict with padded binary annotations and ESM representations
        label_dict[msa_acc] = (padded_binary_annotations, annotated_seq, padded_esm_embeddings, seq_len) 

        row = [msa_acc]
        row.extend(list(annotated_seq))
        annotated_msa_accs_and_seqs[j] = row

    return annotated_msa_accs_and_seqs, label_dict

def save_data(X, y, save_dirname, filename = "foo"):
    
    try:
        save_dirname.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print("Save directory was already there. Saving data there.")
    else:
        print("Save directory not found. Made new one. ")
    
    save_path = save_dirname / filename
    np.savez_compressed(save_path, saved_X = X, saved_y = y)


def create_crossvalidation_fold(X_train, y_train, X_val, y_val, save_path):
    np.savez_compressed(save_path, saved_X_train=X_train, saved_X_val=X_val,
                                        saved_y_train=y_train, saved_y_val=y_val)

def get_from_global_pairwise_alignment(all_cluster_acc_and_seqs):

    cluster_reps = set(all_cluster_acc_and_seqs.keys())
    positional_esm_embeddings = list()
    binary_labels = list()

    for cluster_rep in cluster_reps:
        cluster_acc_and_seqs = all_cluster_acc_and_seqs[cluster_rep]

        #write to temporary file 
        data_to_fasta_format(cluster_acc_and_seqs, "cluster_accs_temp.fasta")
        #get epitope annotations, binary sequence labels and ESM embeddings as dictionary
        label_dict = get_labels_dict(cluster_acc_and_seqs, ESM_EMBEDDINGS)
#       #do global multiple sequence alignment (Muscle is not case sensitive)
        muscle_cline = MuscleCommandline(input="cluster_accs_temp.fasta", 
                                 out="cluster_accs_temp_alignment.fasta",  
                                 maxiters = 1)
        muscle_cline()
 
#       #read multiple sequence alignments from file
        msa_accs_and_seqs = read_accs_and_sequences_from_fasta(Path("cluster_accs_temp_alignment.fasta"))
        #put annotations back into MSA alignments. Also pad binary labels and ESM embeddings all be length of MSA. 
        annotated_msa_accs_and_seqs, label_dict = put_annotation_back_into_msa(msa_accs_and_seqs, label_dict)
        
        msa_acc_names = annotated_msa_accs_and_seqs[:,0]
        msa_len = np.shape(annotated_msa_accs_and_seqs)[1]
      
        for pos in range(1, msa_len):
            position_values = list(annotated_msa_accs_and_seqs[:, pos])
            #check if there is an epitope
            pos_epitopes_idx = [idx for idx, value in enumerate(position_values) if value.isupper()]
            #if positive epitope uncovered
            if pos_epitopes_idx:
                #randomly pick acession msa, which had a epitope annotation 
                acc_idx = random.choice(pos_epitopes_idx)
            #if no positive epitopes at that position at all
            else:
                #randomly pick accesion in msa, that is not a gap.
                acc_idx = random.choice([idx for idx, value in enumerate(position_values) if value != "-"])

            #check for error
            if label_dict[msa_acc_names[acc_idx]][0][pos-1] == 2:
                print("Something is wrong. Padding is going into ESM embedding")
  
            #take esm embedding from corresponding accesion and position
            esm_embedding = label_dict[msa_acc_names[acc_idx]][2][pos-1]
            binary_label = label_dict[msa_acc_names[acc_idx]][0][pos-1]
            seq_len = torch.tensor([label_dict[msa_acc_names[acc_idx]][3]])
            esm_embedding = torch.cat((esm_embedding, seq_len))
  
            positional_esm_embeddings.append(esm_embedding)
            binary_labels.append(binary_label)


    #convert to tensor (num_of_positional_embeddings, 1281)
    positional_esm_embeddings = torch.stack(positional_esm_embeddings).numpy()
    binary_labels = np.asarray(binary_labels)

    return positional_esm_embeddings, binary_labels 

### MAIN ###

#Fold1
all_cluster_acc_and_seqs, _ = get_cluster_accs_from_mmseq_output(PARTITION_PATH / "Fold1.fasta")
partition_1_esm_rep, partition_1_labels = get_from_global_pairwise_alignment(all_cluster_acc_and_seqs)
#Fold2
all_cluster_acc_and_seqs, _ = get_cluster_accs_from_mmseq_output(PARTITION_PATH / "Fold2.fasta")
partition_2_esm_rep, partition_2_labels = get_from_global_pairwise_alignment(all_cluster_acc_and_seqs)
##Fold3
all_cluster_acc_and_seqs, _ = get_cluster_accs_from_mmseq_output(PARTITION_PATH / "Fold3.fasta")
partition_3_esm_rep, partition_3_labels = get_from_global_pairwise_alignment(all_cluster_acc_and_seqs)
#Fold4
all_cluster_acc_and_seqs, _ = get_cluster_accs_from_mmseq_output(PARTITION_PATH / "Fold4.fasta")
partition_4_esm_rep, partition_4_labels = get_from_global_pairwise_alignment(all_cluster_acc_and_seqs)
#Fold5
all_cluster_acc_and_seqs, _ = get_cluster_accs_from_mmseq_output(PARTITION_PATH / "Fold5.fasta")
partition_5_esm_rep, partition_5_labels = get_from_global_pairwise_alignment(all_cluster_acc_and_seqs)
#External Test data
all_cluster_acc_and_seqs, _ = get_cluster_accs_from_mmseq_output(TEST_DATA_PATH)
test_esm_rep, test_labels = get_from_global_pairwise_alignment(all_cluster_acc_and_seqs)

print(np.shape(partition_1_esm_rep), np.shape(partition_1_labels))
print(np.shape(partition_2_esm_rep), np.shape(partition_2_labels))
print(np.shape(partition_3_esm_rep), np.shape(partition_3_labels))
print(np.shape(partition_4_esm_rep), np.shape(partition_4_labels))
print(np.shape(partition_5_esm_rep), np.shape(partition_5_labels))
print(np.shape(test_esm_rep), np.shape(test_labels))

### CREATE CROSSVALIDATION DATASETS ###

#Fold1 
val_embeddings = partition_1_esm_rep
val_labels = partition_1_labels
train_embeddings = np.concatenate( (partition_2_esm_rep, partition_3_esm_rep, partition_4_esm_rep, partition_5_esm_rep) )
train_labels = np.concatenate( (partition_2_labels, partition_3_labels, partition_4_labels, partition_5_labels) )
create_crossvalidation_fold(train_embeddings, train_labels, val_embeddings, val_labels, CROSS_VALIDATION_SAVE_PATH / "Fold1")
##Fold2
val_embeddings = partition_2_esm_rep
val_labels = partition_2_labels
train_embeddings = np.concatenate( (partition_1_esm_rep, partition_3_esm_rep, partition_4_esm_rep, partition_5_esm_rep) ) 
train_labels = np.concatenate( (partition_1_labels, partition_3_labels, partition_4_labels, partition_5_labels) )
create_crossvalidation_fold(train_embeddings, train_labels, val_embeddings, val_labels, CROSS_VALIDATION_SAVE_PATH / "Fold2")
#Fold3
val_embeddings = partition_3_esm_rep
val_labels = partition_3_labels
train_embeddings = np.concatenate( (partition_1_esm_rep, partition_2_esm_rep, partition_4_esm_rep, partition_5_esm_rep) )
train_labels = np.concatenate( (partition_1_labels, partition_2_labels, partition_4_labels, partition_5_labels) )
create_crossvalidation_fold(train_embeddings, train_labels, val_embeddings, val_labels, CROSS_VALIDATION_SAVE_PATH / "Fold3")
#Fold4
val_embeddings = partition_4_esm_rep
val_labels = partition_4_labels
train_embeddings = np.concatenate( (partition_1_esm_rep, partition_3_esm_rep, partition_2_esm_rep, partition_5_esm_rep) )
train_labels = np.concatenate( (partition_1_labels, partition_3_labels, partition_2_labels, partition_5_labels) )
create_crossvalidation_fold(train_embeddings, train_labels, val_embeddings, val_labels, CROSS_VALIDATION_SAVE_PATH / "Fold4")
#Fold5
val_embeddings = partition_5_esm_rep
val_labels = partition_5_labels
train_embeddings = np.concatenate( (partition_1_esm_rep, partition_3_esm_rep, partition_4_esm_rep, partition_2_esm_rep) )
train_labels = np.concatenate( (partition_1_labels, partition_3_labels, partition_4_labels, partition_2_labels) )
create_crossvalidation_fold(train_embeddings, train_labels, val_embeddings, val_labels, CROSS_VALIDATION_SAVE_PATH / "Fold5")
#save test data
save_data(test_esm_rep, test_labels, CROSS_VALIDATION_SAVE_PATH, "test")
