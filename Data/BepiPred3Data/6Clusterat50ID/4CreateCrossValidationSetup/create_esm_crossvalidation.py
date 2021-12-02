### IMPORTS ###
from pathlib import Path
import numpy as np
import sys
import torch
from sklearn.model_selection import KFold
import pickle
import matplotlib.pyplot as plt
import numpy as np
import random

### STATICS PATHS ###
ROOT_DIR = Path.cwd()
ESM_EMBEDDINGS = ROOT_DIR / "../3ESMEncoding/"
EPITOPE_ANNOTATIONS = ROOT_DIR / "../2CreateTestSet" / "bepipred3_test_set_removed.fasta" 
CROSS_VALIDATION_SAVE_PATH = ROOT_DIR / "../5CrossValidationClusterAt50Id"
TEST_EPITOPE_ANNOTATIONS = ROOT_DIR / "../2CreateTestSet" / "bepipred3_test.fasta"
ESM_EMBEDDINGS = ESM_EMBEDDINGS.resolve()
EPITOPE_ANNOTATIONS = EPITOPE_ANNOTATIONS.resolve()
CROSS_VALIDATION_SAVE_PATH = CROSS_VALIDATION_SAVE_PATH.resolve()

### FUNCTIONS ###
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

def get_epitope_labels_as_binary_array_from_fasta(fasta_file_path):
    """
    Inputs: fasta_file_path: Valid fasta file path containing antigens where eptiope are annotated with capitlized letters.
    Outputs: binary_epitope_annotated_arrays: Dict, containing accesions and binarized epitope annotation, e.g.
             {acc:00011110, acc1:01001010,}
    
    """
    
    fasta_file_path = Path(fasta_file_path)
    if not fasta_file_path.is_file():
        sys.exit(f"Invalid fasta file given. Fasta file given: {fasta_file_path}")
        
    pdb_accs_and_sequences = read_accs_and_sequences_from_fasta(fasta_file_path)
    binary_epitope_annotated_arrays = dict()
    
    for pdb_acc_and_sequence in pdb_accs_and_sequences:
        pdb_acc = pdb_acc_and_sequence[0]
        pdb_sequence = pdb_acc_and_sequence[1]
        binary_epitope_annotated_array = get_epitope_labels_as_binary_array_from_sequence(pdb_sequence)
        binary_epitope_annotated_arrays[pdb_acc] = binary_epitope_annotated_array
    
    return binary_epitope_annotated_arrays

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

def save_data(X, y, save_dirname, filename = "foo"):
    
    try:
        save_dirname.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print("Save directory was already there. Saving data there.")
    else:
        print("Save directory not found. Made new one. ")
    
    save_path = save_dirname / filename
    np.savez_compressed(save_path, saved_X = X, saved_y = y)

def get_esm_representations_and_labels(esm_encoded_data_path, antigen_with_epitopes_path):
    
    all_binary_labels = list()
    all_esm_representations = list()
    #get labels from annotated data
    binary_annotations = get_epitope_labels_as_binary_array_from_fasta(antigen_with_epitopes_path)
    esm_encoded_acc_paths = list( esm_encoded_data_path.glob("*.pt") )
    for esm_encoded_acc_path in esm_encoded_acc_paths:
        esm_encoded_acc = torch.load(esm_encoded_acc_path)
        acc = esm_encoded_acc["label"]
        #33 is a dict key they use in the esm representations
        esm_representations = esm_encoded_acc["representations"][33]
        #add accesion into labeling for now, so protein names can be retrieved after crossvalidation shuffling.
        binary_labels = np.append(binary_annotations[acc], acc)
        all_binary_labels.append(binary_labels)
        all_esm_representations.append(esm_representations)
    
    all_esm_representations = np.asarray(all_esm_representations, dtype=object)
    all_binary_labels = np.asarray(all_binary_labels, dtype=object)

    return all_esm_representations, all_binary_labels

def crossvalidation_folds(X, y, num_of_folds, save_dirname):
    """
    Inputs: 
        X: Dataset of sequences.
        y: labels
        num_of_fold: Number crossvalidation folds
        save_dirname: save directory
    """
    
    try:
        save_dirname.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print("Save directory was already there. Saving data there.")
    else:
        print("Save directory not found. Made new one. ")
    
    #cross validation setup
    Folds = 1
    kf = KFold(n_splits=num_of_folds)
    kf.get_n_splits(X, y)   
    for train_index, val_index in kf.split(X):
        fold_path = save_dirname / f"Fold{Folds}"
        X_train, X_val = X[train_index], X[val_index]
        y_train, y_val = y[train_index], y[val_index]

        np.savez_compressed(fold_path, saved_X_train=X_train, saved_X_val=X_val,
                            saved_y_train=y_train, saved_y_val=y_val)      
        Folds += 1

def create_crossvalidation_fold(X_train, y_train, X_val, y_val, save_path):
    np.savez_compressed(save_path, saved_X_train=X_train, saved_X_val=X_val,
                                        saved_y_train=y_train, saved_y_val=y_val)

### CLUSTER50ID DATASET ###
test_esm_rep, test_labels = get_esm_representations_and_labels(ESM_EMBEDDINGS / "Test", TEST_EPITOPE_ANNOTATIONS)
cluster50id_embeddings, cluster50id_labels = get_esm_representations_and_labels(ESM_EMBEDDINGS / "ClusterId50", EPITOPE_ANNOTATIONS)
crossvalidation_folds(cluster50id_embeddings, cluster50id_labels, 5, CROSS_VALIDATION_SAVE_PATH)
#save same test data in this crossvalidation
save_data(test_esm_rep, test_labels, CROSS_VALIDATION_SAVE_PATH, "test")
