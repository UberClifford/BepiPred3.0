## IMPORTS ###
from pathlib import Path
import numpy as np
import sys
import torch
import argparse
import pickle
import numpy as np

### STATICS PATHS ###
ROOT_DIR = Path.cwd()

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Takes pickle formated files, single chain ab-ag complexes, converts to probable light and heavy chain. And joins with ag-ab light and heavy chain separated files.")
parser.add_argument("-i", required=True, action="store", dest="in_dir", type=Path, help="Input directroy containing antibody-antigen interactions. Should be pickle formatted files")
parser.add_argument("-e", required=True, action="store", dest="esm_embeddings_dir", type=Path, help="Directory containing light, heavy and antigen ESM embeddings for all 5 partiions and test set")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target directory. Will write results to this directory.")
parser.add_argument('-simple_approach', action='store_true', dest="simple_approach", help="Simple approach: Average esm embeddings for light and heavy chain. Else: Positional esm embeddings for light and heavy chain")

args = parser.parse_args()
in_dir = args.in_dir
target_dir = args.target_dir
esm_embeddings_dir = args.esm_embeddings_dir
simple_approach = args.simple_approach

### STATIC STUFF AND SOME ERROR HANDLING ###

fold1_ag_ab_pickle_file = in_dir / "Fold1_antigens_antibody_interactions.pickle"
fold2_ag_ab_pickle_file = in_dir / "Fold2_antigens_antibody_interactions.pickle"
fold3_ag_ab_pickle_file = in_dir / "Fold3_antigens_antibody_interactions.pickle"
fold4_ag_ab_pickle_file = in_dir / "Fold4_antigens_antibody_interactions.pickle"
fold5_ag_ab_pickle_file = in_dir / "Fold5_antigens_antibody_interactions.pickle"
test_ag_ab_pickle_file = in_dir / "test.pickle"
ag_ab_pickle_files = (fold1_ag_ab_pickle_file, fold2_ag_ab_pickle_file, fold3_ag_ab_pickle_file, fold4_ag_ab_pickle_file, fold5_ag_ab_pickle_file, test_ag_ab_pickle_file)

#Need to complete previous step and have have ESM encoded data to run script. 
if not all(ag_ab_pickle_file.is_file() for ag_ab_pickle_file in ag_ab_pickle_files):
    sys.exit("AG-AB interaction pickle data missing. Complete previous step of partition ag-ab interaction data into pickle files.")

fold1_esm_embeddings_dir = esm_embeddings_dir / "Fold1_antigens_antibody_interactions"
fold2_esm_embeddings_dir = esm_embeddings_dir / "Fold2_antigens_antibody_interactions"
fold3_esm_embeddings_dir = esm_embeddings_dir / "Fold3_antigens_antibody_interactions"
fold4_esm_embeddings_dir = esm_embeddings_dir / "Fold4_antigens_antibody_interactions"
fold5_esm_embeddings_dir = esm_embeddings_dir / "Fold5_antigens_antibody_interactions"
test_esm_embeddings_dir = esm_embeddings_dir / "test"
esm_directories = (fold1_esm_embeddings_dir, fold2_esm_embeddings_dir, fold3_esm_embeddings_dir, fold4_esm_embeddings_dir, fold5_esm_embeddings_dir, test_esm_embeddings_dir)

#Need to complete previous step and have have ESM encoded data to run script. 
if not all(esm_dir.is_dir() for esm_dir in esm_directories):
    sys.exit("ESM data missing. Complete previous step of ESM encoding light, heavy and antigen chains of 5 partitions and test, before running this script.")

if not target_dir.is_dir():
    target_dir.mkdir(exist_ok=False, parents=True)

num_of_partitions = len(esm_directories)

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

def save_data(X, save_dirname, filename = "foo"):
    """
    In this case, 
    X: (test_ids, test_labels,  test_l_chain_esm, test_h_chain_esm, test_antigen_esm)
    """

    try:
        save_dirname.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print("Save directory was already there. Saving data there.")
    else:
        print("Save directory not found. Made new one. ")
    
    save_path = save_dirname / filename
    np.savez_compressed(save_path, saved_test_ids = X[0], saved_test_labels = X[1],
                        saved_test_lchain=X[2], saved_test_hchain=X[3], saved_test_antigen=X[4])

def create_crossvalidation_fold(train, val, save_path):
    """
    Inputs: train: (train_ids, train_labels, train_l_chain_esm, train_h_chain_esm, train_antigen_esm )
            val: (val_ids, val_labels, val_l_chain_esm, val_h_chain_esm, val_antigen_esm)
    """

    np.savez_compressed(save_path, saved_train_ids=train[0], saved_train_labels=train[1],
                        saved_train_lchain= train[2],saved_train_hchain= train[3],
                        saved_train_antigen= train[4], saved_val_ids=val[0], saved_val_labels=val[1],
                        saved_val_lchain= val[2],saved_val_hchain= val[3],
                        saved_val_antigen= val[4]) 


def add_seq_len(X, seq_len, avg_embedding=False):

    if avg_embedding:
        seq_len = torch.tensor([seq_len])
        new_X = torch.cat((X, seq_len))
        
    else:
        seq_len_v = torch.ones(seq_len)*seq_len
        seq_len_v = seq_len_v.unsqueeze(dim=1)
        new_X = torch.cat((X, seq_len_v), axis=1)

    return new_X


def get_light_heavy_antigen_esm_and_labels(ag_ab_pickle_file, esm_data_path, simple_approach_data=False):

    ag_ab_interactions = pickle.load(open(ag_ab_pickle_file, "rb"))
    light_chain_esm_embedding_path = esm_data_path / "ESM_light_chains"
    heavy_chain_esm_embeddings_path = esm_data_path / "ESM_heavy_chains"
    antigen_esm_embeddings_path = esm_data_path / "ESM_antigens"
    all_ids = list()
    all_binary_labels =list()
    all_light_chain_esm_representations = list()
    all_heavy_chain_esm_representations = list()
    all_antigen_esm_representations = list()

    for ag_ab_interaction in ag_ab_interactions:
        light_chain_id = ag_ab_interaction[0]
        light_chain = ag_ab_interaction[1]
        heavy_chain_id = ag_ab_interaction[2]
        heavy_chain = ag_ab_interaction[3]
        ag_acc = ag_ab_interaction[4]
        ag_chain = ag_ab_interaction[5]
        unique_complex_id = ag_ab_interaction[6]

        #get sequence lengths
        lc_len = len(light_chain)
        hc_len = len(heavy_chain) 
        ag_len = len(ag_chain)

        #get binary labeling for epitopes
        binary_label = get_epitope_labels_as_binary_array_from_sequence(ag_chain)
        
        #get light, heavy and antigen chain ESM representations
        light_chain_esm_representation = torch.load(light_chain_esm_embedding_path / f"{unique_complex_id}.pt")
        heavy_chain_esm_representation = torch.load(heavy_chain_esm_embeddings_path / f"{unique_complex_id}.pt")
        antigen_esm_representation = torch.load(antigen_esm_embeddings_path / f"{unique_complex_id}.pt")

        if simple_approach_data:
            light_chain_esm_representation = light_chain_esm_representation["mean_representations"][33]
            heavy_chain_esm_representation = heavy_chain_esm_representation["mean_representations"][33]
            #add sequence length to embedding
            light_chain_esm_representation = add_seq_len(light_chain_esm_representation, lc_len, avg_embedding=simple_approach_data)
            heavy_chain_esm_representation = add_seq_len(heavy_chain_esm_representation, hc_len, avg_embedding=simple_approach_data)

        else:
            light_chain_esm_representation = light_chain_esm_representation["representations"][33]
            heavy_chain_esm_representation = heavy_chain_esm_representation["representations"][33]
            #add sequence length to embedding
            light_chain_esm_representation = add_seq_len(light_chain_esm_representation, lc_len)
            heavy_chain_esm_representation = add_seq_len(heavy_chain_esm_representation, hc_len)

        antigen_esm_representation = antigen_esm_representation["representations"][33]
        #add sequence length to embedding
        antigen_esm_representation = add_seq_len(antigen_esm_representation, ag_len)

        #appened info, binary labels and esm representations
        all_ids.append([light_chain_id, heavy_chain_id, ag_acc, unique_complex_id])
        all_binary_labels.append(binary_label)
        all_light_chain_esm_representations.append(light_chain_esm_representation)
        all_heavy_chain_esm_representations.append(heavy_chain_esm_representation)
        all_antigen_esm_representations.append(antigen_esm_representation)


    all_ids = np.asarray(all_ids, dtype=object)
    all_binary_labels = np.asarray(all_binary_labels, dtype=object)
    all_light_chain_esm_representations = np.asarray(all_light_chain_esm_representations, dtype=object)
    all_heavy_chain_esm_representations = np.asarray(all_heavy_chain_esm_representations, dtype=object)
    all_antigen_esm_representations = np.asarray(all_antigen_esm_representations, dtype=object)

    return all_ids, all_binary_labels,  all_light_chain_esm_representations, all_heavy_chain_esm_representations, all_antigen_esm_representations

### MAIN ###

#Partition1
print("Getting (light, heavy, antigen) ESM embeddings and epitope labels for Partition1/Fold1")
partition1_ids, partition1_labels, partition1_l_chain_esm, partition1_h_chain_esm, partition1_antigen_esm = get_light_heavy_antigen_esm_and_labels(fold1_ag_ab_pickle_file, fold1_esm_embeddings_dir, simple_approach_data=simple_approach)
print("Done")
#Partition2
print("Getting (light, heavy, antigen) ESM embeddings and epitope labels for Partition2/Fold2")
partition2_ids, partition2_labels, partition2_l_chain_esm, partition2_h_chain_esm, partition2_antigen_esm = get_light_heavy_antigen_esm_and_labels(fold2_ag_ab_pickle_file, fold2_esm_embeddings_dir, simple_approach_data=simple_approach)
print("Done")
##Partition3
print("Getting (light, heavy, antigen) ESM embeddings and epitope labels for Partition3/Fold3")
partition3_ids, partition3_labels, partition3_l_chain_esm, partition3_h_chain_esm, partition3_antigen_esm = get_light_heavy_antigen_esm_and_labels(fold3_ag_ab_pickle_file, fold3_esm_embeddings_dir, simple_approach_data=simple_approach)
print("Done")
##Partition4
print("Getting (light, heavy, antigen) ESM embeddings and epitope labels for Partition4/Fold4")
partition4_ids, partition4_labels, partition4_l_chain_esm, partition4_h_chain_esm, partition4_antigen_esm = get_light_heavy_antigen_esm_and_labels(fold4_ag_ab_pickle_file, fold4_esm_embeddings_dir, simple_approach_data=simple_approach)
print("Done")
##Partition5
print("Getting (light, heavy, antigen) ESM embeddings and epitope labels for Partition5/Fold5")
partition5_ids, partition5_labels, partition5_l_chain_esm, partition5_h_chain_esm, partition5_antigen_esm = get_light_heavy_antigen_esm_and_labels(fold5_ag_ab_pickle_file, fold5_esm_embeddings_dir, simple_approach_data=simple_approach)
print("Done")

###Test
print("Getting (light, heavy, antigen) ESM embeddings and epitope labels for external test set")
test_ids, test_labels,  test_l_chain_esm, test_h_chain_esm, test_antigen_esm = get_light_heavy_antigen_esm_and_labels(test_ag_ab_pickle_file, test_esm_embeddings_dir, simple_approach_data=simple_approach)
print("Done")
print("Saving external test set")
test = (test_ids, test_labels,  test_l_chain_esm, test_h_chain_esm, test_antigen_esm)
save_data(test, target_dir, filename = "test")
print("Done")

#Fold1 
print("Creating crossvalidation fold 1")
train_ids = np.concatenate( (partition2_ids, partition3_ids, partition4_ids, partition5_ids) )
train_labels = np.concatenate( (partition2_labels, partition3_labels, partition4_labels, partition5_labels) )
train_l_chain_esm = np.concatenate( (partition2_l_chain_esm, partition3_l_chain_esm, partition4_l_chain_esm, partition5_l_chain_esm) )
train_h_chain_esm = np.concatenate( (partition2_h_chain_esm, partition3_h_chain_esm, partition4_h_chain_esm, partition5_h_chain_esm) )
train_antigen_esm = np.concatenate( (partition2_antigen_esm, partition3_antigen_esm, partition4_antigen_esm, partition5_antigen_esm) )

train = (train_ids, train_labels, train_l_chain_esm, train_h_chain_esm, train_antigen_esm) 
val = (partition1_ids, partition1_labels,  partition1_l_chain_esm, partition1_h_chain_esm, partition1_antigen_esm)
create_crossvalidation_fold(train, val, target_dir / "Fold1")
print("Done")

#Fold2
print("Creating crossvalidation fold 2")
train_ids = np.concatenate( (partition1_ids, partition3_ids, partition4_ids, partition5_ids) )
train_labels = np.concatenate( (partition1_labels, partition3_labels, partition4_labels, partition5_labels) )
train_l_chain_esm = np.concatenate( (partition1_l_chain_esm, partition3_l_chain_esm, partition4_l_chain_esm, partition5_l_chain_esm) )
train_h_chain_esm = np.concatenate( (partition1_h_chain_esm, partition3_h_chain_esm, partition4_h_chain_esm, partition5_h_chain_esm) )
train_antigen_esm = np.concatenate( (partition1_antigen_esm, partition3_antigen_esm, partition4_antigen_esm, partition5_antigen_esm) )

train = (train_ids, train_labels, train_l_chain_esm, train_h_chain_esm, train_antigen_esm) 
val = (partition2_ids, partition2_labels,  partition2_l_chain_esm, partition2_h_chain_esm, partition2_antigen_esm)
create_crossvalidation_fold(train, val, target_dir / "Fold2")
print("Done")

#Fold3
print("Creating crossvalidation fold 3")
train_ids = np.concatenate( (partition2_ids, partition1_ids, partition4_ids, partition5_ids) )
train_labels = np.concatenate( (partition2_labels, partition1_labels, partition4_labels, partition5_labels) )
train_l_chain_esm = np.concatenate( (partition2_l_chain_esm, partition1_l_chain_esm, partition4_l_chain_esm, partition5_l_chain_esm) )
train_h_chain_esm = np.concatenate( (partition2_h_chain_esm, partition1_h_chain_esm, partition4_h_chain_esm, partition5_h_chain_esm) )
train_antigen_esm = np.concatenate( (partition2_antigen_esm, partition1_antigen_esm, partition4_antigen_esm, partition5_antigen_esm) )

train = (train_ids, train_labels, train_l_chain_esm, train_h_chain_esm, train_antigen_esm) 
val = (partition3_ids, partition3_labels,  partition3_l_chain_esm, partition3_h_chain_esm, partition3_antigen_esm)
create_crossvalidation_fold(train, val, target_dir / "Fold3")
print("Done")

#Fold4
print("Creating crossvalidation fold 4")
train_ids = np.concatenate( (partition2_ids, partition3_ids, partition1_ids, partition5_ids) )
train_labels = np.concatenate( (partition2_labels, partition3_labels, partition1_labels, partition5_labels) )
train_l_chain_esm = np.concatenate( (partition2_l_chain_esm, partition3_l_chain_esm, partition1_l_chain_esm, partition5_l_chain_esm) )
train_h_chain_esm = np.concatenate( (partition2_h_chain_esm, partition3_h_chain_esm, partition1_h_chain_esm, partition5_h_chain_esm) )
train_antigen_esm = np.concatenate( (partition2_antigen_esm, partition3_antigen_esm, partition1_antigen_esm, partition5_antigen_esm) )

train = (train_ids, train_labels, train_l_chain_esm, train_h_chain_esm, train_antigen_esm) 
val = (partition4_ids, partition4_labels,  partition4_l_chain_esm, partition4_h_chain_esm, partition4_antigen_esm)
create_crossvalidation_fold(train, val, target_dir / "Fold4")
print("Done")

#Fold5
print("Creating crossvalidation fold 5")
train_ids = np.concatenate( (partition2_ids, partition3_ids, partition4_ids, partition1_ids) )
train_labels = np.concatenate( (partition2_labels, partition3_labels, partition4_labels, partition1_labels) )
train_l_chain_esm = np.concatenate( (partition2_l_chain_esm, partition3_l_chain_esm, partition4_l_chain_esm, partition1_l_chain_esm) )
train_h_chain_esm = np.concatenate( (partition2_h_chain_esm, partition3_h_chain_esm, partition4_h_chain_esm, partition1_h_chain_esm) )
train_antigen_esm = np.concatenate( (partition2_antigen_esm, partition3_antigen_esm, partition4_antigen_esm, partition1_antigen_esm) )

train = (train_ids, train_labels, train_l_chain_esm, train_h_chain_esm, train_antigen_esm) 
val = (partition5_ids, partition5_labels,  partition5_l_chain_esm, partition5_h_chain_esm, partition5_antigen_esm)
create_crossvalidation_fold(train, val, target_dir / "Fold5")
print("Done")
