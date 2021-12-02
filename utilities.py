### IMPORTS ###
import numpy as np
import torch
import torch.nn as nn
import math
from torch.nn import CrossEntropyLoss, Softmax, BatchNorm1d, BCEWithLogitsLoss, Sigmoid
import pickle
from scipy import stats
from torch.utils.data import DataLoader, Dataset
from torch.nn import functional as F
from sklearn import metrics
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from torch.nn.utils.rnn import pad_sequence, pack_padded_sequence, pad_packed_sequence
import numpy as np
from torch.optim import Adam
from pathlib import Path
import random
import warnings



### SETTING SEEDS ###
random.seed(10)
np.random.seed(10)
torch.manual_seed(10)

### STATIC PATHS ###
ROOT_DIR = Path.cwd()
DATA_DIR = ROOT_DIR / "Data"  
BEPIPRED3_NON_HOM_REDUCED = DATA_DIR / "BepiPred3Data" / "5DatasetPreparation" / "7CrossValidationNonHomologyReduced"
BEPIPRED3_HOM_REDUCED = DATA_DIR / "BepiPred3Data" / "5DatasetPreparation" / "8CrossValidationHomologyReduced"
BEPIPRED3_CLUS50ID = DATA_DIR / "BepiPred3Data" / "6Clusterat50ID" / "5CrossValidationClusterAt50Id"

RESULTS_DIR = ROOT_DIR / "Results"
FIGURE_DIR = RESULTS_DIR / "NeuralNetworks" / "Figures"

### LOADING DATA AND BATCH LOADER FUNCTIONS ###
def load_data(data_path, test_data=False):
    """
    """
    
    loaded_data = np.load(data_path, allow_pickle=True)
    
    if test_data:
        X_test = loaded_data["saved_X"]
        y_test = loaded_data["saved_y"]
        return X_test, y_test
    else:
        X_train = loaded_data["saved_X_train"]
        y_train = loaded_data["saved_y_train"]
        X_val = loaded_data["saved_X_val"]
        y_val = loaded_data["saved_y_val"]
        return X_train, y_train, X_val, y_val

def separate_acc_names_from_labels(y_data):
    """
    """
    acc_names = [label[-1] for label in y_data]
    y_with_no_acc_names = np.asarray([label[:-1] for label in y_data], dtype="object")
    
    return y_with_no_acc_names, acc_names

def compute_no_pred_power_performance(labels):
    """
    No prediction power --> Just predicting epitope residue frequency, which is approximately
    0.1.
    """
    no_pred_power_dict = {0:0.9, 1:0.1}
    
    labels = np.concatenate(labels, axis=0)
    labels = np.asarray(labels, dtype=np.uint8)
    #summed negative log likelihood
    no_pred_power_peformance = sum([-math.log(no_pred_power_dict[label]) for label in labels ]) / len(labels)
    
    return no_pred_power_peformance
    
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


def blosum50_or_one_hot_pad_sequence(sequence, window_size = 9):
    """
    Simple padding scheme used for blosum62 or sparse encodings.
    """
    
    #add padding
    seq_windows = list()
    padding = "".join(["-" for i in range(int(window_size / 2))])
    padded_sequence = padding + sequence + padding
    seq_len = len(padded_sequence)
    for i in range(seq_len - window_size + 1):
        window = padded_sequence[i:i + window_size]
        seq_windows.append(window)
    
    return seq_windows

def blosum50_or_one_hot_encode_data(y_accs_names, sequence_and_accs, encode_dict):
    
    i = 0
    accs = [sequence_and_acc[0] for sequence_and_acc in sequence_and_accs]
    sequences = [sequence_and_acc[1] for sequence_and_acc in sequence_and_accs]
    seq_encodings = list()
    y_binary_labels = list()
    
    for y_acc_name in y_accs_names:
        try:
            idx = accs.index(y_acc_name)
        except ValueError:
            idx = None
        if idx is not None:
            acc_name = accs[idx]
            sequence = sequences[idx]        
            seq_windows = blosum50_or_one_hot_pad_sequence(sequence.upper())
            seq_encoding = list()
            for seq_window in seq_windows:
                position_encoding = [torch.tensor(encode_dict[res]).float() for res in seq_window]
                position_encoding = torch.cat(position_encoding)
                seq_encoding.append(position_encoding)
                
            seq_encoding = torch.stack(seq_encoding)
            seq_encodings.append(seq_encoding)
            y_binary_label = get_epitope_labels_as_binary_array_from_sequence(sequence)
            y_binary_labels.append(y_binary_label)
            
    return seq_encodings, y_binary_labels 

def add_seq_len_feature(X_data):
    #get sequence lengths
    new_X_data = list()
    for X in X_data:
        seq_len = X.size()[0]
        seq_len_v = torch.ones(seq_len)*seq_len
        seq_len_v = seq_len_v.unsqueeze(dim=1)
        new_X = torch.cat((X, seq_len_v), axis=1)
        new_X_data.append(new_X)
        
    return new_X_data

def get_seq_indexes(whole_file):
    """
    Inputs: Whole file in binary.
    Outputs: sequence start and end indices. List of tuples. 
    """

    #the first entry header start and end indices
    header_starts = [whole_file.find(b">")]
    header_ends = [whole_file.find(b"\n", header_starts[-1]) - 1] # the next newline character after header start, is the end of the header.
    #remaining entry header start and end indices
    while True:
        header_start = whole_file.find(b">", header_ends[-1])
        #if there are no more entries, break
        if header_start == -1:
            break
        header_starts.append(header_start)
        header_ends.append(whole_file.find(b"\n", header_starts[-1]) - 1)

    ##infer entry sequence start and end indices, from header starts and ends indices
    number_of_entries = len(header_starts)
    seq_starts = list()
    seq_ends = list()

    #infer all entries
    for i in range(number_of_entries):
        seq_starts.append(header_ends[i] + 2) #skipping newline characters

        #if at last entry
        if i == number_of_entries - 1:
            seq_ends.append(len(whole_file))
        #other entries
        else:
            seq_ends.append(header_starts[i + 1] - 1)
        

    return list( zip(seq_starts, seq_ends) )

#def blosum50_or_one_hot_encode_data(y_accs_names, sequence_and_accs, encode_dict):
#    
#    y_accs_names = set(y_accs_names)
#    i = 0
#    seq_encodings = list()
#    y_binary_labels = list()
#    
#    for sequence_and_acc in sequence_and_accs:
#        acc_name = sequence_and_acc[0]
#        sequence = sequence_and_acc[1]
#        if acc_name in y_accs_names:
#            i += 1
#            seq_windows = blosum50_or_one_hot_pad_sequence(sequence.upper())
#            seq_encoding = list()
#            for seq_window in seq_windows:
#                position_encoding = [torch.tensor(encode_dict[res]).float() for res in seq_window]
#                position_encoding = torch.cat(position_encoding)
#                seq_encoding.append(position_encoding)
#                
#            seq_encoding = torch.stack(seq_encoding)
#            seq_encodings.append(seq_encoding)
#            y_binary_label = get_epitope_labels_as_binary_array_from_sequence(sequence)
#            y_binary_labels.append(y_binary_label)
#            
#    return seq_encodings, y_binary_labels    
            
def get_class_weights(y_train, y_val):
    """
    Inputs: combined_y: Array of array with different length. Each array contains non-epitope (0) and epitope (1) annotations.
    Outputs: class_weights: List: [non_epitope_weight, epitope_weight] for CrossEntropyLoss
             pos_weight: float:  For BCEWithLogitLoss
    """
    
    #counting epitopes and non-epitopes
    non_epitopes = 0
    epitopes = 0
    for arr in y_train:
        count_array = np.asarray(arr, dtype="uint8")
        count_array = np.bincount(count_array)
        non_epitopes += count_array[0]
        epitopes += count_array[1]
    
    for arr in y_val:
        count_array = np.asarray(arr, dtype="uint8")
        count_array = np.bincount(count_array)
        non_epitopes += count_array[0]
        epitopes += count_array[1]
        
    #creating class weights
    class_weights = torch.tensor([non_epitopes, epitopes])
    class_weights = class_weights / class_weights.sum()
    class_weights = 1 / class_weights
    class_weights = class_weights / class_weights.sum()
    #creating pos_weight
    pos_weight = torch.tensor([non_epitopes / epitopes])
    
    return class_weights, pos_weight
        
### EVALUATION FUNCTIONS ###

def get_fpr_tpr_auc_and_opt_threshold(y_true, y_pos_prob):
    fpr, tpr, thresh = metrics.roc_curve(y_true, y_pos_prob)
    auc = metrics.auc(fpr, tpr)
    optimal_thresh = thresh[np.argmax(tpr - fpr)]    
    return fpr, tpr, auc, optimal_thresh

def metrics_with_threshold(y_true, y_pos_prob, optimal_thresh, ensemble = None):
    
    if ensemble == None:
        y_preds = [1 if res >= optimal_thresh else 0 for res in y_pos_prob]
    else:
        y_preds = ensemble
        
    acc = metrics.accuracy_score(y_true, y_preds)
    mcc = metrics.matthews_corrcoef(y_true, y_preds)
    
    recall = metrics.recall_score(y_true, y_preds)
    precision = metrics.precision_score(y_true, y_preds)
    f1_score = metrics.f1_score(y_true, y_preds)
    
    return acc, mcc, recall, precision, f1_score, y_preds

def save_pytorch_model(model, model_save_path):
    model_save_dir = model_save_path.parent
    
    try:
        model_save_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print("Save directory was already there. Saving it there.")
    else:
        print("Save directory not found. Made new one. ")
    torch.save(model.state_dict(), model_save_path)

    
def create_figure_save_path(figure_save_path):
    try:
        figure_save_path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print("Save directory was already there. Saving it there.")
    else:
        print("Save directory not found. Made new one. ")
    

def get_optimal_mcc_score(y_true, y_pos_prob):
    
    #trying 50 different mcc thresholds
    mcc_thresholds = np.linspace(0, 0.95, 50)
    best_mcc_threshold = 0.5
    opt_mcc = 0
    acc = metrics.accuracy_score(y_true, [1 if res >= 0.5 else 0 for res in y_pos_prob])
    
    for mcc_threshold in mcc_thresholds:
        y_preds = [1 if res >= mcc_threshold else 0 for res in y_pos_prob]
        mcc = metrics.matthews_corrcoef(y_true, y_preds)
        
        if mcc > opt_mcc:
            opt_mcc = mcc
            best_mcc_threshold = mcc_threshold
            acc = metrics.accuracy_score(y_true, y_preds)
            
    
    
    return best_mcc_threshold, opt_mcc, acc
    
def get_performance_metrics(model_output, labels):
    """
    Inputs: model_output
    Outputs: Accuracy, AUC and MCC using unoptimized threshold. 
    """
    y_true, y_pos_prob_no_padding = get_labels_preds_and_posprob_without_padding(model_output, labels)
   #accuracy, mcc and AUC
    fpr, tpr, thresh = metrics.roc_curve(y_true, y_pos_prob_no_padding)
    auc = metrics.auc(fpr, tpr)
    best_mcc_threshold, opt_mcc, acc = get_optimal_mcc_score(y_true, y_pos_prob_no_padding)
    
#   sigmoid = Sigmoid()
#
#    y_true = labels.cpu().detach().numpy()
#    y_pos_prob = sigmoid(model_output).cpu().detach().numpy()
#    y_preds = [1 if res >= 0.5 else 0 for res in y_pos_prob]
#
#    acc = metrics.accuracy_score(y_true, y_preds)
#    mcc = metrics.matthews_corrcoef(y_true, y_preds)
#    fpr, tpr, thresh = metrics.roc_curve(y_true, y_pos_prob)
#    auc = metrics.auc(fpr, tpr)
    
    return acc, auc, opt_mcc, best_mcc_threshold
    
    
def compute_auc10(fpr, tpr):
    
    #slice fpr for all values where fpr < 0.1
    fpr_auc10 = [x for x in fpr if x <= 0.1]
    v_len = len(fpr_auc10)
    #get interval between fprs (+ interval from 0 to first fpr) 
    fpr_auc10_intervals = [fpr_auc10[0]] + [fpr_auc10[i+1] - fpr_auc10[i] for i in range(v_len - 1)]
    #get same size from tpr
    tpr_auc10 = tpr[:v_len]
    #sum of products of these is the AUC10 score
    auc10 = sum([fpr_auc10_intervals[i]*tpr_auc10[i] for i in range(v_len)])
    #normalize
    auc10 = auc10 / 0.1
    
    return auc10 



def get_labels_preds_and_posprob_without_padding(model_output, labels):
    
    y_true = list()
    y_pos_prob_no_padding = list()
    y_preds_no_padding = list()
    
    labels = labels.cpu().detach().numpy()
    num_of_labels = len(labels)
    softmax_function = Softmax(dim=1)
    class_probs = softmax_function(model_output)
    all_pos_probs = class_probs[:, 1].cpu().detach().numpy()
#    _, preds = torch.max(class_probs, 1)
    #epitope probs (need to ignoore padding from these)
    
    #get true labels, predictions and positive probs. without padding = 2
    for i in range(num_of_labels):
        #if padding
        if labels[i] == 2:
            pass
        #if non-epitope or epitope
        else:
            y_true.append(labels[i])
            y_pos_prob_no_padding.append(all_pos_probs[i])
#            y_preds_no_padding.append(preds[i])
            
    return y_true, y_pos_prob_no_padding #, y_preds_no_padding
    
  
    