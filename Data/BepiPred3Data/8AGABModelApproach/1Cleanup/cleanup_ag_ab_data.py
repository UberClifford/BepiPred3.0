import pickle
from pathlib import Path
import sys
import argparse


### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Takes pickle formated files, removes complexes with DNA or unknown characters, completely redundant (case-senstive) ag-ag complexes")
parser.add_argument("-i", required=True, action="store", dest="infile", type=Path, help="Input file. Should be pickle formatted file")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target directory. Will write results to this directory.")
parser.add_argument("-pl", required=False, action="store", dest="peptide_len", type=int, default = 40, help="Peptide length. Entries with sequence length lower than this, will be regarded as a peptide and separated from main results.")

args = parser.parse_args()
infile = args.infile
target_dir = args.target_dir
peptide_len = args.peptide_len

if not infile.is_file():
    sys.exit(f"The input file was not a valid file.\nInput file given:{infile.name})")
if not target_dir.is_dir():
    target_dir.mkdir(exist_ok=False, parents=True)

ag_ab_interactions = pickle.load(open(infile, "rb"))
### Data input ###
#pickle,
#list of tuples with dim: [(light_chain_id, light_chain_annotated_string, heavy_chain_id, heavy_chain_annotated_string, antigen_id, antigen_annotated_string)]
#                     or: [(heavy_light_chain_id, heavy_light_chain_annotated_seq, antigen_id, antigen_annotated_seq)]


### FUNCTIONS ###

def detect_sequence_with_unknowns(pdb_seq):
    lowercase_pdb_seq = pdb_seq.lower()
    seq_with_unknowns = None
    
    if "x" in lowercase_pdb_seq:
        seq_with_unknowns = pdb_seq
    else:
        pass

    return seq_with_unknowns

def detect_peptide(pdb_seq, allowed_protein_len = 40):
    seq_len = len(pdb_seq)
    peptide = None
    if seq_len < allowed_protein_len:
        peptide = pdb_seq
    else:
        pass
    return peptide

def filter_unknowns_and_peptides(ag_ab_interactions):

    filtered_ag_ab_interactions = list()
    filtered_ag_single_chain_ab_interactions = list()
    for ag_ab_interaction in ag_ab_interactions:
        if len(ag_ab_interaction) == 6:
            light_chain_id = ag_ab_interaction[0]
            light_chain_anno_string = ag_ab_interaction[1]
            heavy_chain_id = ag_ab_interaction[2]
            heavy_chain_anno_string = ag_ab_interaction[3]
            antigen_id = ag_ab_interaction[4]
            antigen_anno_string = ag_ab_interaction[5]
    
            #check for unknown seuqences
            check1 = detect_sequence_with_unknowns(light_chain_anno_string)
            check2 = detect_sequence_with_unknowns(heavy_chain_anno_string)
            check3 = detect_sequence_with_unknowns(antigen_anno_string)
            #check for peptides
            check4 = detect_peptide(light_chain_anno_string)
            check5 = detect_peptide(heavy_chain_anno_string)
            check6 = detect_peptide(antigen_anno_string)
    
            #if there aren't unknowns in light, heavy or ag chain
            if not any(check != None for check in (check1, check2, check3, check4, check5, check6)):
                filtered_ag_ab_interactions.append((ag_ab_interaction))
    
    
        elif len(ag_ab_interaction) == 4:
            heavy_light_chain_id = ag_ab_interaction[0]
            heavy_light_chain_anno_string = ag_ab_interaction[1]
            antigen_id = ag_ab_interaction[2]
            antigen_anno_string = ag_ab_interaction[3]
    
            #check for unknown sequences
            check1 = detect_sequence_with_unknowns(heavy_light_chain_anno_string)
            check2 = detect_sequence_with_unknowns(antigen_anno_string)
            #check for peptides
            check3 = detect_peptide(heavy_light_chain_anno_string)
            check4 = detect_peptide(antigen_anno_string)
    
            if not any(check != None for check in (check1, check2, check3, check4)):
                filtered_ag_single_chain_ab_interactions.append((ag_ab_interaction))
    
        else:
            print("Something is wrong!")

    return filtered_ag_ab_interactions, filtered_ag_single_chain_ab_interactions

def filter_redundant_agab_complexes(ag_ab_interactions, case_sensitive=False):

    redundant_check_list = list()
    filtered_ag_ab_interactions = list()

    for ag_ab_interaction in ag_ab_interactions:
        light_chain_id = ag_ab_interaction[0]
        light_chain_anno_string = ag_ab_interaction[1]
        heavy_chain_id = ag_ab_interaction[2]
        heavy_chain_anno_string = ag_ab_interaction[3]
        antigen_id = ag_ab_interaction[4]
        antigen_anno_string = ag_ab_interaction[5]
        if case_sensitive:
            interaction_concatenated = light_chain_anno_string + heavy_chain_anno_string + antigen_anno_string
        else:
            interaction_concatenated = light_chain_anno_string.lower() + heavy_chain_anno_string.lower() + antigen_anno_string.lower()

        if interaction_concatenated not in redundant_check_list:
            filtered_ag_ab_interactions.append(ag_ab_interaction)
            redundant_check_list.append(interaction_concatenated)


    return filtered_ag_ab_interactions

def filter_redundant_ag_singlechain_ab_complexes(ag_ab_interactions, case_sensitive=False):

    redundant_check_list = list()
    filtered_ag_ab_interactions = list()

    for ag_ab_interaction in ag_ab_interactions:
        heavy_light_chain_id = ag_ab_interaction[0]
        heavy_light_chain_anno_string = ag_ab_interaction[1]
        antigen_id = ag_ab_interaction[2]
        antigen_anno_string = ag_ab_interaction[3]

        if case_sensitive:
            interaction_concatenated = heavy_light_chain_anno_string + antigen_anno_string
        else:
            interaction_concatenated = heavy_light_chain_anno_string.lower() + antigen_anno_string.lower()

        if interaction_concatenated not in redundant_check_list:
            filtered_ag_ab_interactions.append(ag_ab_interaction)
            redundant_check_list.append(interaction_concatenated)


    return filtered_ag_ab_interactions

### MAIN ###

print(f"Number ag-ab interactions before filtering anything {len(ag_ab_interactions)}")
#filter out unknowns and peptides
filtered_ag_ab_interactions, filtered_ag_single_chain_ab_interactions = filter_unknowns_and_peptides(ag_ab_interactions) 
print(f"Number ag-ab interactions before filtering out unknowns and peptides {len(filtered_ag_ab_interactions) + len(filtered_ag_single_chain_ab_interactions)}")
#filter redundandt ag-ab complexes
filtered_ag_ab_interactions = filter_redundant_agab_complexes(filtered_ag_ab_interactions)
filtered_ag_single_chain_ab_interactions = filter_redundant_ag_singlechain_ab_complexes(filtered_ag_single_chain_ab_interactions)
print(f"Number ag-ab interactions before filtering out unknowns and peptides {len(filtered_ag_ab_interactions) + len(filtered_ag_single_chain_ab_interactions)}")

### Save data ###

#save filtered antigen-antibody interactions
with open(target_dir / "filtered_antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(filtered_ag_ab_interactions, outfile)
#save filterd s antigen-singlechainantibody interactions
with open(target_dir / "filtered_antigens_single_chain_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(filtered_ag_single_chain_ab_interactions, outfile)