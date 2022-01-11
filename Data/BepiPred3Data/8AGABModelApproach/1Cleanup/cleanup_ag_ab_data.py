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


def check_for_epitopes(pdb_seq):
    """
    Checking that there is at least one epitope annotation.
    """
    pos = "foo"
    for res in pdb_seq:
        check = res.isupper()
        if check:
            pos = None
            break
    return pos


def filter_unknowns_and_peptides(ag_ab_interactions):

    filtered_ag_ab_interactions = list()
    filtered_ag_single_chain_ab_interactions = list()
    for ag_ab_interaction in ag_ab_interactions:
        #if light and heavy chain are separately annotated
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
            #check for maximum length
            check7 = check_maximum_length(light_chain_anno_string)
            check8 = check_maximum_length(heavy_chain_anno_string)
            check9 = check_maximum_length(antigen_anno_string)

            #check that there is at least one epitope residue
            check10 = check_for_epitopes(antigen_anno_string)
    
            #if there aren't unknowns in light, heavy or ag chain
            if not any(check != None for check in (check1, check2, check3, check4, check5, check6, check7, check8, check9, check10)):
                filtered_ag_ab_interactions.append((ag_ab_interaction))


        #if light and heavy chain are not separately annoated
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
            #check for maximum length
            check5 = check_maximum_length(heavy_light_chain_anno_string)
            check6 = check_maximum_length(antigen_anno_string)
            #check that there is at least one epitope residue
            check7 = check_for_epitopes(antigen_anno_string)

            if not any(check != None for check in (check1, check2, check3, check4, check5, check6, check7)):
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
print(f"Number ag-ab interactions before filtering out unknowns, peptides and complete negatives {len(filtered_ag_ab_interactions) + len(filtered_ag_single_chain_ab_interactions)}")
#filter redundandt ag-ab complexes
filtered_ag_ab_interactions = filter_redundant_agab_complexes(filtered_ag_ab_interactions)
filtered_ag_single_chain_ab_interactions = filter_redundant_ag_singlechain_ab_complexes(filtered_ag_single_chain_ab_interactions)
print(f"Number ag-ab interactions before filtering completely redundant ag-ab complexes {len(filtered_ag_ab_interactions) + len(filtered_ag_single_chain_ab_interactions)}")

### Save data ###

#save filtered antigen-antibody interactions
with open(target_dir / "filtered_antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(filtered_ag_ab_interactions, outfile)
#save filterd s antigen-singlechainantibody interactions
with open(target_dir / "filtered_antigens_single_chain_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(filtered_ag_single_chain_ab_interactions, outfile)
