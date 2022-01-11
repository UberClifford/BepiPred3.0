### IMPORTS ###
import pickle
from pathlib import Path
import sys
import subprocess
import argparse

### STATIC PATHS ###
ROOT_DIR = Path.cwd()
HMM = ROOT_DIR / "HMM"

#antigens we want to have in the test set (5 antigens from earlier studies that did not include antibody in model)
test_antigens = set(["7lj4_B", "4xak_A", "4ypg_D", "7chz_I", "5d8j_A",
                     "3pnw_O", "6y6c_A","5f72_K","6u6u_R","4qci_D",
                     "7jum_A","5th9_A","6hga_B","2xwt_C", "6vtw_A"])


### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Takes pickle formated files, single chain ab-ag complexes, converts to probable light and heavy chain. And joins with ag-ab light and heavy chain separated files.")
parser.add_argument("-i", required=True, action="store", dest="infile", type=Path, help="Input file. Should be pickle formatted file")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target directory. Will write results to this directory.")
parser.add_argument("-a", required=True, action="store", dest="ag_ab_interactions", type=Path, help="Antibody and antigens interactions where light and heavy, are separately annotated.")

args = parser.parse_args()
infile = args.infile
target_dir = args.target_dir
ag_ab_interactions = args.ag_ab_interactions

if not infile.is_file():
    sys.exit(f"The input file was not a valid file.\nInput file given:{infile.name})")
if not target_dir.is_dir():
    target_dir.mkdir(exist_ok=False, parents=True)

# list of tuples = [(heavy_light_acc, heavy_light_chain, antigen_acc, antigen chain)..]
single_chain_ab_interactions = pickle.load(open(infile, "rb"))
ag_ab_interactions = pickle.load(open(ag_ab_interactions, "rb"))

### FUNCTIONS ###
def add_unique_agab_complex_identifiers(ag_ab_interactions, test_antigens):
    """
    Adding some ag-ab complex identifers to keep track of clustering in next step. 
    """
    i=0
    new_ag_ab_interactions = list()

    for ag_ab_interaction in ag_ab_interactions:
        light_chain_id = ag_ab_interaction[0]
        light_chain = ag_ab_interaction[1]
        heavy_chain_id = ag_ab_interaction[2]
        heavy_chain = ag_ab_interaction[3]
        ag_acc = ag_ab_interaction[4]
        ag_chain = ag_ab_interaction[5]

        if ag_acc in test_antigens:
            unique_id = f"Test_ComplexID_{i}"
            print(ag_acc)
        else:
            unique_id = f"ComplexID_{i}"
        new_ag_ab_interactions.append((light_chain_id, light_chain,
                                       heavy_chain_id, heavy_chain,
                                       ag_acc, ag_chain, unique_id))

        i += 1

    return new_ag_ab_interactions

def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

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

def convert_single_chain_ab_to_l_and_h_chains(single_chain_ab_interactions):
    """
    Inputs : single_chain_ab_interactions: list of tuples. [(light_heavy_acc,light_heavy_chain, ab_acc, ab_chain)]
    Outputs: list of tuples [(light_chain_id, light_chain_annotated_string, heavy_chain_id, heavy_chain_annotated_string, antigen_id, antigen_annotated_string)]
    """

    new_ag_ab_interactions = list()
    l_and_h_chain_not_discerned = 0

    #run HMM search to heavy, lambda and kappa chains
    for single_chain_ab_interaction in single_chain_ab_interactions:
       
        #write single chain AB to fasta
        light_heavy_acc = single_chain_ab_interaction[0]
        light_heavy_chain = single_chain_ab_interaction[1]
        fasta_content = [(light_heavy_acc, light_heavy_chain)]
        data_to_fasta_format(fasta_content, "heavy_light_chain_temp.fasta")
        ag_acc = single_chain_ab_interaction[2]
        ag_chain = single_chain_ab_interaction[3]
        
        #do heavy and kappa HMM profile matching 
        subprocess.run(["hmmalign", "--outformat", "A2M", "-o", "heavy_chain_match.fasta", HMM / "heavy.hmm", "heavy_light_chain_temp.fasta"])
        subprocess.run(["hmmalign", "--outformat", "A2M", "-o", "light_chain_match.fasta", HMM / "kappa.hmm", "heavy_light_chain_temp.fasta"])
    
        #read heavy and light chain matches
        with open("heavy_chain_match.fasta", "r") as readfile:
            h_chain_match = read_accs_and_sequences_from_fasta(readfile)
            h_seq_match = h_chain_match[0][1].replace("-", "")
        with open("light_chain_match.fasta", "r") as readfile:
            l_chain_match = read_accs_and_sequences_from_fasta(readfile)
            l_seq_match = l_chain_match[0][1].replace("-", "")
    
    
        #should be of same length
        if len( set( (len(h_seq_match), len(l_seq_match), len(light_heavy_chain)) ) ) != 1:
            print("Something is wrong. Length of light chain match, heavy chain match and original chain don't match.")
    
        chain_length = len(light_heavy_chain)
        light_chain = ""
        heavy_chain = ""

        l_and_h_chain_match = False
    
        for i in range(chain_length):
            #matched both chains
            if h_seq_match[i].isupper() and l_seq_match[i].isupper():
                l_and_h_chain_match = True
                print("Warning: AA in single chain matched both heavy and light chain.")

            #matches heavy or light chain
            elif h_seq_match[i].isupper():
                heavy_chain += light_heavy_chain[i]
            elif l_seq_match[i].isupper():
                light_chain += light_heavy_chain[i]
            
            #doesn't match heavy or light chain
            elif h_seq_match[i].islower() and l_seq_match[i].islower():
                print("Warning: AA not matched light or heavy chain. Probable linker AA, not including this in light or heavy chain.")
                heavy_chain += light_heavy_chain[i]
                light_chain += light_heavy_chain[i]

        if l_and_h_chain_match:
            print("Chain contained AA's, where light and heavy could not be properly discerned. Not including this antibody.")
            l_and_h_chain_not_discerned += 1
        else:
            light_chain_id = light_heavy_acc + "_L"
            heavy_chain_id = light_heavy_acc + "_H"
            new_ag_ab_interactions.append((light_chain_id, light_chain,
                                           heavy_chain_id, heavy_chain,
                                           ag_acc, ag_chain)) 

    return new_ag_ab_interactions, l_and_h_chain_not_discerned

### MAIN ###

#get heavy and light chain HMM's
#get heavy, lambda kappa chains HMM's from lyra github
print("Downloading HMM models from LYRA github")
subprocess.run(["wget", "-P", HMM, "https://raw.githubusercontent.com/paolomarcatili/lyra/master/bcr_models/data/heavy.hmm"])
subprocess.run(["wget", "-P", HMM, "https://raw.githubusercontent.com/paolomarcatili/lyra/master/bcr_models/data/kappa.hmm"])
new_ag_ab_interactions, l_and_h_chain_not_discerned = convert_single_chain_ab_to_l_and_h_chains(single_chain_ab_interactions)

#combine converted ag-ab interactions ag-ab interactions from earlier
combined_ag_ab_interactions = ag_ab_interactions + new_ag_ab_interactions
print(f"Detected single chain AG-AB interactions: {len(single_chain_ab_interactions)}")
print(f"AG-AB interactions from earlier (light and heavy chain annotated were separately): {len(ag_ab_interactions)}")
print(f"Converted single chain AG-AB interaction to separate light and heavy chains: {len(new_ag_ab_interactions)}")
print(f"Single chain AG-AB complexes where AB could not be properly separated into light and heavy chain: {l_and_h_chain_not_discerned}")
print(f"New total AG-AB interactions: {len(combined_ag_ab_interactions)}")

#adding some unique ag-ab complex identifiers
combined_ag_ab_interactions = add_unique_agab_complex_identifiers(combined_ag_ab_interactions, test_antigens)
#print(combined_ag_ab_interactions)
### Save data ###
#save combined antigen-antibody interactions
with open(target_dir / "antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(combined_ag_ab_interactions, outfile)
