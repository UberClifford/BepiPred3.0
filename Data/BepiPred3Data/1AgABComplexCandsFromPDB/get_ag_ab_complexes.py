### REQUIREMENTS ###

#requires hmmer
#conda install -c bioconda hmmer

#created in Python 3.8.8. But other Python3 versions should work. 

### IMPORTS ###
import subprocess
from collections import Counter
from pathlib import Path
import argparse
import shutil
import pickle

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Check PDB database for possible antigen-antibody complexes. Matches are found by comparing PDB sequence data to"\
        " kappa, lambda and heavy chain hidden markov models (HMMs). Both the PDB database sequence data and HMM's are automatically downloaded."\
        " The final output is a file containing possible antibody-antigen complex pdb accesions. These are separated by newline character.")

parser.add_argument("-e", required=False, action="store", dest="HMM_expected_value", default=1e-18, type=float,
        help="Set the HMM expected value. This sets how significant a match has to be, before it is considered significant")

parser.add_argument("-c", required=False, action="store_true", dest="clean_up", help="If specified, the temporary results, downloaded PDB sequence data and"\
        " hidden markov models will be deleted on completion.")

args = parser.parse_args()
HMM_expected_value = args.HMM_expected_value
clean_up = args.clean_up

#Example 1: Run with default options
#python ag_ab_complexes.py

#Example 2: Run with different HMM expected value and clean up option
#python ag_ab_complexes.py -e 1e-30 -c

### STATIC PATHS ###
ROOT_DIR = Path.cwd()
PDB_DATA = ROOT_DIR / "PDBData" 
HMM = ROOT_DIR / "HMM"
TEMP_RESULTS = ROOT_DIR / "TemporaryResults"
ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES = ROOT_DIR / "Results"

for paths in (PDB_DATA, HMM, TEMP_RESULTS, ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES):
    try:
        paths.mkdir(exist_ok=False)
    except FileExistsError:
        print(f"Directory: {paths}\nAlready existed. Overwriting results.")
        shutil.rmtree(paths)
        paths.mkdir(exist_ok=False)
    else:
        print(f"Directory: {paths}\nNot found. Creating one.")

### GET ALL SEQUENCES ON PDB DATABASE IN FASTA FORMAT ###

#get all pdb sequences 
subprocess.run(["wget", "-P", PDB_DATA, "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"])
#decompress file
with open(PDB_DATA / "pdb_seqres.fasta", "wb") as outfile:
    subprocess.run(["gzip", "-dc", PDB_DATA / "pdb_seqres.txt.gz"], stdout=outfile)
#get all accs
with open(TEMP_RESULTS / "all_accs", "w") as outfile:
    subprocess.Popen(f"grep -oP '>\K\w+' {str(PDB_DATA)}/pdb_seqres.fasta", shell=True, stdout=outfile)

### RUN HMM_MODELS ###

#get heavy, lambda kappa chains HMM's from lyra githu
subprocess.run(["wget", "-P", HMM, "https://raw.githubusercontent.com/paolomarcatili/lyra/master/bcr_models/data/heavy.hmm"])
subprocess.run(["wget", "-P", HMM, "https://raw.githubusercontent.com/paolomarcatili/lyra/master/bcr_models/data/lambda.hmm"])
subprocess.run(["wget", "-P", HMM, "https://raw.githubusercontent.com/paolomarcatili/lyra/master/bcr_models/data/kappa.hmm"])

#run HMM search to heavy, lambda and kappa chains
subprocess.run(["hmmsearch", "--noali", "-E", str(HMM_expected_value), "-o", TEMP_RESULTS / "heavy_chains", HMM / "heavy.hmm", PDB_DATA / "pdb_seqres.fasta"])
subprocess.run(["hmmsearch", "--noali", "-E", str(HMM_expected_value), "-o", TEMP_RESULTS / "lambda_chains", HMM / "lambda.hmm", PDB_DATA / "pdb_seqres.fasta"])
subprocess.run(["hmmsearch", "--noali", "-E", str(HMM_expected_value), "-o", TEMP_RESULTS / "kappa_chains", HMM / "kappa.hmm", PDB_DATA / "pdb_seqres.fasta"])

#get chain accessions 
with open(TEMP_RESULTS / "heavy_accs", "w") as outfile:
    subprocess.Popen(f"grep -oP '>> \K\w+' {str(TEMP_RESULTS)}/heavy_chains", stdout = outfile, shell=True)
with open(TEMP_RESULTS / "lambda_accs", "w") as outfile:
    subprocess.Popen(f"grep -oP '>> \K\w+' {str(TEMP_RESULTS)}/lambda_chains", stdout = outfile, shell=True)
with open(TEMP_RESULTS / "kappa_accs", "w") as outfile:
    subprocess.Popen(f"grep -oP '>> \K\w+' {str(TEMP_RESULTS)}/kappa_chains", stdout = outfile, shell=True)

### GET ALL FASTA ACCs WHICH HAVE LAMBDA/KAPPA CHAIN, HEAVY CHAIN AND NON HMM HIT ###

all_chain_accs = list()
heavy_chain_accs = list()
lambda_chain_accs = list()
kappa_chain_accs = list()

#read all chain accs list
with open(TEMP_RESULTS / "all_accs", "r") as infile:
    for line in infile:
        all_chain_accs.append(line.strip())
#read all heavy chains to list
with open(TEMP_RESULTS / "heavy_accs", "r") as infile:
    for line in infile:
        heavy_chain_accs.append(line.strip())
#read all lambda chains to list
with open(TEMP_RESULTS / "lambda_accs", "r") as infile:
    for line in infile:
        lambda_chain_accs.append(line.strip())
#read all kappa_chains to list
with open(TEMP_RESULTS / "kappa_accs", "r") as infile:
    for line in infile:
        kappa_chain_accs.append(line.strip())

#convert everything to sets
all_chain_accs = set(all_chain_accs)
lambda_chain_accs = set(lambda_chain_accs)
kappa_chain_accs = set(kappa_chain_accs)
heavy_chain_accs = set(heavy_chain_accs)
#union of kappa and lambda chain accs = light chain accs
light_chain_accs = lambda_chain_accs | kappa_chain_accs
lambda_kappa_intersections = lambda_chain_accs & kappa_chain_accs
#union of heavy and light chains accs
heavy_and_light_chain_accs = light_chain_accs | heavy_chain_accs
#no hmm hit accs = accs not found in heavy chain + light chain list
no_hmm_hit_chain_accs = all_chain_accs - heavy_and_light_chain_accs
#heavy and light chain intersection 
heavy_light_chain_intersection = heavy_chain_accs & light_chain_accs

#some stats
num_all_chains = len(all_chain_accs)
num_light_chains = len(light_chain_accs)
num_lambda_chains = len(lambda_chain_accs)
num_kappa_chains = len(kappa_chain_accs)
num_heavy_chains = len(heavy_chain_accs)
num_heavy_and_light_chains = len(heavy_and_light_chain_accs)
num_no_hmm_hit_chains = len(no_hmm_hit_chain_accs)
num_heavy_light_intersections = len(heavy_light_chain_intersection)
num_lambda_kappa_intersections = len(lambda_kappa_intersections)
print(f"A total of {num_heavy_and_light_chains}/{num_all_chains} were assigned to light or heavy chains.")
print(f"Total chain accs assigned by as light chains: {num_light_chains}.")
print(f"Of these, the lambda HMM exclusively detected: {num_light_chains - num_kappa_chains}")
print(f"And the kappa HMM exclusively detected: {num_light_chains - num_lambda_chains}")
print(f"Light chain accs detected by both lambda and kappa HMM's: {num_lambda_kappa_intersections}\n")
print(f"Total chain accs assigned by heavy chain HMM: {num_heavy_chains}")
print(f"Total accs assigned as both heavy and light chain by HMMs: {num_heavy_light_intersections}.\nIf this number is very high, perhaps set the E-value to higher significance."\
        " Some cases exists perhaps due to some accesions not divided into heavy and light chains.")

### GET MAIN ACCESSIONS FOR AG-AB COMPLEXES ### 
##check for add. acc. chains (presumably antigens) for chain accs., assigned exclusively to either light or heavy chain.
light_chain_accs = light_chain_accs - heavy_light_chain_intersection 
#print(light_chain_accs) 
heavy_chain_accs = heavy_chain_accs - heavy_light_chain_intersection
light_accs = set([acc.split("_")[0] for acc in light_chain_accs])
heavy_accs = set([acc.split("_")[0] for acc in heavy_chain_accs])
no_hmm_hit_accs = set([acc.split("_")[0] for acc in no_hmm_hit_chain_accs])

#antibody-antigen complexes are accesions which contained a light, heavy and no hmm hit chain 
ag_ab_accs = light_accs & heavy_accs & no_hmm_hit_accs
print(f"Antibody-antigens complexes where light, heavy chain were exclusively assigned to heavy or light chain HMMs: {len(ag_ab_accs)}")
with open(ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES / "AgAbComplexCandidates", "w") as outfile:
    for ag_ab_acc in ag_ab_accs:
        outfile.write(f"{ag_ab_acc}\n")

#check for add. acc chains (presumably antigens) for chain accs. assigned by both light and heavy chain HMM.
heavy_light_intersection = set([acc.split("_")[0] for acc in heavy_light_chain_intersection])
ag_ab_accs_extra = heavy_light_intersection & no_hmm_hit_accs
print(f"Antibody-antigens complexes where sequences were assigned to both heavy and light chain HMMs: {len(ag_ab_accs_extra)}")
with open(ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES / "AgAbComplexExtraCandidates", "w") as outfile:
    for ag_ab_acc in ag_ab_accs_extra:
        outfile.write(f"{ag_ab_acc}\n")

### SAVE RESULTS NEEDED FOR LATER GETTING RESIDUE CONTACTS ###

#save light chain
lc_outfile = open(ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES / "LightChainAcc.pickle", "wb")
pickle.dump(light_chain_accs, lc_outfile)
#save heavy chain
hc_outfile = open(ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES / "HeavyChainAcc.pickle", "wb")
pickle.dump(heavy_chain_accs, hc_outfile)
#save light and heavy chain intersection
lc_hc_chain_intersection_outfile = open(ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES / "LightHeavyChainIntersectionAcc.pickle", "wb")
pickle.dump(heavy_light_chain_intersection, lc_hc_chain_intersection_outfile)
#save no_hmm_hit_outfile
no_hmm_hit_outfile = open(ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES / "NoHMMHitChainAcc.pickle", "wb")
pickle.dump(no_hmm_hit_chain_accs, no_hmm_hit_outfile)
#save all_chain_accs (having this is nice as a sanity lookup table for later).
#Because when downloading PDB files there will be multiple models. And the chains in these other models will have different chain accessions
all_chain_accs_outfile = open(ANTIGEN_ANTIBODY_COMPLEX_CANDIDATES / "AllChainAccs.pickle", "wb")

pickle.dump(all_chain_accs, all_chain_accs_outfile)

#close files
lc_outfile.close()
hc_outfile.close()
lc_hc_chain_intersection_outfile.close()
no_hmm_hit_outfile.close()
all_chain_accs_outfile.close()

### CLEAN UP ###
if clean_up:
    shutil.rmtree(PDB_DATA)
    shutil.rmtree(HMM)
    shutil.rmtree(TEMP_RESULTS)
