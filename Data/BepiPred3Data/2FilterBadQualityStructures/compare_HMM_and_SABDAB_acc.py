## IMPORTS ###
import subprocess
from collections import Counter
from pathlib import Path
import argparse
import shutil
import pickle

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Compare two files containing accesions.")
parser.add_argument("-f1", required=True, action="store", dest="hmm_accs_path", help="Accesions resulting from light, heavy chain HMM profiling.")
parser.add_argument("-f2", required=True, action="store", dest="hmm_extra_accs_path", help="Extra accesions resulting from light, heavy chain HMM profiling. These extra candidates are typically instances\
        where the light and heavy are not separated into indiviudal chains, but are combined in one chain.")
parser.add_argument("-d", required=True, action="store", dest="sabdab_accs_path", help="Directory containing pdb files resulting from SABDAB database quality query.")

args = parser.parse_args()
hmm_accs_path = args.hmm_accs_path
hmm_extra_accs_path = args.hmm_extra_accs_path
sabdab_accs_path = args.sabdab_accs_path

ROOT_DIR = Path.cwd()
target_dir = ROOT_DIR / "Results"
sabdab_accs_path = Path(sabdab_accs_path)
sabdab_accs = [pdb_file.stem for pdb_file in sabdab_accs_path.glob("*.pdb")]

try:
    target_dir.mkdir(exist_ok=False)
except FileExistsError:
    print(f"Directory: {target_dir}\nAlready existed. Overwriting results.")
    shutil.rmtree(target_dir)
    target_dir.mkdir(exist_ok=False)
else:
    print(f"Directory: {target_dir}\nNot found. Creating one.")
 
hmm_accs = list()
hmm_extra_accs = list()

#read all hmm accs to list
with open(hmm_accs_path, "r") as infile:
    for line in infile:
        hmm_accs.append(line.strip())
#read all hmm extra accs to list
with open(hmm_extra_accs_path, "r") as infile:
    for line in infile:
        hmm_extra_accs.append(line.strip())

#convert to sets
hmm_accs = set(hmm_accs)
hmm_extra_accs = set(hmm_extra_accs)
sabdab_accs = set(sabdab_accs)

num_hmm_accs = len(hmm_accs)
num_hmm_extra_accs = len(hmm_extra_accs)
num_sabdab_accs = len(sabdab_accs)
print(f"Accesions in hmm profiling file {num_hmm_accs}")
print(f"Accesions in extra hmm profiling file {num_hmm_extra_accs}")
print(f"Accesions in sabdab acc file {num_sabdab_accs}")


#get intersection of sets
intersection = hmm_accs & sabdab_accs
extra_intersection = hmm_extra_accs & sabdab_accs
#In the case where we compare SABDABDatabase and HMM profiling results, 
#these are accesions which were identified as AGAB complex HMM profiling 
#and also seem to be of sufficient quality by the SABDABDatabase.
print(f"Number accesions shared between hmm profiling file and sabdab acc file {len(intersection)}")
print(f"Number accesions shared between extra hmm profiling file and sabdab acc file {len(extra_intersection)}")

intersection_outfile = open(target_dir / "hmm_sabdab_intersection.pickle", "wb")
extra_intersection_outfile = open(target_dir / "hmm_sabdab_extra_intersection.pickle", "wb")

hmm_set_outfile = open(target_dir / "hmm_set.pickle", "wb")
hmm_extra_set_outfile = open(target_dir / "hmm_extra_set.pickle", "wb")
sabdab_set_outfile = open(target_dir / "sabdab_set.pickle", "wb")

pickle.dump(intersection, intersection_outfile)
pickle.dump(extra_intersection, extra_intersection_outfile)

pickle.dump(hmm_accs, hmm_set_outfile)
pickle.dump(hmm_extra_accs, hmm_extra_set_outfile)
pickle.dump(sabdab_accs, sabdab_set_outfile)

#close files
intersection_outfile.close()
extra_intersection_outfile.close()

hmm_set_outfile.close()
hmm_extra_set_outfile.close()
sabdab_set_outfile.close()
