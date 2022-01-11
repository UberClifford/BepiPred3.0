### IMPORTS ###
import pickle
from collections import Counter
from pathlib import Path
import sys
import random
import subprocess
import argparse

### STATIC PATHS AND VARIABLES ###
ROOT_DIR = Path.cwd()

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Takes pickle formated files, single chain ab-ag complexes, converts to probable light and heavy chain. And joins with ag-ab light and heavy chain separated files.")
parser.add_argument("-i", required=True, action="store", dest="infile", type=Path, help="Input file containing antibody-antigen interactions. Should be pickle formatted file")
parser.add_argument("-ab", required=False, default=0.9, type=float, action="store", dest="ab_clusID", help="Sequence identity to cluster antibodies at.")
parser.add_argument("-ag", required=False, default=0.7, type=float, action="store", dest="ag_clusID", help="Sequence identity to cluster antigens at.")
parser.add_argument("-t", required=True, action="store", dest="target_dir", type=Path, help="Target directory. Will write results to this directory.")

args = parser.parse_args()
infile = args.infile
target_dir = args.target_dir
ab_clusID = args.ab_clusID
ag_clusID = args.ag_clusID
ABClus_filename = f"ABClus{int(ab_clusID*100)}ID"
AGClus_filename = f"AGClus{int(ag_clusID*100)}ID"

antibody_target = target_dir / ABClus_filename
antigen_target = target_dir / AGClus_filename

if not infile.is_file():
    sys.exit(f"The input file was not a valid file.\nInput file given:{infile.name})")

if not target_dir.is_dir():
    target_dir.mkdir(exist_ok=False, parents=True)
if not antibody_target.is_dir():
    antibody_target.mkdir(exist_ok=False, parents=True)
if not antigen_target.is_dir():
    antigen_target.mkdir(exist_ok=False, parents=True)

ag_ab_interactions = pickle.load(open(infile, "rb"))

### FUNCTIONS ###

def data_to_fasta_format(pdb_accs, outfile_path):

    with open(outfile_path, "w") as outfile:
        output = str()
        for pdb_acc in pdb_accs:
            output += f">{pdb_acc[0]}\n{pdb_acc[1]}\n"

        output = output[:-1]
        outfile.write(output)

def get_sorted_cluster_sizes(all_cluster_acc_and_seq):
    cluster_sizes = list()
    for cluster_rep in all_cluster_acc_and_seq.keys():
        cluster_sizes.append( (cluster_rep, len(all_cluster_acc_and_seq[cluster_rep])) )
    cluster_sizes = sorted(cluster_sizes, key=lambda pair: pair[1], reverse=True)
    
    return cluster_sizes

def get_cluster_accs_from_mmseq_output(infile):

    infile = open(infile, "r")
    readfile = infile.readlines()
    infile.close()
    
    cluster_accs = list()
    cluster_sequences = list()
    all_cluster_acc_and_seqs = dict()
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

    #if last cluster only has one sequence
    if len(cluster_accs) == 1:
        all_cluster_acc_and_seqs[cluster_accs[0]] = [ (cluster_accs[0], cluster_sequence) ]
    #if last cluster has multiple sequences
    else:
        #append last sequence 
        cluster_sequences.append(cluster_sequence)
        all_cluster_acc_and_seqs[cluster_accs[0]] = list( zip(cluster_accs, cluster_sequences) )

    return all_cluster_acc_and_seqs, total_acc_count

def get_antibodies_and_antigens(ag_ab_interactions):
    antibodies = list()
    antigens = list()
    for ag_ab_interaction in ag_ab_interactions:
        light_chain = ag_ab_interaction[1]
        heavy_chain = ag_ab_interaction[3]
        antibody = light_chain + heavy_chain
        antigen = ag_ab_interaction[5]
        unique_agab_complex_id = ag_ab_interaction[6]

        antigens.append((unique_agab_complex_id, antigen))
        antibodies.append((unique_agab_complex_id, antibody))

    return antigens, antibodies



def cluster_for_accs(cluster_dict, check_against_ids):

    cluster_chunks = list()
    cluster_reps = list()
    for cluster_rep in cluster_dict.keys():
        cluster = cluster_dict[cluster_rep]
        cluster_ids = [c[0] for c in cluster]
        cluster_ids_check = set(cluster_ids)

        #check if antigens belong to antibodies, collect it and its entire cluster
        intersection = cluster_ids_check.intersection(check_against_ids)
        if intersection:
            cluster_reps.append(cluster_rep)
            cluster_chunks.extend(cluster_ids)

    cluster_chunks = set(cluster_chunks)

    return cluster_chunks, cluster_reps

def remove_clusters_from_cluster(cluster_dict, keys_to_del):

    for key_to_del in keys_to_del:
        del cluster_dict[key_to_del]

    return cluster_dict

def extract_test_and_train_val(abag_chunks):

    #extract test set chunks
    test_accs = list()
    train_val_accs = list()
    cluster_count = 0
    for abag_chunk in abag_chunks:
        check = [cluster_acc for cluster_acc in abag_chunk if cluster_acc.startswith("Test")]
        cluster_size = len(abag_chunk)

        #avoiding really big chunk to be test set. Doing this, we miss complexes  ['Test_ComplexID_349', 'Test_ComplexID_695', 'Test_ComplexID_762', 'Test_ComplexID_1812']
        if check and cluster_size <= 200:
            cluster_count += 1
            cluster_accs = list(abag_chunk)
            test_accs.extend(cluster_accs)
    
        else:
            train_val_accs.append(abag_chunk)

    print(f"Antibody-antigen clusters in test set {cluster_count}")

    return train_val_accs, test_accs


def cluster_accs_to_partitions(abag_chunks):
    i = 1
    partition1 = list()
    partition2 = list()
    partition3 = list()
    partition4 = list()
    partition5 = list()

    train_val, test = extract_test_and_train_val(abag_chunks)
    #sort cluster chunks from large to small 
    chunk_sizes = [len(chunk) for chunk in train_val]
    train_val = [abag_chunk for _, abag_chunk in sorted(zip(chunk_sizes, train_val), key=lambda pair: pair[0], reverse=True)]

    for tv in train_val:
        cluster_accs = list(tv)

        if i == 1:
            partition1.extend(cluster_accs)
        elif i == 2:
            partition2.extend(cluster_accs)
        elif i == 3:
            partition3.extend(cluster_accs)
        elif i == 4:
            partition4.extend(cluster_accs)
        elif i == 5:
            partition5.extend(cluster_accs)

        i += 1
        if i == 6:
            i = 1

    return partition1, partition2, partition3, partition4, partition5, test


### MAIN ###

#write antibodies (light and heavy chain concatenated) and antigens to fastas files 
antigens, antibodies = get_antibodies_and_antigens(ag_ab_interactions)
data_to_fasta_format(antibodies, "antibodies.fasta")
data_to_fasta_format(antigens, "antigens.fasta")

#Do mmseqs 90% and 70% sequence identity clustering for antibody and antigen respectively
subprocess.run(["mmseqs", "easy-cluster", "antibodies.fasta", antibody_target / ABClus_filename, antibody_target / "tmp", "--min-seq-id", str(ab_clusID), "-c", "0.8", "--cov-mode", "1"])
subprocess.run(["mmseqs", "easy-cluster", "antigens.fasta", antigen_target / AGClus_filename, antigen_target / "tmp", "--min-seq-id", str(ag_clusID), "-c", "0.8", "--cov-mode", "1"])

ab_mmseq_cluster_file = antibody_target / f"{ABClus_filename}_all_seqs.fasta" 
ag_mmseq_cluster_file = antigen_target / f"{AGClus_filename}_all_seqs.fasta"

#First extract anything to do with test set
all_ab_clusters, total_acc_count = get_cluster_accs_from_mmseq_output(ab_mmseq_cluster_file)
all_ag_clusters, total_acc_count = get_cluster_accs_from_mmseq_output(ag_mmseq_cluster_file)

abag_cluster_chunks = list()

### MAIN ###
while all_ab_clusters:
    print(f"Making (AB{int(ab_clusID*100)}ID, AG{int(ag_clusID*100)}ID) cluster chunk\n")
    #randomly pick a 90% seq identity AB cluster representative
    ab_cluster_rep = random.choice(list(all_ab_clusters.keys()))
    ab_cluster = all_ab_clusters[ab_cluster_rep]
    #get all sequence identifiers for cluster
    ab_cluster_ids = [c[0] for c in ab_cluster]
    abag_cluster_chunk = set(ab_cluster_ids)

    #find corresponding antigens
    ab_keys_to_del = [ab_cluster_rep]
    ag_keys_to_del = []
    chunk_done = False
    print(f"Initial chunk size: {len(abag_cluster_chunk)}")

    while not chunk_done:
        #get size of current cluster chunk
        chunk_size = len(abag_cluster_chunk)
        #propagate chunk
        chunk_before = abag_cluster_chunk
        abag_cluster_chunk, new_ag_keys_to_del = cluster_for_accs(all_ag_clusters, abag_cluster_chunk)
        print(f"{len(chunk_before)}-->{len(abag_cluster_chunk)} when checking corresponding antigen clusters")
        chunk_before = abag_cluster_chunk
        abag_cluster_chunk, new_ab_keys_to_del = cluster_for_accs(all_ab_clusters, abag_cluster_chunk)
        print(f"{len(chunk_before)}-->{len(abag_cluster_chunk)} when checking corresponding antibody clusters")

        ab_keys_to_del.append(new_ab_keys_to_del)
        ag_keys_to_del.append(new_ag_keys_to_del)

        if chunk_size > len(abag_cluster_chunk):
            print("Cluster chunk size decreased, something is wrong")
        elif chunk_size < len(abag_cluster_chunk):
            print(f"Cluster chunk size increased, {chunk_size}->{len(abag_cluster_chunk)}")
        elif chunk_size == len(abag_cluster_chunk):
            chunk_done = True
            print(f"No new AG-AB complexes added to chunk. Stopping to collect complex IDs. Chunk size: {len(abag_cluster_chunk)}")
        else:
            print("Something is wrong")
    
    #Collecting cluster chunk
    new_ag_keys_to_del = set(new_ag_keys_to_del)
    new_ab_keys_to_del = set(new_ab_keys_to_del)
    abag_cluster_chunks.append(tuple(abag_cluster_chunk)) 
    #remove collected complex IDS 
    all_ab_clusters = remove_clusters_from_cluster(all_ab_clusters, new_ab_keys_to_del)
    all_ag_clusters = remove_clusters_from_cluster(all_ag_clusters, new_ag_keys_to_del)

#collected chunks should be same size as input AG-AB complex.
chunk_sizes = [len(chunk) for chunk in abag_cluster_chunks]
num_of_collected_agab_complexes = sum(chunk_sizes) 
chunk_size_distribution = Counter(chunk_sizes)
chunk_size_distribution = sorted([(chunk_key, chunk_size_distribution[chunk_key]) for chunk_key in chunk_size_distribution.keys()], key=lambda pair: pair[1], reverse=True)

print(f"Antibody-antigen complexes in input file {len(ag_ab_interactions)}")
print(f"Total collected antigen-antibody complexes {num_of_collected_agab_complexes}")
if num_of_collected_agab_complexes != len(ag_ab_interactions):
    print(f"Collected clustered antibody-antigen complexes not the same as input antibody-antigen complex. Something is wrong!")
print(f"Amount of (AB{int(ab_clusID*100)}ID, AG{int(ag_clusID*100)}ID) of clusters: {len(abag_cluster_chunks)}")


#sort cluster chunks from large to small 
#sorted_abag_chunks = [abag_chunk for _, abag_chunk in sorted(zip(chunk_sizes, abag_cluster_chunks), key=lambda pair: pair[0], reverse=True)]
#this chunk was very large, therefore interesting
#very_large_chunk = set(sorted_abag_chunks[0])



partition1, partition2, partition3, partition4, partition5, test = cluster_accs_to_partitions(abag_cluster_chunks)
partition1_set = set(partition1)
partition2_set = set(partition2)
partition3_set = set(partition3)
partition4_set = set(partition4)
partition5_set = set(partition5)
test_set = set(test)
#
print(f"Partition size 1: {len(partition1_set)}\nPartition size 2: {len(partition2_set)}\nPartition size 3: {len(partition3_set)}\nPartition size 4: {len(partition4_set)}\nPartition size 5: {len(partition5_set)}")
print(f"Test set size: {len(test_set)}")
comparison = partition1_set.intersection(partition1_set, partition2_set, partition3_set, partition4_set, test_set)
print(f"Check that all partitions are dont contain same complex IDs: {comparison}. Should be an empty set")

print("Overall chunk size distribution")
for chunk in chunk_size_distribution:
    chunk_size = chunk[0]
    chunk_size_count = chunk[1]
    print(f"Chunk size of {chunk_size} counts: {chunk_size_count}")

### Save partitions as 5 pickle files ### 

partition1 = list()
partition2 = list()
partition3 = list()
partition4 = list()
partition5 = list()
test = list()
test_fasta = str()
partition1_fasta = str()
partition2_fasta= str()
partition3_fasta = str()
partition4_fasta = str()
partition5_fasta = str()

for ag_ab_interaction in ag_ab_interactions:
    unique_agab_complex_id = ag_ab_interaction[6]

    if unique_agab_complex_id in partition1_set:
        partition1.append(ag_ab_interaction)
        partition1_fasta  += f">{ag_ab_interaction[4]}\n{ag_ab_interaction[5]}\n"
    elif unique_agab_complex_id in partition2_set:
        partition2.append(ag_ab_interaction)
        partition2_fasta  += f">{ag_ab_interaction[4]}\n{ag_ab_interaction[5]}\n"
    elif unique_agab_complex_id in partition3_set:
        partition3.append(ag_ab_interaction)
        partition3_fasta  += f">{ag_ab_interaction[4]}\n{ag_ab_interaction[5]}\n"
    elif unique_agab_complex_id in partition4_set:
        partition4.append(ag_ab_interaction)
        partition4_fasta  += f">{ag_ab_interaction[4]}\n{ag_ab_interaction[5]}\n"
    elif unique_agab_complex_id in partition5_set:
        partition5.append(ag_ab_interaction)
        partition5_fasta  += f">{ag_ab_interaction[4]}\n{ag_ab_interaction[5]}\n"
    elif unique_agab_complex_id in test_set:
        test.append(ag_ab_interaction)
        test_fasta += f">{ag_ab_interaction[4]}\n{ag_ab_interaction[5]}\n"

    else:
        print("Something is wrong!")

with open(target_dir / "Fold1_antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(partition1, outfile)
with open(target_dir / "Fold2_antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(partition2, outfile)
with open(target_dir / "Fold3_antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(partition3, outfile)
with open(target_dir / "Fold4_antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(partition4, outfile)
with open(target_dir / "Fold5_antigens_antibody_interactions.pickle", "wb") as outfile:
    pickle.dump(partition5, outfile)
with open(target_dir / "test.pickle", "wb") as outfile:
    pickle.dump(test, outfile)

#also partition and test antigens as fasta files
with open(target_dir / "partition1.fasta", "w") as outfile:
    outfile.write(partition1_fasta[:-1])
with open(target_dir / "partition2.fasta", "w") as outfile:
    outfile.write(partition2_fasta[:-1])
with open(target_dir / "partition3.fasta", "w") as outfile:
    outfile.write(partition3_fasta[:-1])
with open(target_dir / "partition4.fasta", "w") as outfile:
    outfile.write(partition4_fasta[:-1])
with open(target_dir / "partition5.fasta", "w") as outfile:
    outfile.write(partition5_fasta[:-1])

with open(target_dir / "test.fasta", "w") as outfile:
    outfile.write(test_fasta[:-1])


