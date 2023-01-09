# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 14:39:53 2022

@author: Max
"""
import subprocess
import tempfile
import numpy as np
import argparse
import os
import glob
import scipy.cluster.hierarchy as hierarchy
from Bio.Align import substitution_matrices
import csv

def parse_peptide_sequences(pdb_files, chain_length):
    """Parses the peptide sequences from the specified PDB files and returns them as a dictionary with the file name as the key and the peptide sequence as the value.

    Arguments:
        pdb_files {list} -- A list of paths to the PDB files
        chain_length {int} -- The length of the protein chain containing the peptides

    Returns:
        dict -- A dictionary with the file name as the key and the peptide sequence as the value
    """

    peptide_sequences = {}

    # parse the peptide sequences from the PDB files
    for pdb_file in pdb_files:
        with open(pdb_file) as f:
            # parse the lines using column numbers
            lines = [line for line in f if line.startswith("ATOM")]
            res_types = [line[17:20].strip() for line in lines]
            chains = [line[21] for line in lines]
            residue_sequences = [line[22:26].strip() for line in lines]

            # check if the atoms are part of the chain with the specified length
            res_types_by_chain = {}
            for rt, c, rs in zip(res_types, chains, residue_sequences):
                if c not in res_types_by_chain:
                    res_types_by_chain[c] = {}
                if rs not in res_types_by_chain[c]:
                    res_types_by_chain[c][rs] = rt
            
            # find the chain with the specified length
            target_chain = None
            for c, rts in res_types_by_chain.items():
                if len(rts) == chain_length:
                    target_chain = c
                    break

            # store the peptide sequence as a string
            if target_chain is not None:
                peptide_sequence = "".join(res_types_by_chain[target_chain].values())
            else:
                peptide_sequence = ""
        
        # Get the file name with the extension
        file_name_with_ext = os.path.basename(pdb_file)
        # Split the file name and extension
        file_name, file_ext = os.path.splitext(file_name_with_ext)
        # Add the file name and peptide sequence to the dictionary only if the peptide sequence is not empty
        if peptide_sequence:
            peptide_sequence = convert_amino_acid_identifiers(peptide_sequence)
            peptide_sequences[file_name] = peptide_sequence

    return peptide_sequences

def convert_amino_acid_identifiers(peptide_sequence):
    # Create a dictionary to map 3 letter amino acid identifiers to single letter identifiers
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    # Initialize an empty result string
    result = ''

    # Iterate through the input string and convert the 3 letter amino acid identifiers to single letter identifiers
    for i in range(0, len(peptide_sequence), 3):
        aa = peptide_sequence[i:i+3]
        result += aa_map.get(aa, 'X')

    # Return the result
    return result

def cluster_peptides_with_cdhit(peptides, n_clusters):

    # open the file in write mode
    with open('peptides.faa', 'w') as f:
        # write the peptide sequences in fasta format to the file
        for i, peptide in enumerate(peptides):
            f.write(f">Peptide_{i+1}\n{peptide}\n")

    # run CD-HIT using the file you just created
    cmd = f"cd-hit -i peptides.faa -o peptides_clusters.txt -c 1.0 -n {n_clusters} -M 16000"
    subprocess.run(cmd.split())

    # parse the cluster representatives from the CD-HIT output file
    representatives = []
    with open('peptides_clusters.txt') as clstr_file:
        for line in clstr_file:
            if line.startswith(">"):
                representatives.append(line.split()[2])

    # return the cluster representatives
    return representatives



def compute_pam250_distance(pam250, pep1, pep2):
    """Computes the distance between two peptide sequences using the PAM250 matrix.

    Arguments:
        pep1 {str} -- The first peptide sequence
        pep2 {str} -- The second peptide sequence

    Returns:
        float -- The distance between the two peptide sequences
    """

    # compute the distance between the two peptides
    distance = sum(pam250[aa1][aa2] for aa1, aa2 in zip(pep1, pep2))

    return distance

def cluster_peptides_pam250(peptides, n_clusters):
    # Fetch the PAM250 matrix using the Bio.Align.SubstitutionMatrix class
    pam250 =  substitution_matrices.load("PAM250")
 
    # Compute the pairwise distances between peptides using the PAM250 matrix
    distances = np.zeros((len(peptides), len(peptides)))
    for i, p1 in enumerate(peptides):
        for j, p2 in enumerate(peptides):
            distances[i, j] = compute_pam250_distance(pam250, p1, p2)

    # Cluster the peptides using hierarchical clustering
    linkage = hierarchy.linkage(distances, method="single")
    clusters = hierarchy.fcluster(linkage, n_clusters, criterion="maxclust")

    # Select the cluster representatives
    representatives = []
    for i in range(1, n_clusters + 1):
        cluster_peptides = [pep for pep, c in zip(peptides, clusters) if c == i]
        representative = select_most_dissimilar_peptide(pam250, cluster_peptides)
        representatives.append(representative)

    return representatives

def select_most_dissimilar_peptide(peptides, pam250):
    """Selects the most dissimilar peptide from a list of peptides.

    Arguments:
        peptides {list} -- A list of peptide sequences
        pam250 {ndarray} -- The PAM250 matrix

    Returns:
        str -- The most dissimilar peptide from the input list
    """
    # find the pair of peptides with the greatest distance
    max_distance = -float("inf")
    max_i, max_j = -1, -1
    for i in range(len(peptides)):
        for j in range(len(peptides)):
            if i != j:
                distance = compute_pam250_distance(pam250, peptides[i], peptides[j])
                if distance > max_distance:
                    max_distance = distance
                    max_i, max_j = i, j

    # return the most dissimilar peptide
    if max_distance == -float("inf"):
        return peptides[0]
    else:
        return peptides[max_i] if max_i > max_j else peptides[max_j]

def write_peptides_to_csv(peptides, peptide_IDs, file_name):
    with open(file_name, 'w', newline='') as csv_file:
        # Create a CSV writer object
        writer = csv.writer(csv_file)

        # Write the headers to the CSV file
        writer.writerow(['peptide_ID', 'peptide'])

        # Zip the peptide IDs and peptide sequences into a list of tuples
        peptide_data = zip(peptide_IDs, peptides)

        # Write the data to the CSV file
        for peptide_ID, peptide in peptide_data:
            writer.writerow([peptide_ID, peptide])

def get_file_names(file_paths):
    file_names = []
    for file_path in file_paths:
        # Get the file name with the extension
        file_name_with_ext = os.path.basename(file_path)
        # Split the file name and extension
        file_name, file_ext = os.path.splitext(file_name_with_ext)
        # Append the file name to the list
        file_names.append(file_name)
    return file_names

if __name__ == "__main__":
    # parse the path to the directory containing PDB files from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_dir", help="Path to the directory containing PDB files")
    args = parser.parse_args()

    # get the list of PDB files in the specified directory
    pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))
    print(pdb_files)

    # parse the peptide sequences from the PDB files
    peptides = parse_peptide_sequences(pdb_files, 9)
    print(peptides)

    # get peptide pdb IDs
    peptide_IDs = get_file_names(pdb_files)

    # convert residues to single letter ID
    peptides = convert_amino_acid_identifiers(peptides)
    
    # Write peptides to csv files
    write_peptides_to_csv(peptides, peptide_IDs, 'cluster_peptides.csv')

    # cluster the peptides with CD-HIT
    cd_hit_representatives = cluster_peptides_with_cdhit(peptides, n_clusters=10)
    print(cd_hit_representatives)
    
    #cluster peptides with Pam250
    Pam250_representatives = cluster_peptides_pam250(peptides, n_clusters=10)
    print(Pam250_representatives)
