# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:10:01 2023

@author: Max
"""
import os
import time
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from tqdm import tqdm
import re
from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

def parse_pdb_information(pdb_file):
    information_lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                break
            else:
                information_lines.append(line.strip())
    return information_lines

def chain_information(information_lines):
    information_string = '\n'.join(information_lines)
    chain_lines = re.findall(r"COMPND\s+\d+\s+CHAIN:\s+([A-Z,\s]+)", information_string)
    molecule_lines = re.findall(r"COMPND\s+\d+\s+MOLECULE:\s+(.+)", information_string)
    
    parsed_chains = []
    for chain, molecule in zip(chain_lines, molecule_lines):
        chains = []
        if any(substring.lower() in molecule.lower() for substring in ["peptide", "epitope", "class I", "Beta"]):
            # extract the chains from the chain line
            for chain_list in chain.split(", "):
                for chain_id in chain_list.strip().split():
                    if not chains or ord(chain_id) == ord(chains[-1]) + 1:
                        chains.append(chain_id)
                    else:
                        break  # ignore this chain list
            separator = ', '
            wanted_chains = separator.join(chains)
            line = re.search(r"COMPND\s+\d+\s+CHAIN:\s+{}".format(chain), information_string).group(0)
            new_line = line.replace(chain, wanted_chains)
            information_string = information_string.replace(line, new_line)
        else:
            line = re.search(r"COMPND\s+\d+\s+CHAIN:\s+{}".format(chain), information_string).group(0)
            # remove the line from the information string
            information_string = information_string.replace(line, '')
            
        # append the parsed chains to the list
        parsed_chains.append(','.join(chains))

    # return the parsed chains list and the updated information string
    return parsed_chains, information_string

def pdb_cleaner(pdb_file):
    # Create a PDBFixer object
    fixer = PDBFixer(pdb_file)
    fixer.removeHeterogens(True)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    try:
        fixer.addMissingAtoms()
    except KeyError:
        print(f"Error: failed to add missing atoms to {pdb_file}. Skipping...")

    # Set the residue names to a standard naming convention
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Write the modified PDB file to the output file
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'))

if __name__ == "__main__":
    # Get directory path from user
    directory_path = input("Enter directory path containing PDB files: ")

    # Find all PDB files in directory
    pdb_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.pdb')]
    # Initialize tqdm progress bar
    progress_bar = tqdm(total=len(pdb_files), desc="Processing PDB files")
    
    start_time = time.time()
    for pdb_file in pdb_files:
        pdb_cleaner(pdb_file)
        progress_bar.update(1)
    end_time = time.time()
    progress_bar.close()
    
    print(f"Finished processing {len(pdb_files)} files in {end_time - start_time:.2f} seconds") 