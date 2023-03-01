# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:10:01 2023

@author: Max
"""
import argparse
import os
import time
from pdbfixer import PDBFixer
from tqdm import tqdm
import re
from openmm.app import PDBFile
from pathlib import Path

def parse_pdb_information(pdb_file):
    """
    Parses the PDB file and extracts the header information from it.

    Args:
        pdb_file (str): The path to the PDB file to be parsed.

    Returns:
        information_lines (list): A list of strings containing the header information extracted from the PDB file.
    """
    information_lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                break
            else:
                information_lines.append(line)
    return information_lines

def chain_information(information_lines, pdb_file):
    """
    Extracts the chain information from the header of a PDB file.

    Args:
        information_lines (list): A list of strings containing the header information of the PDB file.
        pdb_file (str): The name of the PDB file being processed.

    Returns:
        parsed_chains (list): A list of characters representing the chain IDs present in the PDB file.
        information_string (str): A string containing the updated header information of the PDB file after parsing the chain information.
    """
    
    information_lines = [line.strip() for line in information_lines]
    information_string = '\n'.join(information_lines)
    chain_lines = re.findall(r"COMPND\s+\d+\s+CHAIN:\s+([a-zA-Z,\s]+);?\n", information_string)
    molecule_lines = re.findall(r"COMPND\s+\d+\s+MOLECULE:\s+(.+)", information_string)
    pdb_ID = Path(pdb_file).resolve().stem.lower()
    peptide_chains = re.findall(r"IMGT\sreceptor\sdescription.+?Peptide.+?Chain\sID.+?{}_(\w)".format(pdb_ID), information_string, re.MULTILINE | re.DOTALL)

    separator = ', '
    parsed_chains = []
    last_chain_id = None

    for chain, molecule in zip(chain_lines, molecule_lines):
        chain_IDs = chain.split(", ")
        chains = []
        
        if any(substring.lower() in molecule.lower() for substring in ["peptide", "epitope", "class I", "globulin", "hla"]) or any(peptide_chain in chain for peptide_chain in peptide_chains):
            try:
                for chain_id in chain_IDs:
                    # Remove any possible white spaces
                    chain_id = chain_id.strip()
                    if not chains or ord(chain_id) == ord(last_chain_id) + 1:
                        chains.append(chain_id)
                        last_chain_id = chain_id
                    else:
                        break
                line = re.search(r"COMPND\s+\d+\s+CHAIN:\s+{}".format(chain), information_string).group(0)
                if chains:
                    new_line = line.replace(chain, separator.join(chains))
                    information_string = information_string.replace(line, new_line)
                    parsed_chains.extend(chains)
                else:
                    information_string = information_string.replace(line, 'COMPND')
            except:
                print("Error processing chain in file {}: {}".format(pdb_file, chain))
                print("chain_id that caused the error: {}".format(chain_id))
        else:
            line = re.search(r"COMPND\s+\d+\s+CHAIN:\s+{}".format(chain), information_string).group(0)
            information_string = information_string.replace(line, 'COMPND')

    missing_chains_file = "pdb_files_with_missing_chains.txt"
    if len(parsed_chains) < 3:
        # Check if the line is already present in the file
        line_to_write = "{} did not contain the minimum expected amount of parsed chains\n".format(pdb_file)
        with open(missing_chains_file, "a+") as f:
            f.seek(0)
            if line_to_write in f.read():
                print(f"{line_to_write.strip()} already present in {missing_chains_file}. Skipping...")
            else:
                f.write(line_to_write)
                print(f"{line_to_write.strip()} written to {missing_chains_file}.")

    return parsed_chains, information_string

def pdb_cleaner(pdb_file, parsed_chains, information_lines=None):
    """
    Cleans the PDB file by removing heterogens, finding missing residues and atoms, and adding missing heavy atoms using PDBFixer.

    Args:
        pdb_file (str): The path to the PDB file to be cleaned.
        parsed_chains (list): A list of characters representing the chain IDs present in the PDB file.
        information_lines (list): A list of strings representing the header information.

    Returns:
        
    """
    # Create a PDBFixer object
    fixer = PDBFixer(pdb_file)

    # Remove heterogens
    fixer.removeHeterogens(True)

    # Find missing residues and atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()

    # Add missing atoms
    try:
        fixer.addMissingAtoms()
    except KeyError:
        print(f"Error: failed to add missing atoms to {pdb_file}. Skipping...")

    # Set the residue names to a standard naming convention
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Remove unwanted chains
    fixer.removeChains([i for i, chain in enumerate(fixer.topology.chains()) if chain.id not in parsed_chains])

    # Write the modified PDB file to the output file
    temp_pdb_file = "temp.pdb"
    PDBFile.writeFile(fixer.topology, fixer.positions, open(temp_pdb_file, 'w'))

    # Combine the header information and the temporary PDB file
    with open(temp_pdb_file, 'r') as f:
        temp_pdb_lines = '\n'.join(f.readlines())
    with open(pdb_file, 'w') as f:
        f.writelines(information_string + "\n" +  temp_pdb_lines)

    # Delete the temporary file
    os.remove(temp_pdb_file)

def is_pdb_cleaned(pdb_file):
    """
    Check if the PDB file has already been cleaned.

    Parameters:
        pdb_file (str): The path to the PDB file to check.

    Returns:
        bool: True if the file has been cleaned (contains a line starting with "REMARK" and containing
              "CREATED WITH OPENMM"), False otherwise.
    """
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('REMARK') and 'CREATED WITH OPENMM' in line:
                return True
    return False



if __name__ == "__main__":
    # Get directory path from user
    # Parse the input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory_path", help="The path containing PDB files.")
    args = parser.parse_args()
    
    # If the user did not provide a directory path, ask for it
    if not args.directory_path:
        args.directory_path = input("Please provide the path containing PDB files: ")
        
    directory_path = args.directory_path

    # Find all PDB files in directory
    pdb_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.pdb')]
    # Initialize tqdm progress bar
    progress_bar = tqdm(total=len(pdb_files), desc="Processing PDB files")
    
    start_time = time.time()
    for pdb_file in pdb_files:
        if is_pdb_cleaned(pdb_file):
            print(f"{pdb_file} has already been cleaned. Skipping to the next file...")
            progress_bar.update(1)
            continue
            
        information_lines = parse_pdb_information(pdb_file)
        parsed_chains, information_string = chain_information(information_lines, pdb_file)
        pdb_cleaner(pdb_file, parsed_chains, information_string)
        progress_bar.update(1)
    end_time = time.time()
    progress_bar.close()
    
    print(f"Finished processing {len(pdb_files)} files in {end_time - start_time:.2f} seconds") 