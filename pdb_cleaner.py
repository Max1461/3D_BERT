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