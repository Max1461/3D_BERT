# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 09:05:20 2023

@author: Max
"""

import os
import math
from Bio.PDB import PDBParser, PDBIO
import numpy as np
from numpy.linalg import norm
from Bio.PDB import Select

# from Bio import PDB
# from Bio.PDB import PDBParser, PDBIO, Selection, Atom, Residue, Chain, Structure, NeighborSearch
#from Bio.PDB.Structure import Structure
#from Bio.PDB.Chain import Chain

#from Bio.PDB.Structure import Structure



def get_peptide_chain_id(pdb_file):
    """
    Determines the ID of the peptide chain in a PDB file by finding the chain with the shortest length.

    Args:
        pdb_file (str): The path to the PDB file.

    Returns:
        The ID of the peptide chain as a string.
    """
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    chains = structure.get_chains()
    peptide_chain_id = None
    min_length = float('inf')
    for chain in chains:
        chain_length = len(list(chain.get_residues()))
        if chain_length < min_length:
            peptide_chain_id = chain.get_id()
            min_length = chain_length
    return peptide_chain_id


def get_alpha_carbons(pdb_file, chain_id):
    # Initialize parser and read PDB file
    parser = PDBParser()
    structure = parser.get_structure('example', pdb_file)

    # Get the alpha carbons of the residues in the specified chain
    alpha_carbons = []
    for model in structure:
        for chain in model:
            if chain.get_id() == chain_id:
                for residue in chain:
                    if residue.has_id('CA'):
                        alpha_carbons.append(residue['CA'])

    return structure, alpha_carbons

def is_within_radius(atom, center_atom, radius):
    distance = norm(atom.get_coord() - center_atom.get_coord())
    if distance < radius:
        return True
    return False

class AtomSelect(Select):
    def __init__(self, selected_atoms):
        self.selected_atoms = selected_atoms

    def accept_atom(self, atom):
        return atom in self.selected_atoms

def select_atoms_within_radius(alpha_carbons, radius, structure, exclude_same_residue=True, exclude_neighbours=False, exclude_chain=None):
    for center_atom in alpha_carbons:
        # Determine the start and end residue indices to exclude neighbours
        if exclude_neighbours:
            residue_ids = [residue.get_id()[1] for residue in center_atom.parent.get_list()]
            start_index = max(residue_ids.index(center_atom.get_parent().get_id()[1]) - 2, 0)
            end_index = min(start_index + 5, len(residue_ids))

        # Iterate over alpha carbons and select the atoms within the radius of each alpha carbon
        selected_atoms = []
        for model in structure:
            for chain in model:
                if exclude_chain and chain.get_id() == exclude_chain:
                    continue
                for residue in chain:
                    if exclude_same_residue and residue.get_id() == center_atom.get_parent().get_id():
                        continue
                    if exclude_neighbours and residue.get_id()[1] in residue_ids[start_index:end_index]:
                        continue
                    for atom in residue:
                        if is_within_radius(atom, center_atom, radius):
                            selected_atoms.append(atom)
        # Write selected atoms to a new PDB file
        io = PDBIO()
        io.set_structure(structure)
        io.save(f'{center_atom.get_id()}.pdb', select=AtomSelect(selected_atoms))




if __name__ == "__main__":
    # Get directory path from user
    directory_path = input("Enter directory path containing PDB files: ")

    # Find all PDB files in directory
    pdb_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.pdb')]

    # Set radius for new PDB files
    radius = 8.5

    # Iterate over PDB files and generate new PDB files for each residue in peptide chain
    for pdb_file in pdb_files:
        peptide_chain_id = get_peptide_chain_id(pdb_file)
        
        structure, alpha_carbons = get_alpha_carbons(pdb_file, peptide_chain_id)
        select_atoms_within_radius(alpha_carbons, radius, structure)



