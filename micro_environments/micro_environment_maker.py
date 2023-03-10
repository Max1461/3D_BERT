# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 10:48:00 2023

@author: Max
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 09:05:20 2023

@author: Max
"""

import os
from Bio.PDB import PDBParser, PDBIO
from numpy.linalg import norm
from Bio.PDB import Select

from tqdm import tqdm
import time

from multiprocessing import Pool, cpu_count

from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.NeighborSearch import NeighborSearch

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


def get_chain_residues(pdb_file, chain_id):
    """
    Returns the parsed PDB structure and a list of residues that make up the specified chain in the PDB file.

    Args:
        pdb_file (str): The path to the PDB file.
        chain_id (str): The identifier of the chain to select.

    Returns:
        A tuple containing the parsed PDB structure and a list of residues that make up the specified chain.
    """
    # Initialize parser and read PDB file
    parser = PDBParser()
    structure = parser.get_structure('example', pdb_file)

    # Get the alpha carbons of the residues in the specified chain
    chain_residues = []
    for model in structure:
        for chain in model:
            if chain.get_id() == chain_id:
                for residue in chain:
                    chain_residues.append(residue)

    return structure, chain_residues

def get_alpha_carbons(pdb_file, chain_id):
    """
    Returns a list of alpha carbons of the residues in the specified chain of the PDB file.

    Args:
        pdb_file (str): The path to the PDB file.
        chain_id (str): The identifier of the chain to select.

    Returns:
        A tuple containing the parsed PDB structure and a list of alpha carbon atoms.
    """
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
    """
    Returns a boolean indicating whether the given atom is within the specified radius of the center atom.

    Args:
        atom (Atom): The atom to check.
        center_atom (Atom): The center atom.
        radius (float): The radius to check against.

    Returns:
        True if the atom is within the radius, False otherwise.
    """
    distance = norm(atom.get_coord() - center_atom.get_coord())
    if distance < radius:
        return True
    return False

class AtomSelect(Select):
    """
    Select object to use with the PDBIO module.

    Args:
        selected_atoms (list): A list of atoms to select.

    Methods:
        accept_atom(atom): Returns True if the given atom is in the list of selected atoms, False otherwise.
    """
    def __init__(self, selected_atoms):
        self.selected_atoms = selected_atoms

    def accept_atom(self, atom):
        return atom in self.selected_atoms


def remove_side_chain(residue, structure):
    """Remove the side chain atoms from an amino acid residue in a protein structure.

    Args:
        residue (Bio.PDB.Residue): A residue object from the Biopython PDB module.
        structure (Bio.PDB.Structure): A protein structure object from the Biopython PDB module.

    Returns:
        Bio.PDB.Residue: The modified residue object with only the main chain atoms.

    Raises:
        TypeError: If the residue is not an amino acid residue.

    """
    # skip if the residue is not an amino acid and return the original residue
    if not is_aa(residue):
        raise TypeError("The residue is not an amino acid.")
    
    # set of main chain atom names
    main_chain_atoms = {"N", "CA", "C", "O"}
    
    main_chain_coordinates = [residue[a.get_name()].get_coord() for a in residue if a.get_name() in main_chain_atoms]

    # Create a NeighborSearch object with all the atoms in the structure
    ns = NeighborSearch(list(structure.get_atoms()))
    
    # Find all the atoms that are within a distance of 1.5 angstroms from each main chain atom
    for main_chain_coordinate in main_chain_coordinates:
        neighbors = ns.search(main_chain_coordinate, 1.5)
        hydrogen_atoms = [a for a in neighbors if a.element == "H"]
        # Add the names of the hydrogen atoms to the main_chain_atoms set
        for a in hydrogen_atoms:
            main_chain_atoms.add(a.get_name())
    
    # iterate over all child atoms of the residue
    for atom in residue.child_list.copy():
        # detach the atom if it is not part of the main chain
        if atom.name not in main_chain_atoms:
            residue.detach_child(atom.id)
    
    # return the modified residue
    return residue



def select_residues_within_radius(chain_residues, radius, structure, pdb_filename, exclude_same_residue=True, exclude_neighbours=False, exclude_chain=None):
    """
    Selects residues within the specified radius of any atom in a list of central residues in a PDB structure, and writes the selection to a new PDB file.

    Args:
        chain_residues (list): A list of central residues for selecting nearby residues.
        radius (float): The radius to use for selecting residues.
        structure (Structure): The parsed PDB structure.
        pdb_filename (str): The path to the PDB file.
        exclude_same_residue (bool): If True, exclude side chain of residues containing atoms in the same residue as the central residues (default True).
        exclude_neighbours (bool): If True, exclude side chain of residues containing atoms in neighbouring residues of the central residues (default False).
        exclude_chain (str): If specified, exclude side chain of residues containing atoms in the specified chain (default None).
    """
    
    # Extract the name of the PDB file without the extension
    base_filename = os.path.basename(os.path.splitext(pdb_filename)[0])

    # Initialize a counter for generating unique file names
    counter = 1

    # Determine output directory based on exclude options
    if exclude_same_residue:
        output_dir=""
        output_dir = os.path.join(output_dir, "exclude_same_residue")
    if exclude_neighbours:
        output_dir=""
        output_dir = os.path.join(output_dir, "exclude_neighbours")
    if exclude_chain:
        output_dir=""
        output_dir = os.path.join(output_dir, "exclude_peptide_chain")

    # Check if output directory exists and create it if not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Iterate over center residues and select the residues within the radius of each center residue
    for center_residue in chain_residues:
        # Determine the start and end residue indices to exclude neighbours
        if exclude_neighbours:
            chain_id = center_residue.get_parent().get_id()
            residues = center_residue.get_parent().get_residues()
            residue_ids = [residue.get_id()[1] for residue in residues if residue.get_parent().get_id() == chain_id]
            start_index = max(residue_ids.index(center_residue.get_id()[1]) - 2, 0)
            end_index = min(start_index + 5, len(residue_ids))

        selected_residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Check if residue has an atom within radius of center residue atoms
                    if any(is_within_radius(atom, center_atom, radius) for atom in residue.get_list() for center_atom in center_residue.get_list()):
                        # Check for exclude chain
                        if exclude_chain and chain.get_id() == exclude_chain:
                            selected_residues.append(remove_side_chain(residue, structure))
                        # Check for exclude neighbours
                        elif exclude_neighbours and residue.get_id()[1] in residue_ids[start_index:end_index] and residue.get_parent().get_id() == chain_id:
                            selected_residues.append(remove_side_chain(residue, structure))
                        # Check for exclude same residue
                        elif exclude_same_residue and residue.get_id() == center_residue.get_id():
                            selected_residues.append(remove_side_chain(residue, structure))
                        else:
                            selected_residues.append(residue)
        
        # Get atoms from selected residues
        selected_atoms = []
        for selected_residue in selected_residues:
            selected_atoms += selected_residue.get_list()
        # Write selected atoms to a new PDB file
        io = PDBIO()
        io.set_structure(structure)
        output_filename = os.path.join(output_dir, f"{base_filename}_{counter}_{center_residue.get_resname()}.pdb")
        io.save(output_filename, select=AtomSelect(selected_atoms))

        # Increment the counter for the next iteration
        counter += 1

def select_atoms_within_radius(alpha_carbons, radius, structure, pdb_filename, exclude_same_residue=True, exclude_neighbours=False, exclude_chain=None):
    """
    Selects the atoms within the specified radius of the alpha carbons of the specified chain in the PDB file, and writes the selection to a new PDB file.

    Args:
        alpha_carbons (list): A list of alpha carbon atoms.
        radius (float): The radius to use for selecting atoms.
        structure (Structure): The parsed PDB structure.
        pdb_filename (str): The path to the PDB file.
        exclude_same_residue (bool): If True, exclude atoms in the same residue as the alpha carbon (default True).
        exclude_neighbours (bool): If True, exclude atoms in neighbouring residues of the alpha carbon (default False).
        exclude_chain (str): If specified, exclude atoms in the specified chain (default None).
    """
    # Extract the name of the PDB file without the extension
    base_filename = os.path.basename(os.path.splitext(pdb_filename)[0])

    # Initialize a counter for generating unique file names
    counter = 1

    # Determine output directory based on exclude options
    if exclude_same_residue:
        output_dir=""
        output_dir = os.path.join(output_dir, "exclude_same_residue")
    if exclude_neighbours:
        output_dir=""
        output_dir = os.path.join(output_dir, "exclude_neighbours")
    if exclude_chain:
        output_dir=""
        output_dir = os.path.join(output_dir, "exclude_peptide_chain")

    # Check if output directory exists and create it if not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Iterate over alpha carbons and select the atoms within the radius of each alpha carbon
    for center_atom in alpha_carbons:
        # Determine the start and end residue indices to exclude neighbours
        if exclude_neighbours:
            chain_id = center_atom.get_parent().get_parent().get_id()
            residues = center_atom.get_parent().get_parent().get_residues()
            residue_ids = [residue.get_id()[1] for residue in residues if residue.get_parent().get_id() == chain_id]
            start_index = max(residue_ids.index(center_atom.get_parent().get_id()[1]) - 2, 0)
            end_index = min(start_index + 5, len(residue_ids))
        
        selected_atoms = []
        for model in structure:
            for chain in model:
                if exclude_chain and chain.get_id() == exclude_chain:
                    continue
                for residue in chain:
                    if exclude_same_residue and residue.get_id() == center_atom.get_parent().get_id():
                        continue
                    if exclude_neighbours and residue.get_id()[1] in residue_ids[start_index:end_index] and residue.get_parent().get_id() == chain_id:
                        continue
                    for atom in residue:
                        if is_within_radius(atom, center_atom, radius):
                            selected_atoms.append(atom)

        # Write selected atoms to a new PDB file
        io = PDBIO()
        io.set_structure(structure)
        output_filename = os.path.join(output_dir, f"{base_filename}_{counter}.pdb")
        io.save(output_filename, select=AtomSelect(selected_atoms))

        # Increment the counter for the next iteration
        counter += 1

def process_pdb_file(pdb_file):
    """
    Processes a single PDB file by selecting the atoms within the specified radius of the alpha carbons of the peptide chain, and writing the selections to new PDB files.

    Args:
        pdb_file (str): The path to the PDB file.
    """
    # Set radius for new PDB files
    radius = 8.5
    peptide_chain_id = get_peptide_chain_id(pdb_file)
    
    structure, alpha_carbons = get_chain_residues(pdb_file, peptide_chain_id)
    select_residues_within_radius(alpha_carbons, radius, structure, pdb_file)
    select_residues_within_radius(alpha_carbons, radius, structure, pdb_file, exclude_chain=peptide_chain_id)
    select_residues_within_radius(alpha_carbons, radius, structure, pdb_file, exclude_neighbours=True)
    
    # structure, alpha_carbons = get_alpha_carbons(pdb_file, peptide_chain_id)
    # select_atoms_within_radius(alpha_carbons, radius, structure, pdb_file)
    # select_atoms_within_radius(alpha_carbons, radius, structure, pdb_file, exclude_chain=peptide_chain_id)
    # select_atoms_within_radius(alpha_carbons, radius, structure, pdb_file, exclude_neighbours=True)


if __name__ == "__main__":
    # Get directory path from user
    directory_path = input("Enter directory path containing PDB files: ")
    
    # Find all PDB files in directory
    pdb_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.pdb')]
    

    
    # Create a multiprocessing pool with the number of worker processes equal to the number of CPUs
    num_workers = int(0.2 * cpu_count())
    pool = Pool(num_workers)
    
    # Use the pool to process all PDB files in parallel
    progress_bar = tqdm(total=len(pdb_files), desc="Processing PDB files")
    start_time = time.time()
    
    for _ in pool.imap_unordered(process_pdb_file, pdb_files):
        progress_bar.update()
    
    pool.close()
    pool.join()
    
    end_time = time.time()
    progress_bar.close()
    print(f"Finished processing {len(pdb_files)} files in {end_time - start_time:.2f} seconds")
    

