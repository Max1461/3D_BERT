# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:51:54 2022

@author: Max
"""
import argparse
import os
import subprocess
import sys
from sys import argv 
import re
import shutil
import string
from openmm.app import PDBFile
from pdbfixer import PDBFixer

def extract_chains_from_pdb(pdb_file):
    """
    Extracts the chain identifiers from the given PDB file and rewrites the file with the chains sorted in alphabetical order.

    Args:
        pdb_file: The path to the PDB file.

    Returns:
        A set containing the chain identifiers found in the PDB file, sorted in alphabetical order.
    """
    # Open the PDB file in read mode
    with open(pdb_file, "r") as f:
        # Create a set to store the chains that are found in the file
        chains = set()

        # Create a temporary list to store the modified lines
        modified_lines = []

        # Loop over the lines in the file
        for line in f:
            # Check if this line is an ATOM or HETATM record
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                # Extract the chain identifier from the record
                chain_id = line[21]
                # Check if the chain ID has been seen before
                if chain_id in chains:
                    # If it has been seen, skip this line
                    continue
                else:
                    # If it has not been seen, add it to the set of chains and append the line to the modified lines list
                    chains.add(chain_id)
                    modified_lines.append(line)
            # If the line is not an atom or hetatom line, add it to the modified lines list
            else:
                modified_lines.append(line)

        # Sort the chains in alphabetical order
        chains = sorted(chains)

        # Initialize a dictionary to store the lines for each chain
        chain_lines = {chain_id: [] for chain_id in chains}

    # Open the PDB file in read mode
    with open(pdb_file, "r") as f:
        # Iterate over each line in the file
        for line in f:
            # Check if the line corresponds to an ATOM or HETATM record for one of the specified chain IDs
            if line.startswith("ATOM  ") or line.startswith("HETATM") and line[21] in chains:
                # If so, add the line to the list for the appropriate chain
                chain_lines[line[21]].append(line)

    # Replace the first occurrence of each chain's lines in the modified_lines list with the lines of the corresponding chain ID
    for chain_id, lines in chain_lines.items():
        if lines:
            modified_lines[modified_lines.index(lines[0])] = lines

    # Initialize a list to store the final lines to be written to the file
    pdb_lines = []
    # Iterate over the modified lines
    for modified_line in modified_lines:
        # If the modified line is a string, add it to the list of final lines
        if isinstance(modified_line, str):
            pdb_lines.append(modified_line)
        # If the modified line is a list, add each element in the list to
        else:
            for line in modified_line:
                pdb_lines.append(line)
                

    # Open the pdb file in write mode
    with open(pdb_file, 'w') as f:
        # Overwrite the file with the modified lines
        f.writelines(pdb_lines)

    petide_ID = min(chain_lines, key=lambda k: len(chain_lines[k]))

    return chains, petide_ID

def process_pdb_file(pdb_file):
    """
    Processes the given PDB file using pdb4amber, gmx, and sed.

    Args:
        pdb_file: The path to the PDB file to be processed.
    """
    pdb_path = pdb_file
    # Get the file name with the extension
    pdb_file = os.path.basename(pdb_file)

    # Clean the PDB file and add missing heavy atoms
    # cleaned_pdb_file = "cleaned_{}".format(pdb_file)
    # run_line = "pdb_delresname -HOH {} > {}".format(pdb_path, cleaned_pdb_file)
    # subprocess.call(run_line, shell=True)

    gmx_cleaned_pdb_file = "gmx_cleaned_{}".format(pdb_file)

    # Create a PDBFixer object
    fixer = PDBFixer(pdb_file)
    fixer.removeHeterogens(False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Set the residue names to a standard naming convention
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Write the modified PDB file to the output file
    PDBFile.writeFile(fixer.topology, fixer.positions, open(gmx_cleaned_pdb_file, 'w'))

    # Use gmx to process the PDB file
    processed_gro_file = "{}_processed.gro".format(pdb_file.split(".")[0])
    proc = subprocess.run(["echo 1 1 1 1 1 1 | gmx pdb2gmx -f {} -o {} -ff charmm36-jul2022 -water tip3p -ter".format(gmx_cleaned_pdb_file, processed_gro_file)],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

    # Read the initial output from the subprocess
    initial_output = proc.stdout
    initial_output_str = str(initial_output, "utf-8")
    # Find the index number for the expected NH3+ start and COO- end terminals
    ter_index_nrs = [match[0] for match in re.findall(r"\s(\d+): (.*?NH3\+|COO-)", initial_output_str)]
    # Join the index numbers so they can be used as input for the final version to be used
    ter_inputs = " ".join(ter_index_nrs)
    # Run final version
    subprocess.run(["echo {} | gmx pdb2gmx -f {} -o {} -ff charmm36-jul2022 -water tip3p -ter".format(ter_inputs, gmx_cleaned_pdb_file, processed_gro_file)],shell=True)
    
    return processed_gro_file

def process_gro_file(processed_gro_file):
    """
    Processes the given processed GRO file by running several GROMACS commands to prepare it for energy minimization.

    Args:
        processed_gro_file: The name of the processed GRO file to be processed.

    """
    run_line = "gmx editconf -f {} -o newbox.gro -bt dodecahedron -d 1.0".format(processed_gro_file)
    os.system(run_line)

    cmd = "gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro"
    subprocess.call(cmd, shell=True)

    cmd = "gmx grompp -f settings.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1"
    subprocess.call(cmd, shell=True)

    cmd = "echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral"
    subprocess.call(cmd, shell=True)

    # cmd = "gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr"
    # subprocess.call(cmd, shell=True)

    # Run temperature and pressure equilibriations
    # md_simulation_preparation(n_cpus=20)


def change_chain(pdb_file, old_chain, new_chain):
    # Open the input file for reading
    with open(pdb_file, "r") as fin:
        # Read the entire file into a list of lines
        lines = fin.readlines()

    # Open the input file for writing
    with open(pdb_file, "w") as fout:
        # Loop over the lines in the list
        for line in lines:
            # Check if this line is a ATOM or HETATM record
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                # Check if the chain identifier matches the old chain
                if line[21] == old_chain:
                    # Replace the old chain identifier with the new chain identifier
                    line = line[:21] + new_chain + line[22:]

            # Write the modified line to the output file
            fout.write(line)

def create_index(chains):
    """
    Create index.ndx file using the chain identifiers in the given set of chains.

    Args:
        chains: A set of chain identifiers (single characters) in the PDB file.

    Returns:
        None
    """

    chain_identifiers = ""
    for chain in chains:
        chain_identifiers += "chain {}\n".format(chain)

    subprocess.run(["gmx make_ndx -f em.pdb -o index.ndx"], input = "{}q\n".format(chain_identifiers), shell=True, text=True)

    # Start the subprocess with the initial input
    proc = subprocess.run(["echo q |gmx make_ndx -n index.ndx -o index.ndx -quiet"],
                            shell=True,
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    # Read the initial output from the subprocess
    initial_output = proc.stdout
    initial_output_str = str(initial_output, "utf-8")

    complete_input = ""
    for chain_ID in chains:
    # for chain_ID in ["A", "B", "P"]:

        # Use a regular expression to parse the initial output
        chain_nr = re.findall(r"(\d+).+?ch{}".format(chain_ID), initial_output_str)[0]
        if chain_nr:
            # Based on the parsed output, determine the new input
            new_input = "4 & {}\n".format(chain_nr)
            complete_input += new_input

    complete_input += "q\n"            
    # Send the new input to the subprocess
    subprocess.run(["gmx make_ndx -n index.ndx -o index.ndx -quiet"], input = complete_input, shell=True, text=True)
    
def md_simulation_preparation(n_cpus):
    """
    Runs a Gromacs MD simulation using the specified number of CPUs.

    Args:
        n_cpus: The number of CPUs to use for the simulation.
    """
    cmd = "gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr"
    subprocess.call(cmd, shell=True)

    cmd = f"gmx mdrun -v -deffnm em -nt {n_cpus}"
    subprocess.call(cmd, shell=True)

    cmd = "echo 0 | gmx trjconv -f em.gro -o em.pdb -s em.tpr"
    subprocess.call(cmd, shell=True)

    # Adjust system pdb file
    change_chain("em.pdb", "C", "P")
    
    chains[chains.index('C')] = "P"
    # Create index file for gromacs
    create_index(chains)

    cmd = "gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr"
    subprocess.call(cmd, shell=True)


    cmd = f"gmx mdrun -v -deffnm nvt -nt {n_cpus}"
    subprocess.call(cmd, shell=True)

    cmd = "gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr"
    subprocess.call(cmd, shell=True)

    cmd = f"gmx mdrun -deffnm npt -nt {n_cpus}"
    subprocess.call(cmd, shell=True)
    
def add_energy_groups():
    """
    Add the 'energy-grps' line to the md file, if it does not already exist.

    """
    # Add the groups for energy calcs later
    # Flag to check if the line was added or not
    line_added = False

    # Open the file in read/write mode
    with open('md.mdp', 'r+') as f:
        # Read the file lines into a list
        lines = f.readlines()

        # Loop over the lines in the list
        for i, line in enumerate(lines):
            # Check if the line starts with "rvdw                    = 1.2,"
            if line.startswith('rvdw                    = 1.2'):
                # Check if the next line is the line to add
                if lines[i+1].strip() != 'energy-grps             = chA chB chP':
                    # Add the line after the current line
                    lines.insert(i+1, 'energy-grps             = chA chB chP\n')
                    # Set the flag to True
                    line_added = True

        # Check if the line was added
        if line_added:
            # Move the file pointer to the beginning of the file
            f.seek(0)
            # Write the modified lines to the file
            f.writelines(lines)
            # Truncate the file to the new size
            f.truncate()

def run_md_simulation(n_cpus):
    """
    Run a molecular dynamics simulation using the specified number of CPU cores.
    
    Args:
        n_cpus (int): The number of CPU cores to use for the simulation.
    """
    cmd = "gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_100.tpr"
    subprocess.call(cmd, shell=True)

    cmd = f"gmx mdrun -deffnm md_0_100 -nt {n_cpus}"
    subprocess.call(cmd, shell=True)


def fix_duplicate_chain_ids(pdb_file_path):
    # Open the PDB file for reading and writing
    with open(pdb_file_path, 'r+') as f:
        # Read in the lines of the file
        lines = f.readlines()

        # Create a mapping from chain identifiers to unique identifiers
        chain_map = {}
        unique_id = 0
        for line in lines:
            if line.startswith('ATOM'):
                # Get the current chain identifier
                chain_id = line[21]
                if chain_id not in chain_map:
                    # If this is a new chain identifier, assign it a unique identifier
                    chain_map[chain_id] = string.ascii_uppercase[unique_id]
                    unique_id += 1

        # Seek to the beginning of the file
        f.seek(0)

        # Write the modified lines to the file
        for line in lines:
            if line.startswith('ATOM'):
                # Replace the chain identifier with the unique identifier
                fixed_line = line[:21] + chain_map[line[21]] + line[22:]
                f.write(fixed_line)
            else:
                f.write(line)

        # Truncate the file to the correct length
        f.truncate()

if __name__ == "__main__":
    pdb_file = argv[1]

    # with open("test.txt", "w") as w:
    #     w.write("test")
    
    #Change duplicate chain IDs
    fix_duplicate_chain_ids(pdb_file)

    # Extract the chains from the PDB file
    chains, petide_ID = extract_chains_from_pdb(pdb_file)
    # if petide_ID != "P":
    #     change_chain(pdb_file, petide_ID, "P")
    #     chains = extract_chains_from_pdb(pdb_file)

    # Process the PDB file
    processed_gro_file = process_pdb_file(pdb_file)

    # Process the processed GRO file
    process_gro_file(processed_gro_file)

    # # Adjust system pdb file
    # change_chain("em.pdb", "C", "P")
    
    # chains[chains.index('C')] = "P"
    # # Create index file for gromacs
    # create_index(chains)

    # Run temperature and pressure equilibriations
    md_simulation_preparation(n_cpus=20)

    # Add the energy group lines to the md.mdp file
    add_energy_groups()

    # Run the energy minimization and MD simulations
    run_md_simulation(n_cpus=20)
