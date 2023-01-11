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
    # Open the PDB file
    with open(pdb_file, "r") as file:
        pdb_file = file.read()

    # Use a regular expression to find the lines that contain the chain IDs
    chain_lines = re.findall(r"COMPND\s+\d+\s+CHAIN:\s+([A-Z,\s]+)", pdb_file)

    # Extract the chain IDs from the lines and removes the spaces and ','
    chains = [i.replace(" ", "").replace(",", "") for i in chain_lines]

    # Find the index of the lines that contain the word "COMPND"
    compnd_index = [i for i, line in enumerate(pdb_file.split('\n')) if 'COMPND' in line]

    # peptide_chains = []
    # # iterate over the index of the lines that contain the word "COMPND"
    # for i in compnd_index:
    #     lines = pdb_file.split("\n")[i:] # get all lines that follows the line contain "COMPND"
    #     for j,line in enumerate(lines):
    #         if "MOL_ID" in line:
    #             break
    #         if "peptide" in line.lower():
    #             peptide_chain = re.findall(r"CHAIN:\s+([A-Z,\s]+)", line)
    #             if peptide_chain:
    #                 peptide_chains.append(peptide_chain[0].replace(" ", "").replace(",", ""))

    # petide_ID = peptide_chains[0]

    regex = r"G.*-.*A.*L.*P.*H.*A.*?\d\s(\(\d+?-\d+?\))\s\[D\d\]"
    g_domains = re.findall(regex, pdb_file, re.S)

    return chains, g_domains

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


if __name__ == "__main__":
    pdb_file = argv[1]

    # Extract the chains from the PDB file
    chains, petide_ID = extract_chains_from_pdb(pdb_file)

    # Process the PDB file
    processed_gro_file = process_pdb_file(pdb_file)

    # Process the processed GRO file
    process_gro_file(processed_gro_file)

    # Run temperature and pressure equilibriations
    md_simulation_preparation(n_cpus=20)

    # Add the energy group lines to the md.mdp file
    add_energy_groups()

    # Run the energy minimization and MD simulations
    run_md_simulation(n_cpus=20)
