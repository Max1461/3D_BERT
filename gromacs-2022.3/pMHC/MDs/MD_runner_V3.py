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
from pathlib import Path
from openmm.app import PDBFile
from pdbfixer import PDBFixer

def extract_information_from_pdb(pdb_file):
    """
    Extracts various information from the given PDB file including:
    - The chain identifiers found in the PDB file, sorted in alphabetical order.
    - The peptide chain ID
    - The G-domain range strings

    Args:
        pdb_file (str): The path to the PDB file.

    Returns:
        A dictionary containing the extracted variables. The keys are 'chains', 'peptide_ID', and 'g_domain_range'
    """

    pdb_ID = Path(pdb_file).resolve().stem.lower()

    # Open the PDB file
    with open(pdb_file, "r") as file:
        pdb_data = file.read()
        lines = pdb_data.splitlines()

    # Use a regular expression to find the lines that contain the chain IDs
    chain_lines = re.findall(r"COMPND\s+\d+\s+CHAIN:\s+([A-Z,\s]+)", pdb_data)
    # Extract the chain IDs from the lines and removes the spaces and ','
    chains = [chain for sublist in chain_lines for chain in sublist.split(", ")]

    peptide_chains = re.findall(r"IMGT\sreceptor\sdescription.+?Peptide.+?Chain\sID.+?{}_(\w)".format(pdb_ID), pdb_data, re.MULTILINE | re.DOTALL)

    chain_lines_dict = {chain_id: [] for chain_id in chains}

    if not peptide_chains:
        print("No peptide chain found in pdb file documentation, finding smallest chain to search for peptide")
        # Initialize a dictionary to store the lines for each chain
        chain_lines = {chain_id: [] for chain_id in chains}
        # Open the PDB file in read mode
        with open(pdb_file, "r") as f:
            # Iterate over each line in the file
            for line in f:
                # Check if the line corresponds to an ATOM or HETATM record for one of the specified chain IDs
                if line.startswith("ATOM  ") or line.startswith("HETATM") and line[21] in chains:
                    # If so, add the line to the list for the appropriate chain
                    chain_lines_dict[line[21]].append(line)
        peptide_chains = min(chain_lines_dict, key=lambda k: len(chain_lines_dict[k]))
        petide_ID = peptide_chains[0]

    if len(peptide_chains)>1:
        print("Found multiple peptides in documentation, selecting first occurence(for now)")
        petide_ID = peptide_chains[0]
    else:
        petide_ID = peptide_chains[0]

    g_domain = ""

    for i, line in enumerate(lines):
            if "G-DOMAIN" in line:
                for j in range(i-11, i):
                    g_domain += lines[j]
                break
    else:
        print("G-DOMAIN not found in file.")

    regex = r"\[\s+?G.*?-.*?A.*?L.*?P.*?H.*?A.*?\d\s\((\d+-\d+)\)\s\[D\d+\]"
    g_domain_range_strings = re.findall(regex, g_domain)

    if len(g_domain_range_strings) != 2:
        g_domain_parts = len(g_domain_range_strings)
        print("Expect to find two parts of the g-domain range, alpha 1 and alpha 2, but found {}".format(g_domain_parts))

    pdb_information = {
    'chains': chains,
    'peptide_ID': petide_ID,
    'g_domain_range': g_domain_range_strings
    }
    return pdb_information

def process_pdb_file(pdb_file, pdb_information):
    """
    Processes the given PDB file. The function uses PDBFixer to remove heterogens, 
    find missing residues and atoms, and add missing heavy atoms.
    The function also uses gmx to convert the PDB file to a GRO file,
    and finds the index numbers for the expected NH3+ start and COO- end terminals
    and uses them to run the final version of gmx pdb2gmx. 

    Args:
        pdb_file (str): The path to the PDB file to be processed.
    Returns:
        processed_gro_file (str) :  The processed gro file.
        gmx_cleaned_pdb_file (str) :  The processed gmx cleaned pdb file.

    """
    chains = pdb_information["chains"]
    chains_count = len(chains)
    pdb_path = pdb_file
    # Get the file name with the extension
    pdb_file = os.path.basename(pdb_file)

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
    gmx_input = " 1 1"
    proc = subprocess.run(["echo{} | gmx pdb2gmx -f {} -o {} -ff charmm36-jul2022 -water tip3p -ter -v -debug".format(gmx_input*chains_count, gmx_cleaned_pdb_file, processed_gro_file)],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

    # Read the initial output from the subprocess
    initial_output = proc.stdout
    initial_output_str = str(initial_output, "utf-8")
    # Find the index number for the expected NH3+ start and COO- end terminals
    ter_index_nrs = [match[0] for match in re.findall(r"\s(\d+): (.*?NH3\+|COO-)", initial_output_str)]
    # Join the index numbers so they can be used as input for the final version to be used
    ter_inputs = " ".join(ter_index_nrs)
    # Run final version
    subprocess.run(["echo {} | gmx pdb2gmx -f {} -o {} -ff charmm36-jul2022 -water tip3p -ter -v -debug".format(ter_inputs, gmx_cleaned_pdb_file, processed_gro_file)],shell=True)
    
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

def create_index(pdb_information):
    """
    Create index.ndx file using the chain identifiers in the given set of chains.

    Args:
        chains: A set of chain identifiers (single characters) in the PDB file.

    Returns:
        None
    """

    chains = pdb_information["chains"]
    peptide_ID = pdb_information["peptide_ID"]
    g_domain_range_strings = pdb_information['g_domain_range']

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
    print(initial_output_str)
    # Get highest current index group number
    counter = max([int(i) for i in re.findall(r"(\d+)\s.+?\s+:\s+\d+\s+atoms", initial_output_str)])

    complete_input = ""
    for chain_ID in chains:

        # Use a regular expression to parse the initial output
        chain_nr = re.findall(r"(\d+).+?ch{}".format(chain_ID), initial_output_str)[0]
        if chain_nr:
            # Generate new index group number
            counter += 1 
            # Based on the parsed output, determine the new input
            new_input = "4 & {}\n".format(chain_nr)
            complete_input += new_input
            if chain_ID == peptide_ID:
                new_input += "name {} Peptide_Backbone\n".format(counter)
                # Rename peptide chain to Peptide for easier downstream working and analysis 
                new_input += "name {} Peptide\n".format(chain_nr)
                complete_input += new_input
            else:
                new_input += "name {} Chain_{}_Backbone\n".format(counter, chain_ID)
                complete_input += new_input

    # Add groove domain to index
    g_domain_groups = []
    g_alpha_count = 1
    for g_domain_range in g_domain_range_strings:
        counter += 1
        g_domain_groups.append(counter)
        new_input = "r {}\nname {} G_Alpha{}\n".format(g_domain_range, counter, g_alpha_count)
        g_alpha_count += 1
        complete_input += new_input

    counter += 1
    new_input = "or {}\nname {} G_domain\n4 & {}\nname {} G_domain_Backbone\n".format(" ".join([str(i) for i in g_domain_groups]), counter, counter, counter+1)
    complete_input += new_input
    complete_input += "q\n"            
    # Send the new input to the subprocess
    subprocess.run(["gmx make_ndx -n index.ndx -o index.ndx -quiet"], input = complete_input, shell=True, text=True)
    
def md_simulation_preparation(pdb_information, n_cpus):
    """
    Runs a Gromacs MD simulation using the specified number of CPUs.

    Args:
        n_cpus: The number of CPUs to use for the simulation.
    """
    cmd = "gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr"
    subprocess.call(cmd, shell=True)

    cmd = f"gmx mdrun -v -deffnm em -pin on -nt {n_cpus}"
    subprocess.call(cmd, shell=True)

    cmd = "echo 0 | gmx trjconv -f em.gro -o em.pdb -s em.tpr"
    subprocess.call(cmd, shell=True)

    # Create index file for gromacs
    # Add the residue range for the G-domain into the index as well
    create_index(pdb_information)

    cmd = "gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr"
    subprocess.call(cmd, shell=True)


    cmd = f"gmx mdrun -v -deffnm nvt -pin on -nt {n_cpus}"
    subprocess.call(cmd, shell=True)

    cmd = "gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr"
    subprocess.call(cmd, shell=True)

    cmd = f"gmx mdrun -deffnm npt -pin on -nt {n_cpus}"
    subprocess.call(cmd, shell=True)
    
def add_energy_groups(pdb_information):
    """
    Add the 'energy-grps' line to the md file, if it does not already exist.

    """
    # Add the groups for energy calcs later
    chains = pdb_information["chains"]
    peptide_ID = pdb_information["peptide_ID"]
    chains = [x for x in chains if x != peptide_ID]
    groups_pattern = " ".join(["ch{}".format(i) for i in chains])
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
                if lines[i+1].strip() != 'energy-grps             = {} Peptide G_domain'.format(groups_pattern):
                    # Add the line after the current line
                    lines.insert(i+1, 'energy-grps             = {} Peptide G_domain\n'.format(groups_pattern))
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

    cmd = f"gmx mdrun -deffnm md_0_100 -pin on -nt {n_cpus}"
    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    pdb_file = argv[1]
    n_cpus = argv[2]

    # Extract the chains from the PDB file
    pdb_information = extract_information_from_pdb(pdb_file)

    # Process the PDB file
    processed_gro_file = process_pdb_file(pdb_file, pdb_information)

    # Process the processed GRO file
    process_gro_file(processed_gro_file)

    # Run temperature and pressure equilibriations
    md_simulation_preparation(pdb_information, n_cpus)

    # Add the energy group lines to the md.mdp file
    add_energy_groups(pdb_information)

    # Run the energy minimization and MD simulations
    run_md_simulation(n_cpus)
