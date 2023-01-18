# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 21:47:58 2022

@author: Max
"""
import os
import argparse
import shutil
import subprocess

def main():
     # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Script for running Gromacs simulations on a directory or pdb file.")
    
    # Add a positional argument for the directory or pdb file path
    parser.add_argument("path", help="Path to a directory containing pdb files or a single pdb file.")

    # Add a new argument for the number of CPUs
    parser.add_argument("--cpu", type=int, help="Number of CPUs to use for the run. Default is 1.", default=1)

    # Parse the command line arguments
    args = parser.parse_args()

    # Get the path and number of CPUs from the parsed arguments
    path = args.path
    cpu = args.cpu
    
    # Get the list of pdb files
    pdb_files = get_pdb_file_list(path)
    
    # If pdb_files is empty, prompt the user for input
    if not pdb_files:
        path = input("Please enter a valid path to a directory or pdb file: ")

    cwd = os.getcwd()
    print(cwd)
    for pdb_file in pdb_files:
        
        
        # Split the file name and extension using os.path.splitext
        file_name, file_ext = os.path.splitext(pdb_file)
        
        # Take the file name from the resulting tuple
        file_name = file_name.split('/')[-1]
        
        output_dir = "{}_MD".format(file_name)
        nwd = os.path.join(cwd, output_dir)
        if not os.path.isdir(nwd):
            os.mkdir(nwd)

        #fetch needed mdp, pdb and MD script files
        get_and_copy_mdp_files("../mdp", nwd)
        copy_file(pdb_file, nwd)
        copy_file("MD_runner_V3.py", nwd)
        MD_script = os.path.join(nwd, "MD_runner_V3.py")
        pdb_file = os.path.join(nwd, os.path.basename(pdb_file))
        os.chdir(os.path.dirname(nwd))
        subprocess.run(["python", "{}".format(MD_script), "{}".format(pdb_file), "{}".format(cpu)], cwd=nwd)

def get_pdb_files(directory):
    # Get all files in the directory
    files = os.listdir(directory)

    # Filter out any files that do not have the .pdb extension
    pdb_files = [f for f in files if f.endswith(".pdb")]

    return pdb_files

def get_pdb_file_list(path):
    # Check if the path is a directory or a file
    if os.path.isdir(path):
        # Get the list of pdb files in the directory
        pdb_files = get_pdb_files(path)
    else:
        # Check if the path is a .pdb file
        if path.endswith(".pdb"):
            # Set the pdb_files variable to a list containing the single .pdb file
            pdb_files = [path]
        else:
            # If the path is not a directory or a .pdb file, set pdb_files to an empty list
            pdb_files = []

    return pdb_files

def get_and_copy_mdp_files(directory, output_directory):
    """
    Get the needed mdp files for Gromacs simulation
    
    Args:
        directory (str): The path of the directory where the files are located.
        output_directory (str): The path of the directory where the files should be copied.
    """    
    # Get all files in the directory
    files = os.listdir(directory)

    # Filter out any files that do not have the .mdp extension
    mdp_files = [f for f in files if f.endswith(".mdp")]

    # Copy the .mdp files to the specified output directory
    for mdp_file in mdp_files:
        shutil.copy(os.path.join(directory, mdp_file), output_directory)

def copy_file(input_file, output_directory):
    """
    Copies the specified input file to the given output directory.
    
    Args:
        input_file (str): The path of the input file to copy.
        output_directory (str): The path of the directory where the file should be copied.
    """
    shutil.copy(input_file, output_directory)

if __name__ == "__main__":
    main()