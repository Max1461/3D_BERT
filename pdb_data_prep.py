# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:20:22 2023

@author: Max
"""

import argparse
import subprocess

def main():
    """
    Main function for the PDB downloader and cleaner.
    This script downloads PDB files from the IMGT PDB database and cleans them by removing heterogens, finding missing residues and atoms, and adding missing heavy atoms using PDBFixer.

    Args:
        None

    Returns:
        None
    """
    # Parse the input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory_path", help="The path to directory where pdb files should be placed.")
    args = parser.parse_args()

    # If the user did not provide a directory path, ask for it
    if not args.directory_path:
        args.directory_path = input("Please provide the path containing PDB files: ")

    # Get the directory path from the user input
    directory_path = args.directory_path

    # Call the pdb_downloader.py script to download PDB files from the IMGT PDB database
    subprocess.run(["python", "pdb_downloader.py", "{}".format(directory_path)])

    # Call the pdb_cleaner.py script to clean the downloaded PDB files
    subprocess.run(["python", "pdb_cleaner.py", "{}".format(directory_path)])

if __name__ == "__main__":
    main()
