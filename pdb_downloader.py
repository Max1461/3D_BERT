# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 14:33:49 2022

@author: Max Luppes
"""
import os
import argparse
import urllib.request
import urllib.parse
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import parse_pdb_header
from Bio import SeqIO
import gzip
import shutil
import PANDORA
from PANDORA import Contacts
from PANDORA import Template
from Bio.PDB import NeighborSearch
from Bio.SeqUtils import seq1
from Bio.PDB import Chain
from string import ascii_uppercase
import re
from PANDORA import Database
import pandas as pd
def download_pdb_files(pdb_ids, download_path):
    """Downloads PDB files for the given PDB IDs.

    Args:
        pdb_ids: A list of PDB IDs.
        download_path: The path where the PDB files will be downloaded.
    """
    # Set the URL for the PDB database
    # pdb_url = "https://www.rcsb.org/pdb/files/"
    

    # Create the download directory if it does not exist
    if not os.path.exists(download_path):
        os.makedirs(download_path)

    # Download the PDB files for each PDB ID
    for pdb_id in pdb_ids:
        pdb_file = pdb_id + ".pdb.gz"
        pdb_url = "https://www.imgt.org/3Dstructure-DB/IMGT-FILE/IMGT-{}".format(pdb_file)
        file_path = os.path.join(download_path, pdb_file)

        # Check if the file already exists in the download path
        if not os.path.exists(file_path) or os.path.exists(file_path.rstrip('.gz')):
            # Download the file if it does not exist
            # urllib.request.urlretrieve(pdb_url + pdb_file, file_path)
            urllib.request.urlretrieve(pdb_url, file_path)
            print(f"PDB file for {pdb_id} successfully downloaded to {download_path}.")
        else:
            # Skip the download if the file already exists
            print(f"PDB file for {pdb_id} already exists in {download_path}. Skipping download.")


def get_pdb_ids(receptor_type):
    """
    Retrieves a list of PDB IDs for the given receptor type.

    Args:
        receptor_type: The type of receptor for which to retrieve PDB IDs.

    Returns:
        A list of PDB IDs.
    """
    # Set the parameters for the IMGT database query
    params = urllib.parse.urlencode({
        'ReceptorType': receptor_type,
        'type-entry': 'PDB'
    })

    # Set the URL for the IMGT database
    url = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi?%s" % params

    # Query the IMGT database
    with urllib.request.urlopen(url) as response:
        text = response.read().decode('utf-8')
        text = text.splitlines()

    # Parse the results from the IMGT database query
    information = [re.sub("\(<a.+?a>\)", "", x) for x in text if re.match("^<t(d|r)",x) if x != "<td></td>"]
    joined_information = ''.join(f"{row}\n" for row in information[1:])
    information_per_ID = re.findall("\<tr\>\n((.+?\n){10})",joined_information)
    data = [re.findall("\">(.*?)<",x[0])[2:] for x in information_per_ID]
    columns = ["IMGT entry ID", "IMGT molecule name", "Species", "IMGT entry type", "IMGT receptor description", "Ligand(s)", "Experimental technique",	"Resolution", "PDB release date"]
    df = pd.DataFrame(data, columns=columns)

    # Filter the results based on the length of the peptide and the HLA allele type
    db = Database.load()
    len_9_PAN_pMHC=([ x.lower() for x in  db.MHCI_data if 'HLA-A*02:01' in db.MHCI_data[x].allele_type and len(db.MHCI_data[x].peptide) == 9])
    df = df[df["IMGT entry ID"].isin(len_9_PAN_pMHC)].sort_values(by=["Resolution"])

    # Return the list of PDB IDs
    return list(df["IMGT entry ID"])

import os
import shutil

def unpack_pdb_gz(directory, remove_zipped = False):
    """
    Unpacks all .pdb.gz files in the specified directory using gzip.open() method

    Args:
        directory (str): The path to the directory where the .pdb.gz files are located.
        remove_zipped (bool, optional): If True remove the original compressed .pdb.gz files after they have been unpacked. Default is False

    Returns:
        None
    """
    for filename in os.listdir(directory):
        if filename.endswith('.pdb.gz'):
            file_path = os.path.join(directory, filename)
            with gzip.open(file_path, 'rb') as f_in:
                with open(file_path.strip('.gz'), 'wb') as f_out:
                    f_out.write(f_in.read())
            if remove_zipped:
                os.remove(file_path)

if __name__ == "__main__":
    # Parse the input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("download_path", help="The path where the PDB files will be downloaded.")
    args = parser.parse_args()

    # If the user did not provide a download path, ask for it
    if not args.download_path:
        args.download_path = input("Please provide the path where the PDB files will be downloaded: ")
    # Read the PDB IDs from IMGT database
    receptor_type='MH1'
    pdb_ids = get_pdb_ids(receptor_type)
    pdb_ids = [i.upper() for i in pdb_ids]

    # Download the PDB files for the PDB IDs
    download_pdb_files(pdb_ids, args.download_path)
    # Unzip PDB files
    unpack_pdb_gz(args.download_path, True)



