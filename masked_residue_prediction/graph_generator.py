import argparse
import os
from deeprankcore.query import QueryCollection, ProteinPeptideAtomicQuery, ProteinProteinInterfaceAtomicQuery
from deeprankcore.query import Query
from pathlib import Path
import torch
import importlib

def get_all_atoms( # pylint: disable=too-many-locals
    pdb_path: str,
) -> List[Atom]:
    """Gets the contact atoms from pdb2sql and wraps them in python objects."""

    pdb = pdb2sql(pdb_path)
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]

    rows = pdb.get("x,y,z,name,element,altLoc,occ,chainID,resSeq,resName,iCode")

    structure = PDBStructure(f"atoms_{pdb_name}")

    for row in rows:
        (
            x,
            y,
            z,
            atom_name,
            element_name,
            altloc,
            occupancy,
            chain_id,
            residue_number,
            residue_name,
            insertion_code
        ) = row

        _add_atom_data_to_structure(structure,
                                    x, y, z,
                                    atom_name,
                                    altloc, occupancy,
                                    element_name,
                                    chain_id,
                                    residue_number,
                                    residue_name,
                                    insertion_code)
    return structure.get_atoms()

def _load_all_atoms(pdb_path: str,
                    include_hydrogens: bool) -> List[Atom]:

    # get the contact atoms
    if include_hydrogens:

        pdb_name = os.path.basename(pdb_path)
        hydrogen_pdb_file, hydrogen_pdb_path = tempfile.mkstemp(
            prefix="hydrogenated-", suffix=pdb_name
        )
        os.close(hydrogen_pdb_file)

        add_hydrogens(pdb_path, hydrogen_pdb_path)

        try:
            all_atoms = get_all_atoms(hydrogen_pdb_path)
        finally:
            os.remove(hydrogen_pdb_path)
    else:
        all_atoms = get_all_atoms(pdb_path)

    if len(all_atoms) == 0:
        raise ValueError("no contact atoms found")

    return all_atoms

class ProteinPeptideAtomicQuery(Query):

    def __init__(  # pylint: disable=too-many-arguments
        self,
        pdb_path: str,
        targets: Optional[Dict[str, float]] = None,
        distance_cutoff: Optional[float] = 20
    ):
        """
        A query that builds atom-based graphs, using the residues at a protein-protein interface.

        Args:
            pdb_path (str): The path to the .PDB file.
            distance_cutoff (Optional[float], optional): Max distance in Ångström between two interacting atoms of the two proteins.
                Defaults to 5.5.
            targets (Optional[Dict(str,float)], optional): Named target values associated with this query. Defaults to None.
        """

        model_id = os.path.splitext(os.path.basename(pdb_path))[0]

        Query.__init__(self, model_id, targets)

        self._pdb_path = pdb_path
        self._distance_cutoff = distance_cutoff

    def get_query_id(self) -> str:
        "Returns the string representing the complete query ID."
        return f"atom-prot_pep-{self.model_id}"

    def __eq__(self, other) -> bool:
        return (
            isinstance(self, type(other))
            and self.model_id == other.model_id
        )

    def __hash__(self) -> hash:
        return hash((self.model_id))

    def build(self, feature_modules: List[ModuleType], include_hydrogens: bool = False) -> Graph:
        """Builds the graph from the .PDB structure.

        Args:
            feature_modules (List[ModuleType]): Each must implement the :py:func:`add_features` function.
            include_hydrogens (bool, optional): Whether to include hydrogens in the :class:`Graph`. Defaults to False.

        Returns:
            :class:`Graph`: The resulting :class:`Graph` object with all the features and targets. 
        """

        all_atoms = _load_all_atoms(self._pdb_path,
                                        include_hydrogens)

        # build the graph
        graph = build_atomic_graph(
            all_atoms, self.get_query_id(), self._distance_cutoff
        )

        # add data to the graph
        self._set_graph_targets(graph)

        # add the features
        for feature_module in feature_modules:
            feature_module.add_features(self._pdb_path, graph)

        graph.center = np.mean([atom.position for atom in all_atoms], axis=0)
        return graph

def graph_generation(pdb_files):

    residue_encoding = {
    'ala': torch.tensor([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'arg': torch.tensor([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'asn': torch.tensor([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'asp': torch.tensor([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'cys': torch.tensor([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'gln': torch.tensor([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'glu': torch.tensor([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'gly': torch.tensor([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'his': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'ile': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'leu': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    'lys': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    'met': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]),
    'phe': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]),
    'pro': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),
    'ser': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]),
    'thr': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
    'trp': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]),
    'tyr': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),
    'val': torch.tensor([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    }


    queries = QueryCollection()

    for pdb_file in pdb_files:
        pdb_ID = Path(pdb_file).resolve().stem.lower()
        masked_residue = pdb_ID.split('_')[2]
        # Append data points
        queries.add(ProteinProteinInterfaceAtomicQuery(
        pdb_path = pdb_file,
        chain_id1 = "A",
        chain_id2 = "C",
        targets = {
            "residue": residue_encoding[masked_residue]
        }
        ))

        queries.add(ProteinPeptideAtomicQuery(
        pdb_path = pdb_file,
        targets = {
            "residue": residue_encoding[masked_residue]
        }
        ))

    output_directory = "/test/test"
    # Generate graphs and save them in hdf5 files
    feature_names = ['components', 'contact', 'exposure', 'surfacearea']
    feature_modules = [importlib.import_module('deeprankcore.features.' + name) for name in feature_names]
    output_paths = queries.process(feature_modules = feature_modules)

def get_predicted_aa(predicted_tensor, residue_encoding):
    for aa, encoding in residue_encoding.items():
        if torch.all(torch.eq(predicted_tensor, torch.tensor(encoding))):
            return aa

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory_path", help="The path to directory where micro-environment PDB files are located.")
    args = parser.parse_args()

    # If the user did not provide a directory path, ask for it
    if not args.directory_path:
        args.directory_path = input("Please provide the path containing micro-environment PDB files: ")

    # Get the directory path from the user input
    directory_path = args.directory_path

    # Find all PDB files in directory
    pdb_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.pdb')]