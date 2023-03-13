from typing import Optional
import numpy as np
from deeprankcore.utils.graph import Graph
from deeprankcore.utils.parsing import atomic_forcefield
from deeprankcore.molstruct.variant import SingleResidueVariant
from deeprankcore.molstruct.atom import Atom
from deeprankcore.molstruct.residue import Residue
from deeprankcore.domain import nodestorage as Nfeat
import logging
from pathlib import Path

_log = logging.getLogger(__name__)

def add_features( # pylint: disable=unused-argument
    pdb_path: str,
    graph: Graph,
    single_amino_acid_variant: Optional[SingleResidueVariant] = None
    ):

    pdb_ID = Path(pdb_path).resolve().stem.lower()
    residue_number = int(pdb_ID.split('_')[1])

    for node in graph.nodes:
        if isinstance(node.id, Residue):
            residue = node.id
            chain = residue._chain
            number = residue._number
            if chain == "C" and number == residue_number:
                node.features[Nfeat.TARGETNODE] = 1
            else:
                node.features[Nfeat.TARGETNODE] = 0
        elif isinstance(node.id, Atom):
            atom = node.id
            residue = atom.residue
            chain = residue._chain
            number = residue._number
            if chain == "C" and number == residue_number and atom.name == "CA":
                node.features[Nfeat.TARGETNODE] = 1
            else:
                node.features[Nfeat.TARGETNODE] = 0
        else:
            raise TypeError(f"Unexpected node type: {type(node.id)}") 
        