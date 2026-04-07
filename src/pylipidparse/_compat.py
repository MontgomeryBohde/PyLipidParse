"""RDKit version compatibility shims.

RDKit's API has changed across versions. This module isolates all
version-specific calls so the rest of the codebase can use a stable interface.
"""
from typing import Optional

from rdkit import Chem


def mol_to_inchi(mol: Chem.Mol) -> Optional[str]:
    """Convert an RDKit Mol to an InChI string.

    Works across RDKit versions (2021.03+).
    """
    try:
        from rdkit.Chem.inchi import MolToInchi  # type: ignore[import]

        return MolToInchi(mol)
    except ImportError:
        pass
    try:
        from rdkit.Chem import MolToInchi  # type: ignore[import]

        return MolToInchi(mol)
    except ImportError:
        pass
    # Fallback for very old RDKit
    try:
        from rdkit.Chem.inchi import MolToInchi as _MolToInchi  # type: ignore[import]

        return _MolToInchi(mol)
    except Exception:
        raise ImportError(
            "Could not import MolToInchi from RDKit. "
            "Ensure your RDKit installation includes InChI support."
        )


def mol_to_inchikey(mol: Chem.Mol) -> Optional[str]:
    """Convert an RDKit Mol to an InChIKey string.

    Works across RDKit versions (2021.03+).
    """
    inchi = mol_to_inchi(mol)
    if inchi is None:
        return None

    try:
        from rdkit.Chem.inchi import InchiToInchiKey  # type: ignore[import]

        return InchiToInchiKey(inchi)
    except ImportError:
        pass
    try:
        from rdkit.Chem import InchiToInchiKey  # type: ignore[import]

        return InchiToInchiKey(inchi)
    except ImportError:
        raise ImportError(
            "Could not import InchiToInchiKey from RDKit. "
            "Ensure your RDKit installation includes InChI support."
        )


def compute_2d_coords(mol: Chem.RWMol) -> None:
    """Add 2D coordinates to a molecule in-place."""
    from rdkit.Chem import AllChem

    AllChem.Compute2DCoords(mol)
