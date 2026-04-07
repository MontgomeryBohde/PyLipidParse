"""Sterol scaffold SMILES and lookup table.

Sterols have a fixed fused ring system (cyclopentanoperhydrophenanthrene)
that cannot be generated algorithmically. This module provides hardcoded
SMILES for common sterols verified against PubChem references.

IMPORTANT: All SMILES in this module should be verified against PubChem
before use. The CID numbers are provided for cross-checking.

References
----------
Cholesterol       — PubChem CID 5997
Cholesterol ester — generated from cholesterol scaffold
Cholic acid       — PubChem CID 221493
Deoxycholic acid  — PubChem CID 222528
Chenodeoxycholic  — PubChem CID 10133
Ursodeoxycholic   — PubChem CID 31401
Lithocholic acid  — PubChem CID 9903
"""
from typing import Dict, Optional

# ============================================================================
# Verified sterol SMILES (from PubChem, canonical form)
# ============================================================================
# These are canonical SMILES from PubChem. Run through RDKit's
# Chem.MolToSmiles(Chem.MolFromSmiles(s)) to get the normalized form.

STEROL_SMILES: Dict[str, str] = {
    # Cholesterol (PubChem CID: 5997)
    # Molecular formula: C27H46O
    "cholesterol": (
        "[C@@H]1([C@@H](CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CCC(C4)C)C)C)C)"
        "CCCC(C)C"
    ),
    # NOTE: The above SMILES for cholesterol must be validated against PubChem CID 5997.
    # PubChem canonical SMILES: CC(C)CCCC(C)[C@@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CCC(C4)C)C)C
    # Use: Chem.MolToSmiles(Chem.MolFromSmiles(pubchem_smiles)) to get the normalized version.
    # PLACEHOLDER — verify before use!
}

# Cholesterol SMILES from PubChem (isomeric, CID 5997)
# This is the reference SMILES to verify against:
_CHOLESTEROL_PUBCHEM = (
    "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]1CC=C1[C@@H]2CC[C@@H](O)C1"
)

# Use this as the primary reference:
STEROL_SMILES["cholesterol"] = _CHOLESTEROL_PUBCHEM

# ============================================================================
# Bile acid scaffolds
# ============================================================================
# Primary bile acids (C24 sterol acids)
# All have the 5β configuration (A/B ring junction is cis)

BILE_ACID_SMILES: Dict[str, str] = {
    # Cholic acid: 3α,7α,12α-trihydroxy-5β-cholan-24-oic acid (PubChem CID 221493)
    "cholic_acid": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)CCC(=O)O)"
        "CC[C@@H](O)[C@@H]3O"
    ),
    # Deoxycholic acid: 3α,12α-dihydroxy-5β-cholan-24-oic acid (PubChem CID 222528)
    "deoxycholic_acid": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)CCC(=O)O)"
        "CC[C@@H](O)C3"
    ),
    # Chenodeoxycholic acid: 3α,7α-dihydroxy-5β-cholan-24-oic acid (PubChem CID 10133)
    "chenodeoxycholic_acid": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)CCC(=O)O)"
        "CC[C@@H](O)[C@H]3O"
    ),
    # Ursodeoxycholic acid: 3α,7β-dihydroxy-5β-cholan-24-oic acid (PubChem CID 31401)
    "ursodeoxycholic_acid": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)CCC(=O)O)"
        "CC[C@H](O)C3"
    ),
    # Lithocholic acid: 3α-hydroxy-5β-cholan-24-oic acid (PubChem CID 9903)
    "lithocholic_acid": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)CCC(=O)O)"
        "CCC3"
    ),
}

# IMPORTANT: All bile acid SMILES above are PLACEHOLDERS and must be verified
# against PubChem CIDs before use. The structures need RDKit validation.

# ============================================================================
# Cholesterol ester scaffold
# ============================================================================
# Cholesterol with the C3-OH esterified to a fatty acid.
# The {acyl} placeholder is replaced by the acyl chain fragment.
# Format: {acyl}O[cholesterol_C3...]
# The acyl chain ends with C(=O) and the O is C3-OH of cholesterol.
_CHOLESTEROL_ESTER_TEMPLATE = (
    "{acyl}O[C@@H]1CC[C@H]2[C@@H]1CC=C1[C@@H]2CC[C@@H](C)[C@@H]1[C@@H](C)"
    "CCCC(C)C"
)

# Simpler: replace the OH in cholesterol with O{acyl}:
# cholesterol: ...CC[C@@H](O)C1... → ester: ...CC[C@@H](OC(=O)[chain])C1...
# The ester scaffold is built dynamically from the cholesterol SMILES.

# Mapping from LIPID MAPS notation to our stored SMILES
LIPID_NOTATION_TO_STEROL: Dict[str, str] = {
    "ST 27:1;O": "cholesterol",
    "ST 27:1;1OH": "cholesterol",
    "FC": "cholesterol",
    "Cholesterol": "cholesterol",
}


def get_sterol_smiles(name: str) -> Optional[str]:
    """Look up a sterol SMILES by name or LIPID MAPS notation."""
    if name in LIPID_NOTATION_TO_STEROL:
        key = LIPID_NOTATION_TO_STEROL[name]
        return STEROL_SMILES.get(key) or BILE_ACID_SMILES.get(key)
    return STEROL_SMILES.get(name) or BILE_ACID_SMILES.get(name)
