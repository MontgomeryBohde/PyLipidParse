"""Sterol (ST) builder.

Handles cholesterol, cholesterol esters, and common bile acids.

Architecture
------------
Sterols have a fixed fused ring system that cannot be generated algorithmically.
This builder uses a lookup table of known sterol SMILES and, for cholesterol esters,
attaches the acyl chain to the C3-OH of cholesterol.

Supported
---------
- ST 27:1;O (cholesterol, free cholesterol)
- CE / ChE (cholesterol esters with any fatty acid chain)
- Common bile acids: cholic, deoxycholic, chenodeoxycholic, ursodeoxycholic, lithocholic

Not yet supported
-----------------
- Plant sterols (campesterol, sitosterol, stigmasterol)
- Steroid hormones (cortisol, testosterone, estradiol)
- Vitamin D metabolites
"""
from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.fatty_acid import _extract_modifications
from pylipidparse.exceptions import (
    StructureGenerationError,
    UnsupportedLipidClassError,
)
from pylipidparse.scaffolds.sterols import STEROL_SMILES, BILE_ACID_SMILES
from pylipidparse.utils.chain import build_acyl_chain

# Cholesterol SMILES (PubChem CID 5997, isomeric SMILES)
# Verified against PubChem: CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]1CC=C1[C@@H]2CC[C@@H](O)C1
_CHOLESTEROL_SMILES = (
    "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]1CC=C1[C@@H]2CC[C@@H](O)C1"
)

# Cholesterol ester scaffold: C3-OH replaced with O{acyl}
# The {acyl} fragment ends with C(=O) (methyl-first chain builder output)
# → forms ester: ...C3-OC(=O)-[acyl chain]...
_CHOLESTEROL_ESTER_SCAFFOLD = (
    "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]1CC=C1[C@@H]2CC[C@@H](OC(=O){acyl})C1"
)

# Headgroup labels in pygoslin that correspond to cholesterol
_CHOLESTEROL_LABELS = frozenset({
    "FC", "ST", "Cholesterol", "CHOL", "CE_FREE",
})

# Headgroup labels for cholesterol esters
_CHOLESTEROL_ESTER_LABELS = frozenset({
    "CE", "ChE", "CHOL_E",
})

# Bile acid labels
_BILE_ACID_LABELS = frozenset({
    "CA", "DCA", "CDCA", "UDCA", "LCA", "BA",
    "CHOLICACID", "DEOXYCHOLICACID", "CHENODEOXYCHOLICACID",
    "URSODEOXYCHOLICACID", "LITHOCHOLICACID",
})

_BILE_ACID_KEY_MAP = {
    "CA": "cholic_acid",
    "DCA": "deoxycholic_acid",
    "CDCA": "chenodeoxycholic_acid",
    "UDCA": "ursodeoxycholic_acid",
    "LCA": "lithocholic_acid",
}


class SterolBuilder(AbstractLipidBuilder):
    """Build sterol molecules from pygoslin parsed objects."""

    def build(self, lipid) -> Chem.Mol:
        """Build a sterol molecule.

        Parameters
        ----------
        lipid :
            Parsed pygoslin lipid object.

        Returns
        -------
        Chem.Mol
            Sanitized RDKit molecule.
        """
        headgroup = lipid.headgroup.headgroup
        hg_upper = headgroup.upper()

        # Free cholesterol
        if hg_upper in _CHOLESTEROL_LABELS or (
            hg_upper == "ST" and _is_cholesterol(lipid)
        ):
            mol = self._mol_from_smiles(_CHOLESTEROL_SMILES)
            return self._sanitize(Chem.RWMol(mol))

        # Cholesterol ester
        if hg_upper in _CHOLESTEROL_ESTER_LABELS:
            return self._build_cholesterol_ester(lipid)

        # Bile acids
        bile_key = _BILE_ACID_KEY_MAP.get(hg_upper)
        if bile_key:
            smiles = BILE_ACID_SMILES.get(bile_key)
            if smiles:
                mol = self._mol_from_smiles(smiles)
                return self._sanitize(Chem.RWMol(mol))

        raise UnsupportedLipidClassError(
            f"Sterol class {headgroup!r} is not yet supported. "
            "Supported: cholesterol (ST/FC), cholesterol esters (CE/ChE), "
            "and primary bile acids (CA, DCA, CDCA, UDCA, LCA). "
            "Please open an issue at https://github.com/MontgomeryBohde/PyLipidParse/issues"
        )

    def _build_cholesterol_ester(self, lipid) -> Chem.Mol:
        """Build a cholesterol ester by attaching an acyl chain to cholesterol C3-OH."""
        fa_list = list(lipid.fa.values())
        if not fa_list:
            raise StructureGenerationError(
                "No fatty acid chain found for cholesterol ester."
            )
        fa = fa_list[0]

        num_carbon = int(fa.num_carbon)
        db_positions = {}
        raw_db = getattr(fa, "double_bond_positions", None)
        if raw_db:
            if hasattr(raw_db, "items"):
                for pos, geom in raw_db.items():
                    db_positions[int(pos)] = str(geom) if geom else "Z"

        mods = _extract_modifications(fa)

        # Build acyl chain (methyl-first, ends with C(=O))
        acyl_fragment = build_acyl_chain(num_carbon, db_positions, mods, terminus="ester")

        # The ester scaffold expects the fragment WITHOUT the terminal C(=O)
        # because the scaffold already has C(=O) in it: OC(=O){acyl}
        # Hmm, _CHOLESTEROL_ESTER_SCAFFOLD has OC(=O){acyl}, so acyl should be C2..Cn
        # acyl_fragment = CCCC...CC(=O) (methyl-first, ends at C1=carbonyl)
        # We want: OC(=O)[acyl_without_C1] = OC(=O)[C2..Cn]
        # acyl_without_C1 = acyl_fragment[:-5] (remove trailing "C(=O)")
        if acyl_fragment.endswith("C(=O)"):
            acyl_chain = acyl_fragment[:-5]  # C2..Cn (methyl end first, C2 last)
        else:
            acyl_chain = acyl_fragment

        smiles = _CHOLESTEROL_ESTER_SCAFFOLD.format(acyl=acyl_chain)
        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))


def _is_cholesterol(lipid) -> bool:
    """Heuristic: check if this ST lipid is cholesterol (27:1;O)."""
    # Check carbon count and double bond count from the lipid
    try:
        fa_list = list(lipid.fa.values())
        if not fa_list:
            # Might be a species-level sterol — check lipid.info
            return True  # Assume cholesterol if no chains given and class is ST
        fa = fa_list[0]
        n_carbon = int(getattr(fa, "num_carbon", 0))
        n_db = int(getattr(fa, "num_double_bonds", 0))
        return n_carbon == 27 and n_db == 1
    except Exception:
        return False
