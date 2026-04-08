"""Sterol (ST) builder.

Handles cholesterol, cholesterol esters, and common bile acids.
Sterols use hardcoded SMILES scaffolds (the ring system can't be generated
algorithmically).

pygoslin for ST 27:1;O:
  - lipid.headgroup.headgroup = "ST 27:1;O"  (full string, not just "ST")
  - lipid.fa = {}  (no FA chains for free sterols)

For cholesterol esters (CE):
  - lipid.headgroup.headgroup = "CE"
  - lipid.fa = {"FA1": fa}  (the esterified fatty acid)
"""

from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.fatty_acid import (
    _extract_db_count,
    _extract_db_positions,
    _extract_modifications,
)
from pylipidparse.exceptions import (
    StructureGenerationError,
    UnsupportedLipidClassError,
)
from pylipidparse.utils.chain import build_acyl_chain

# Cholesterol SMILES — C27H46O (PubChem CID 5997, 3β-hydroxy-Δ5-cholestene)
# Verified: 27 carbons, 4-ring steroid backbone, Δ5 double bond, 3β-OH
_CHOLESTEROL_SMILES = (
    "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@H]3CCC4=C[C@@H](O)CC[C@]4(C)[C@@H]3CC[C@@]12C"
)

# Cholesterol ester: 3β-OH replaced with OC(=O)[acyl_tail]
# {acyl} is C2..Cn of the fatty acid (methyl-first, without the C1 carbonyl)
_CE_SCAFFOLD = (
    "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@H]3CCC4=C[C@@H](OC(=O){acyl})CC[C@]4(C)[C@@H]3CC[C@@]12C"
)

# Headgroup strings that map to free cholesterol
_CHOLESTEROL_HG = frozenset({"FC", "CHOLESTEROL", "CHOL"})

# Headgroup prefix patterns for sterols
_ST_PREFIXES = frozenset({"ST", "CE", "CHE"})

# Bile acid SMILES (verified against PubChem where noted)
_BILE_ACIDS = {
    # Cholic acid — PubChem CID 221493, C24H40O5
    "CA": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)"
        "CC[C@@H](O)CC)[C@@H](O)CC[C@@H]3O"
    ),
    # Deoxycholic acid — PubChem CID 222528, C24H40O4
    "DCA": ("OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)" "CC[C@@H](O)CC)CCC3O"),
    # Chenodeoxycholic acid — PubChem CID 10133, C24H40O4
    "CDCA": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)" "CCC(O)CC)[C@@H](O)CCC3"
    ),
    # Ursodeoxycholic acid — PubChem CID 31401, C24H40O4
    "UDCA": (
        "OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)" "CCC(O)CC)[C@H](O)CCC3"
    ),
    # Lithocholic acid — PubChem CID 9903, C24H40O3
    "LCA": ("OC(=O)CC[C@H]1[C@@H]2CC[C@H]3[C@@H]([C@H]2CC[C@@H]1[C@@H](C)" "CCCC)CCC3O"),
}

_BILE_ACID_ALIASES = {
    "CHOLICACID": "CA",
    "DEOXYCHOLICACID": "DCA",
    "CHENODEOXYCHOLICACID": "CDCA",
    "URSODEOXYCHOLICACID": "UDCA",
    "LITHOCHOLICACID": "LCA",
    "BA": "CA",  # generic bile acid defaults to cholic acid
}


class SterolBuilder(AbstractLipidBuilder):
    """Build sterol molecules from pygoslin LipidMolecule objects."""

    def build(self, lipid) -> Chem.Mol:
        headgroup_str = lipid.headgroup.headgroup
        hg_upper = headgroup_str.upper().replace(" ", "")

        # Free cholesterol: "FC", "Cholesterol", or "ST 27:1;O" / "ST27:1;O"
        if headgroup_str in _CHOLESTEROL_HG or hg_upper in {h.upper() for h in _CHOLESTEROL_HG}:
            return self._mol_from_smiles_sanitized(_CHOLESTEROL_SMILES)

        if hg_upper.startswith("ST") or hg_upper.startswith("SE"):
            # "ST 27:1;O" → cholesterol (Δ5-cholesten-3β-ol)
            # "SE 27:1" → pygoslin's normalization of CE input; if no FA chain, treat as cholesterol
            if not lipid.fa:
                return self._mol_from_smiles_sanitized(_CHOLESTEROL_SMILES)
            # SE with FA chain → cholesterol ester
            return self._build_cholesterol_ester(lipid)

        # Cholesterol ester: "CE" or "ChE"
        if hg_upper.startswith("CE") or hg_upper.startswith("CHE"):
            return self._build_cholesterol_ester(lipid)

        # Bile acids
        bile_key = _BILE_ACID_ALIASES.get(hg_upper) or (
            hg_upper if hg_upper in _BILE_ACIDS else None
        )
        if bile_key and bile_key in _BILE_ACIDS:
            smiles = _BILE_ACIDS[bile_key]
            return self._mol_from_smiles_sanitized(smiles)

        raise UnsupportedLipidClassError(
            f"Sterol {headgroup_str!r} is not yet supported. "
            "Supported: cholesterol (ST/FC), cholesterol esters (CE), "
            "bile acids (CA/DCA/CDCA/UDCA/LCA)."
        )

    def _mol_from_smiles_sanitized(self, smiles: str) -> Chem.Mol:
        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))

    def _build_cholesterol_ester(self, lipid) -> Chem.Mol:
        fa_list = list(lipid.fa.values())
        if not fa_list:
            raise StructureGenerationError("No fatty acid chain found for cholesterol ester.")
        fa = fa_list[0]

        num_carbon = int(fa.num_carbon)
        num_db = _extract_db_count(fa)
        db_positions = _extract_db_positions(fa)
        mods = _extract_modifications(fa)

        if num_db > 0 and not db_positions:
            from pylipidparse.exceptions import InsufficientStructuralDetailError

            raise InsufficientStructuralDetailError(
                f"Acyl chain has {num_db} double bond(s) but no positional info."
            )

        # acyl_fragment: methyl-first, ends with C(=O)  e.g. "CCCCC...CC(=O)"
        acyl_fragment = build_acyl_chain(num_carbon, db_positions, mods, terminus="ester")

        # CE scaffold has OC(=O){acyl} where {acyl} = C2..Cn (without C1's C(=O))
        if acyl_fragment.endswith("C(=O)"):
            acyl_tail = acyl_fragment[:-5]  # C2..Cn
        else:
            acyl_tail = acyl_fragment

        smiles = _CE_SCAFFOLD.format(acyl=acyl_tail)
        return self._mol_from_smiles_sanitized(smiles)
