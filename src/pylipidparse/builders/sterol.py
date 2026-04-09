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
# SMILES taken directly from PubChem CID 5997 canonical isomeric SMILES.
# 4-ring ABCD steroid backbone: A(6)-B(6)-C(6)-D(5), Δ5 double bond, 3β-OH.
_CHOLESTEROL_SMILES = (
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
)

# Cholesterol ester scaffold: 3β-OH replaced with OC(=O){acyl}.
# {acyl} is C2..Cn of the fatty acid (methyl-first, without the C1 carbonyl).
# Ring system is identical to _CHOLESTEROL_SMILES (PubChem CID 5997).
_CE_SCAFFOLD = (
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)OC(=O){acyl})C)C"
)

# Headgroup strings that map to free cholesterol
_CHOLESTEROL_HG = frozenset({"FC", "CHOLESTEROL", "CHOL"})

# Headgroup prefix patterns for sterols
_ST_PREFIXES = frozenset({"ST", "CE", "CHE"})

# Bile acid SMILES — taken from PubChem isomeric SMILES and verified by formula
_BILE_ACIDS = {
    # Cholic acid — PubChem CID 221493, C24H40O5, 3α,7α,12α-trihydroxy-5β-cholanic acid
    "CA": (
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2"
        "[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C"
    ),
    # Deoxycholic acid — PubChem CID 222528, C24H40O4, 3α,12α-dihydroxy-5β-cholanic acid
    "DCA": (
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2"
        "CC[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C"
    ),
    # Chenodeoxycholic acid — PubChem CID 10133, C24H40O4, 3α,7α-dihydroxy-5β-cholanic acid
    "CDCA": (
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2"
        "[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C"
    ),
    # Ursodeoxycholic acid — PubChem CID 31401, C24H40O4, 3α,7β-dihydroxy-5β-cholanic acid
    "UDCA": (
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2"
        "[C@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C"
    ),
    # Lithocholic acid — PubChem CID 9903, C24H40O3, 3α-hydroxy-5β-cholanic acid
    "LCA": (
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2"
        "CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C"
    ),
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

        # Build in C1-first order so the chain reads C2→Cn (left to right after the
        # scaffold's OC(=O)).  Methyl-first order would reverse the directional bonds
        # for unsaturated chains, placing the double bond at the wrong Δ position.
        acyl_c1 = build_acyl_chain(num_carbon, db_positions, mods, terminus="ester", c1_first=True)

        # CE scaffold has OC(=O){acyl} where {acyl} = C2..Cn (strip the leading C(=O))
        if acyl_c1.startswith("C(=O)"):
            acyl_tail = acyl_c1[5:]  # C2..Cn
        else:
            acyl_tail = acyl_c1

        smiles = _CE_SCAFFOLD.format(acyl=acyl_tail)
        return self._mol_from_smiles_sanitized(smiles)
