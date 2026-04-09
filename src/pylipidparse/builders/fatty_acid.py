"""Fatty acid (FA) lipid builder.

pygoslin API (confirmed against pygoslin 2.x):
  - lipid_mol.headgroup.headgroup  → str, e.g. "FA"
  - lipid_mol.fa                   → dict, keys like "FA1"
  - fa.num_carbon                  → int
  - fa.double_bonds                → int OR dict {pos: geom_str}
  - fa.db_num()                    → int (always the count)
  - fa.lipid_FA_bond_type          → LipidFaBondType enum
  - fa.functional_groups           → dict {name: [FunctionalGroup...]}
"""

from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    StructureGenerationError,
)
from pylipidparse.utils.chain import build_acyl_chain


def _extract_db_positions(fa) -> dict:
    """Extract double bond positions from a pygoslin FA object.

    ``fa.double_bonds`` is either:
    - an ``int`` (total count, no positions) → return {}
    - a ``dict`` {position: geometry_str}    → return as-is with int keys

    Returns {} when positions are unknown (triggers InsufficientStructuralDetail
    upstream if there are any double bonds).
    """
    raw = fa.double_bonds
    if isinstance(raw, dict):
        return {int(pos): str(geom) if geom else "Z" for pos, geom in raw.items()}
    return {}


def _extract_db_count(fa) -> int:
    """Return the total double bond count from a pygoslin FA object."""
    if callable(fa.db_num):
        return int(fa.db_num())
    raw = fa.double_bonds
    if isinstance(raw, int):
        return raw
    if isinstance(raw, dict):
        return len(raw)
    return 0


def _extract_modifications(fa) -> dict:
    """Extract modifications (OH, oxo, Me) from a pygoslin FA object.

    Returns ``{position: mod_type}`` for use with ``build_acyl_chain``.
    ``fa.functional_groups`` is a dict like ``{'OH': [FG...], 'oxo': [FG...]}``
    where each FG has ``.position`` and ``.count``.
    """
    mods = {}
    raw = getattr(fa, "functional_groups", None)
    if not raw:
        return mods

    for name, fg_list in raw.items():
        name_upper = str(name).upper()
        for fg in fg_list:
            pos = getattr(fg, "position", -1)
            if pos is None or int(pos) <= 0:
                continue  # unlocalized modification — skip
            pos = int(pos)
            if name_upper in ("OH", "O"):
                mods[pos] = "OH"
            elif name_upper in ("OXO", "KET", "KETO"):
                mods[pos] = "oxo"
            elif name_upper in ("ME", "METHYL"):
                mods[pos] = "Me"

    return mods


class FattyAcidBuilder(AbstractLipidBuilder):
    """Build fatty acid molecules from pygoslin LipidMolecule objects."""

    def build(self, lipid) -> Chem.Mol:
        """Build a fatty acid molecule."""
        fa_list = list(lipid.fa.values())
        if not fa_list:
            raise StructureGenerationError(f"No fatty acid chains found in parsed lipid: {lipid}")
        fa = fa_list[0]

        num_carbon = int(fa.num_carbon)
        num_db = _extract_db_count(fa)
        db_positions = _extract_db_positions(fa)

        if num_db > 0 and not db_positions:
            raise InsufficientStructuralDetailError(
                f"Fatty acid has {num_db} double bond(s) but no positional info. "
                "Provide positions, e.g. 'FA 18:1(9Z)' not 'FA 18:1'."
            )

        modifications = _extract_modifications(fa)

        headgroup = lipid.headgroup.headgroup
        if headgroup in ("FAL",):
            terminus = "aldehyde"
        elif headgroup in ("FOH",):
            terminus = "alcohol"
        else:
            terminus = "free_acid"

        smiles = build_acyl_chain(
            num_carbon=num_carbon,
            double_bond_positions=db_positions,
            modifications=modifications,
            terminus=terminus,
        )

        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))
