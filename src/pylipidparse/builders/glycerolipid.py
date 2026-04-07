"""Glycerolipid (GL) builder: MAG, DAG, TAG.

Handles mono-, di-, and triacylglycerols with ester and ether linkages.
"""
from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.fatty_acid import _extract_modifications
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    StructureGenerationError,
    UnsupportedLipidClassError,
)
from pylipidparse.scaffolds.headgroups import get_glycerolipid_scaffold
from pylipidparse.utils.chain import build_acyl_chain, build_alkyl_chain

# pygoslin FA bond type enum values (may vary by version)
_ETHER_PLASMALOGEN = "ETHER_PLASMALOGEN"
_ETHER_ALKYL = "ETHER_ALKYL"
_ESTER = "ESTER"
_NO_FA = "NO_FA"


def _get_bond_type(fa) -> str:
    """Get the bond type string from a pygoslin FA object."""
    bt = getattr(fa, "lipid_FA_bond_type", None)
    if bt is None:
        return _ESTER
    bt_str = str(bt).upper()
    if "PLASMALOGEN" in bt_str or bt_str.endswith("_P"):
        return _ETHER_PLASMALOGEN
    if "ALKYL" in bt_str or bt_str.endswith("_O"):
        return _ETHER_ALKYL
    if "NO_FA" in bt_str:
        return _NO_FA
    return _ESTER


def _get_db_positions(fa) -> dict:
    """Extract double bond positions from a pygoslin FA object."""
    db_positions = {}
    num_db = int(getattr(fa, "num_double_bonds", 0))
    raw_db = getattr(fa, "double_bond_positions", None)

    if raw_db:
        if hasattr(raw_db, "get_positions"):
            try:
                for pos, geom in raw_db.get_positions().items():
                    db_positions[int(pos)] = str(geom) if geom else "Z"
            except AttributeError:
                pass
        elif hasattr(raw_db, "items"):
            for pos, geom in raw_db.items():
                db_positions[int(pos)] = str(geom) if geom else "Z"

    if num_db > 0 and not db_positions:
        raise InsufficientStructuralDetailError(
            f"FA chain has {num_db} double bond(s) but no positional information. "
            "Provide full structural notation with double bond positions."
        )
    return db_positions


def _build_chain_fragment(fa, position_label: str = "") -> str:
    """Build the chain SMILES fragment for one FA chain.

    Returns the fragment ending at C1 (appropriate terminus for scaffold
    substitution).
    """
    bond_type = _get_bond_type(fa)
    num_carbon = int(fa.num_carbon)
    db_positions = _get_db_positions(fa)
    mods = _extract_modifications(fa)

    if bond_type == _ETHER_PLASMALOGEN:
        return build_alkyl_chain(num_carbon, db_positions, mods, plasmalogen=True)
    elif bond_type == _ETHER_ALKYL:
        return build_alkyl_chain(num_carbon, db_positions, mods, plasmalogen=False)
    elif bond_type == _ESTER:
        return build_acyl_chain(num_carbon, db_positions, mods, terminus="ester")
    elif bond_type == _NO_FA:
        return ""  # This position is empty (lyso)
    else:
        raise UnsupportedLipidClassError(
            f"Unsupported FA bond type: {bond_type!r} at position {position_label}"
        )


class GlycerolipidBuilder(AbstractLipidBuilder):
    """Build glycerolipid (MAG, DAG, TAG) molecules from pygoslin parsed objects."""

    def build(self, lipid) -> Chem.Mol:
        """Build a glycerolipid molecule.

        Parameters
        ----------
        lipid :
            Parsed pygoslin lipid object.

        Returns
        -------
        Chem.Mol
            Sanitized RDKit molecule.
        """
        headgroup = lipid.headgroup.headgroup.upper()
        fa_dict = lipid.fa  # dict like {"FA1": fa1, "FA2": fa2, ...} or {"sn1": ...}

        # Determine which sn-positions are present and their chains
        sn1_fa, sn2_fa, sn3_fa = _extract_sn_chains(fa_dict, headgroup)

        sn1_present = sn1_fa is not None and _get_bond_type(sn1_fa) != _NO_FA
        sn2_present = sn2_fa is not None and _get_bond_type(sn2_fa) != _NO_FA
        sn3_present = sn3_fa is not None and _get_bond_type(sn3_fa) != _NO_FA

        # Require known sn-positions for structure generation
        # (underscore-separated unknown positions → raise error)
        if not _has_known_positions(lipid):
            raise InsufficientStructuralDetailError(
                f"Cannot generate a unique structure for {lipid.headgroup.headgroup} "
                "with unknown sn-positions (underscore notation). "
                "Use slash notation to specify positions, e.g. 'TG 16:0/18:1/18:2'."
            )

        scaffold = get_glycerolipid_scaffold(
            headgroup, sn1_present, sn2_present, sn3_present
        )

        # Build chain fragments
        subs = {}
        if sn1_present:
            subs["sn1"] = _build_chain_fragment(sn1_fa, "sn1")
        if sn2_present:
            subs["sn2"] = _build_chain_fragment(sn2_fa, "sn2")
        if sn3_present:
            subs["sn3"] = _build_chain_fragment(sn3_fa, "sn3")

        smiles = scaffold.format(**subs)
        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))


def _has_known_positions(lipid) -> bool:
    """Return True if all chain positions are known (slash notation in input)."""
    # pygoslin tracks this via the 'level' or via the FA position info
    # We check if any FA has position = 0 or unset, which indicates unknown
    for fa in lipid.fa.values():
        pos = getattr(fa, "position", -1)
        # If position is -1 or 0 and there are multiple chains, positions are unknown
        if pos is not None and int(pos) == -1:
            return False
    return True


def _extract_sn_chains(fa_dict, headgroup):
    """Map the pygoslin FA dictionary to (sn1, sn2, sn3) FA objects.

    pygoslin uses keys like "FA1", "FA2", "FA3" for glycerolipids.
    For DG/TAG with known positions, the keys may have position metadata.
    """
    values = list(fa_dict.values())
    keys = list(fa_dict.keys())

    # Build a mapping from sn-position to FA object
    sn_map = {1: None, 2: None, 3: None}

    for key, fa in zip(keys, values):
        key_upper = str(key).upper()
        # Try to determine position from key name
        if "1" in key_upper or key_upper in ("FA1", "SN1", "FA_1"):
            sn_map[1] = fa
        elif "2" in key_upper or key_upper in ("FA2", "SN2", "FA_2"):
            sn_map[2] = fa
        elif "3" in key_upper or key_upper in ("FA3", "SN3", "FA_3"):
            sn_map[3] = fa
        else:
            # Fall back: assign sequentially
            for i in range(1, 4):
                if sn_map[i] is None:
                    sn_map[i] = fa
                    break

    return sn_map[1], sn_map[2], sn_map[3]
