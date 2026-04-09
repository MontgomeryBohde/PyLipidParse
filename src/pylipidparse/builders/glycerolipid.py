"""Glycerolipid (GL) builder: MAG, DAG, TAG."""

from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.fatty_acid import (
    _extract_db_count,
    _extract_db_positions,
    _extract_modifications,
)
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
)
from pylipidparse.scaffolds.headgroups import get_glycerolipid_scaffold
from pylipidparse.utils.chain import build_acyl_chain, build_alkyl_chain

# LipidFaBondType enum name fragments (checked via str comparison for version safety)
_BT_PLASMENYL = "PLASMENYL"  # ETHER_PLASMENYL = P- (plasmalogen)
_BT_PLASMANYL = "PLASMANYL"  # ETHER_PLASMANYL = O- (alkyl ether)
_BT_ETHER = "ETHER"  # generic ether
_BT_ESTER = "ESTER"
_BT_NO_FA = "NO_FA"
_BT_LCB = "LCB"  # LCB_REGULAR or LCB_EXCEPTION


def _get_bond_type(fa) -> str:
    """Return a normalised bond-type string from a pygoslin FA object."""
    bt = getattr(fa, "lipid_FA_bond_type", None)
    if bt is None:
        return _BT_ESTER
    bt_str = str(bt).upper()
    if _BT_NO_FA in bt_str:
        return _BT_NO_FA
    if _BT_LCB in bt_str:
        return _BT_LCB
    if _BT_PLASMENYL in bt_str:
        return _BT_PLASMENYL  # P- prefix
    if _BT_PLASMANYL in bt_str:
        return _BT_PLASMANYL  # O- prefix
    if _BT_ETHER in bt_str:
        return _BT_PLASMANYL  # generic ether treated as alkyl ether
    return _BT_ESTER


def _build_chain_fragment(fa, position_label: str = "", c1_first: bool = False) -> str:
    """Build the chain SMILES fragment for one FA chain.

    Parameters
    ----------
    c1_first : bool
        If True, assemble in C1-first (carboxyl-first) order. Required for
        sn-2 and sn-3 positions where the scaffold has ``O{sn}`` so that
        the ester linkage reads ``O-C(=O)-chain`` rather than ``O-chain-C(=O)``.
    """
    bond_type = _get_bond_type(fa)
    num_carbon = int(fa.num_carbon)
    num_db = _extract_db_count(fa)
    db_positions = _extract_db_positions(fa)
    mods = _extract_modifications(fa)

    if num_db > 0 and not db_positions:
        raise InsufficientStructuralDetailError(
            f"FA chain at {position_label} has {num_db} double bond(s) but "
            "no positional info. Provide positions, e.g. '18:1(9Z)'."
        )

    if bond_type == _BT_PLASMENYL:
        # Plasmalogen vinyl-ether: fragment starts with /C=C\ regardless of direction
        return build_alkyl_chain(num_carbon, db_positions, mods, plasmalogen=True)
    elif bond_type in (_BT_PLASMANYL,):
        return build_alkyl_chain(
            num_carbon, db_positions, mods, plasmalogen=False, c1_first=c1_first
        )
    elif bond_type == _BT_NO_FA:
        return ""
    else:  # ESTER or LCB
        return build_acyl_chain(num_carbon, db_positions, mods, terminus="ester", c1_first=c1_first)


class GlycerolipidBuilder(AbstractLipidBuilder):
    """Build glycerolipid molecules from pygoslin LipidMolecule objects."""

    def build(self, lipid) -> Chem.Mol:
        fa_dict = lipid.fa
        headgroup = lipid.headgroup.headgroup.upper()

        # Detect unknown positions (underscore notation → pygoslin may not set positions)
        sn1_fa, sn2_fa, sn3_fa = _extract_sn_chains(fa_dict)

        def _is_present(fa):
            return fa is not None and _get_bond_type(fa) != _BT_NO_FA and int(fa.num_carbon) > 0

        sn1_present = _is_present(sn1_fa)
        sn2_present = _is_present(sn2_fa)
        sn3_present = _is_present(sn3_fa)

        if not any([sn1_present, sn2_present, sn3_present]):
            raise InsufficientStructuralDetailError(
                f"No sn-positioned chains found for {headgroup}. "
                "Use slash notation, e.g. 'TG 16:0/18:1(9Z)/18:2(9Z,12Z)'."
            )

        scaffold = get_glycerolipid_scaffold(headgroup, sn1_present, sn2_present, sn3_present)

        subs = {}
        if sn1_present:
            # sn-1 scaffold pattern: {sn1}O — chain is LEFT of O, methyl-first ends C(=O)O ✓
            subs["sn1"] = _build_chain_fragment(sn1_fa, "sn1", c1_first=False)
        if sn2_present:
            # sn-2 scaffold pattern: O{sn2} — chain is RIGHT of O, needs C1-first so C(=O) leads ✓
            subs["sn2"] = _build_chain_fragment(sn2_fa, "sn2", c1_first=True)
        if sn3_present:
            # sn-3 scaffold pattern: O{sn3} — same as sn-2, needs C1-first ✓
            subs["sn3"] = _build_chain_fragment(sn3_fa, "sn3", c1_first=True)

        smiles = scaffold.format(**subs)
        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))


def _extract_sn_chains(fa_dict):
    """Map pygoslin FA dict to (sn1, sn2, sn3) FA objects.

    pygoslin uses keys "FA1", "FA2", "FA3" for glycerolipids.
    """
    sn_map = {1: None, 2: None, 3: None}

    for key, fa in fa_dict.items():
        key_str = str(key).upper().replace("FA", "").replace("SN", "").strip("_- ")
        try:
            idx = int(key_str)
            if idx in sn_map:
                sn_map[idx] = fa
                continue
        except ValueError:
            pass
        # Fallback sequential assignment
        for i in range(1, 4):
            if sn_map[i] is None:
                sn_map[i] = fa
                break

    return sn_map[1], sn_map[2], sn_map[3]
