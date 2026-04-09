"""Glycerophospholipid (GP) builder."""

from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.glycerolipid import (
    _BT_NO_FA,
    _build_chain_fragment,
    _extract_sn_chains,
    _get_bond_type,
)
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    UnsupportedLipidClassError,
)
from pylipidparse.scaffolds.headgroups import get_gp_scaffold

_FULL_GP_CLASSES = frozenset({"PC", "PE", "PA", "PI", "PS", "PG"})
_LYSO_PREFIX = frozenset({"LPC", "LPE", "LPA", "LPI", "LPS", "LPG"})


class GlycerophospholipidBuilder(AbstractLipidBuilder):
    """Build glycerophospholipid molecules from pygoslin LipidMolecule objects."""

    def build(self, lipid) -> Chem.Mol:
        headgroup = lipid.headgroup.headgroup.upper()
        fa_dict = lipid.fa

        if headgroup in _LYSO_PREFIX:
            base_hg = headgroup[1:]
        else:
            base_hg = headgroup

        if base_hg not in _FULL_GP_CLASSES:
            raise UnsupportedLipidClassError(
                f"Glycerophospholipid class {headgroup!r} is not yet supported. "
                f"Supported: {sorted(_FULL_GP_CLASSES | _LYSO_PREFIX)}"
            )

        sn1_fa, sn2_fa, _ = _extract_sn_chains(fa_dict)

        def _is_present(fa):
            return fa is not None and _get_bond_type(fa) != _BT_NO_FA and int(fa.num_carbon) > 0

        sn1_present = _is_present(sn1_fa)
        sn2_present = _is_present(sn2_fa)

        if not sn1_present and not sn2_present:
            raise InsufficientStructuralDetailError(
                f"No fatty acid chains found for {headgroup}. "
                "Provide full structural notation, e.g. 'PC 16:0/18:1(9Z)'."
            )

        try:
            scaffold = get_gp_scaffold(base_hg, sn1_present, sn2_present)
        except KeyError as exc:
            raise UnsupportedLipidClassError(f"No scaffold available for {headgroup}.") from exc

        subs = {}
        if sn1_present:
            # sn-1 scaffold pattern: {sn1}O — chain LEFT of O, methyl-first ends C(=O)O ✓
            subs["sn1"] = _build_chain_fragment(sn1_fa, "sn1", c1_first=False)
        if sn2_present:
            # sn-2 scaffold pattern: O{sn2} — chain RIGHT of O, needs C1-first ✓
            subs["sn2"] = _build_chain_fragment(sn2_fa, "sn2", c1_first=True)

        smiles = scaffold.format(**subs)
        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))
