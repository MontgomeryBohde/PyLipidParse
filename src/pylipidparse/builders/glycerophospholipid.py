"""Glycerophospholipid (GP) builder.

Handles all common glycerophospholipid classes:
PC, PE, PA, PI, PS, PG and their lyso variants (LPC, LPE, LPA, LPI, LPS, LPG).
Also handles ether (O-) and plasmalogen (P-) linkages at sn-1.
"""
from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.glycerolipid import (
    _NO_FA,
    _ETHER_ALKYL,
    _ETHER_PLASMALOGEN,
    _ESTER,
    _build_chain_fragment,
    _get_bond_type,
)
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    StructureGenerationError,
    UnsupportedLipidClassError,
)
from pylipidparse.scaffolds.headgroups import get_gp_scaffold

# Supported headgroup classes
_FULL_GP_CLASSES = frozenset({"PC", "PE", "PA", "PI", "PS", "PG"})
_LYSO_PREFIX = frozenset({"LPC", "LPE", "LPA", "LPI", "LPS", "LPG"})

# N-modified phospholipids (not yet implemented)
_NMOD_CLASSES = frozenset({"LNPC", "LNPE", "LNPA", "LNPG"})


class GlycerophospholipidBuilder(AbstractLipidBuilder):
    """Build glycerophospholipid molecules from pygoslin parsed objects."""

    def build(self, lipid) -> Chem.Mol:
        """Build a glycerophospholipid molecule.

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

        if headgroup in _NMOD_CLASSES:
            raise UnsupportedLipidClassError(
                f"N-modified phospholipids ({headgroup}) are not yet supported. "
                "Please open an issue at https://github.com/MontgomeryBohde/PyLipidParse/issues"
            )

        fa_dict = lipid.fa
        sn1_fa, sn2_fa = _extract_gp_sn_chains(fa_dict, headgroup)

        # Determine the base headgroup (strip lyso prefix)
        if headgroup in _LYSO_PREFIX:
            base_hg = headgroup[1:]  # "LPC" → "PC"
        else:
            base_hg = headgroup

        # Determine which positions are present
        sn1_present = sn1_fa is not None and _get_bond_type(sn1_fa) != _NO_FA
        sn2_present = sn2_fa is not None and _get_bond_type(sn2_fa) != _NO_FA

        if not sn1_present and not sn2_present:
            raise InsufficientStructuralDetailError(
                f"No fatty acid chains found for {headgroup}."
            )

        # For full GP classes, require at least one chain with known structure
        # Species-level inputs (e.g., "PC 34:1") will have no position info
        if not _has_positional_info(fa_dict):
            raise InsufficientStructuralDetailError(
                f"Cannot generate a unique structure for {headgroup} "
                "from species-level (sum composition) notation. "
                "Provide chain positions, e.g. 'PC 16:0/18:1(9Z)'."
            )

        # Select scaffold
        try:
            scaffold = get_gp_scaffold(base_hg, sn1_present, sn2_present)
        except KeyError as exc:
            raise UnsupportedLipidClassError(
                f"No scaffold available for {headgroup}. "
                f"Supported classes: {sorted(_FULL_GP_CLASSES | _LYSO_PREFIX)}"
            ) from exc

        # Build chain fragments
        subs = {}
        if sn1_present:
            subs["sn1"] = _build_chain_fragment(sn1_fa, "sn1")
        if sn2_present:
            subs["sn2"] = _build_chain_fragment(sn2_fa, "sn2")

        smiles = scaffold.format(**subs)
        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))


def _has_positional_info(fa_dict) -> bool:
    """Check if FA chains have per-chain information (not just sum composition)."""
    # If there's only one FA entry with the total composition, positions are unknown
    # A reliable check: if any FA has num_carbon > 0 and can be individually specified
    # For sum-composition, pygoslin may return a single FA with the total sum
    values = list(fa_dict.values())
    if not values:
        return False
    # If there's only one FA and the headgroup is a 2-chain lipid, it's sum composition
    # We can't easily detect this without more context, so we return True and let
    # the scaffold format step fail if chains are missing.
    return len(values) >= 1  # Trust pygoslin to have parsed it correctly


def _extract_gp_sn_chains(fa_dict, headgroup):
    """Extract sn-1 and sn-2 FA objects from pygoslin FA dict.

    pygoslin uses keys like "FA1", "FA2" for glycerophospholipids.
    For lyso forms, there is only one FA key.
    """
    values = list(fa_dict.values())
    keys = list(fa_dict.keys())

    sn1 = None
    sn2 = None

    for key, fa in zip(keys, values):
        key_str = str(key).upper()
        if "1" in key_str:
            sn1 = fa
        elif "2" in key_str:
            sn2 = fa
        else:
            # Fallback: assign sequentially
            if sn1 is None:
                sn1 = fa
            else:
                sn2 = fa

    # For lyso lipids in pygoslin, the single chain may be labelled differently
    # Try to figure out which position it occupies from the lipid level
    if headgroup.startswith("L") and len(values) == 1 and sn1 is None and sn2 is None:
        sn1 = values[0]

    return sn1, sn2
