"""Fatty acid (FA) lipid builder.

Handles free fatty acids and fatty acid derivatives parsed from LIPID MAPS
shorthand notation via pygoslin.

Supported classes
-----------------
- FA: Fatty acids (saturated, mono/polyunsaturated, branched, hydroxylated)
- Partial support for FAL (fatty aldehydes), FOH (fatty alcohols)

Notes on pygoslin integration
-------------------------------
The pygoslin ``LipidAdduct`` for a fatty acid has exactly one FA chain.
The chain is accessed as the first (and only) value in ``lipid.fa``.

The ``double_bond_positions`` attribute in pygoslin is a dict-like object
mapping integer position to geometry string ('Z', 'E', or '').
"""
from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    StructureGenerationError,
)
from pylipidparse.utils.chain import build_acyl_chain


class FattyAcidBuilder(AbstractLipidBuilder):
    """Build fatty acid molecules from pygoslin parsed objects."""

    def build(self, lipid) -> Chem.Mol:
        """Build a fatty acid molecule.

        Parameters
        ----------
        lipid :
            Parsed pygoslin lipid object for a fatty acid.

        Returns
        -------
        Chem.Mol
            Sanitized RDKit molecule.

        Examples
        --------
        Palmitic acid (16:0), oleic acid (18:1(9Z)), arachidonic acid (20:4), etc.
        """
        fa_list = list(lipid.fa.values())
        if not fa_list:
            raise StructureGenerationError(
                f"No fatty acid chains found in parsed lipid: {lipid}"
            )
        fa = fa_list[0]

        # Extract chain parameters from pygoslin object
        num_carbon = int(fa.num_carbon)
        num_db = int(fa.num_double_bonds)

        # Double bond positions: pygoslin uses {int: str} mapping
        # If no position info but there are double bonds, we can't place them
        db_positions = {}
        raw_db = getattr(fa, "double_bond_positions", None)
        if raw_db:
            # pygoslin may return a DoubleBonds object or a plain dict
            if hasattr(raw_db, "get_num_db"):
                # It's a DoubleBonds wrapper — iterate its positions
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
                f"Lipid has {num_db} double bond(s) but no positional information. "
                "Provide full structural notation with double bond positions, "
                "e.g. 'FA 18:1(9Z)' instead of 'FA 18:1'."
            )

        # Extract modifications
        modifications = _extract_modifications(fa)

        # Determine terminus
        headgroup = lipid.headgroup.headgroup
        if headgroup in ("FAL",):
            terminus = "aldehyde"
        elif headgroup in ("FOH",):
            # Fatty alcohol: C1 is just CH2OH → but we build as alkyl with OH at C1
            # Treat as free acid chain with C1 = C, then add OH mod at C1
            terminus = "alkyl"
            modifications[1] = "OH"
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


def _extract_modifications(fa) -> dict:
    """Extract modifications (OH, oxo, Me, etc.) from a pygoslin FA object.

    Returns a dict ``{position: mod_type}`` for use with ``build_acyl_chain``.
    """
    mods = {}
    raw_mods = getattr(fa, "modifications", None)
    if not raw_mods:
        return mods

    # pygoslin modifications can be a list, dict, or FattyAcidModifications object
    if hasattr(raw_mods, "items"):
        items = raw_mods.items()
    elif hasattr(raw_mods, "__iter__"):
        items = raw_mods
    else:
        return mods

    for item in items:
        # Item may be (position, mod_name) tuple or a Modification object
        if isinstance(item, (tuple, list)) and len(item) == 2:
            pos, mod_name = item
        elif hasattr(item, "position") and hasattr(item, "modification_name"):
            pos = item.position
            mod_name = item.modification_name
        else:
            continue

        pos = int(pos)
        mod_str = str(mod_name).upper() if mod_name else ""

        if "OH" in mod_str or mod_str == "OH":
            mods[pos] = "OH"
        elif "OXO" in mod_str or "KET" in mod_str:
            mods[pos] = "oxo"
        elif mod_str in ("ME", "METHYL"):
            mods[pos] = "Me"

    return mods
