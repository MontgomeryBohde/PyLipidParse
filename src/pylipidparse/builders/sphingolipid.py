"""Sphingolipid (SP) builder.

pygoslin API for ceramides (confirmed):
  - lipid.fa keys: 'LCB' (sphingoid base) and 'FA1' (N-acyl chain)
  - LCB.num_carbon             → int
  - LCB.double_bonds           → int (count only — no positions from pygoslin for LCB!)
  - LCB.db_num()               → int
  - LCB.functional_groups      → {'O': [FG(count=2, position=-1)]} for d-type (;O2)
  - LCB.lipid_FA_bond_type     → LipidFaBondType.LCB_EXCEPTION or LCB_REGULAR
  - FA1.num_carbon             → int
  - FA1.double_bonds           → int or dict
  - FA1.lipid_FA_bond_type     → ESTER

Since pygoslin does NOT give double bond positions for the LCB, we apply
biological defaults:
  - 1 double bond → position 4, geometry E (standard d-erythro sphingosine)
  - 0 double bonds → no double bond (dihydrosphingosine)
"""

from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.fatty_acid import (
    _extract_db_count,
    _extract_db_positions,
    _extract_modifications,
)
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    StructureGenerationError,
    UnsupportedLipidClassError,
)
from pylipidparse.scaffolds.headgroups import SPHINGOLIPID_HEADGROUPS
from pylipidparse.utils.chain import build_acyl_chain


def _get_lcb_n_oh(lcb_fa) -> int:
    """Return number of hydroxyl groups on the LCB.

    pygoslin encodes this via functional_groups: {'O': [FG(count=N)]}.
    """
    raw = getattr(lcb_fa, "functional_groups", {}) or {}
    for name, fg_list in raw.items():
        if str(name).upper() in ("O", "OH"):
            total = sum(int(getattr(fg, "count", 1)) for fg in fg_list)
            if total > 0:
                return total

    # Fallback: check num_oxygens if available
    try:
        n = lcb_fa.num_oxygens()
        if n and int(n) > 0:
            return int(n)
    except (TypeError, AttributeError):
        pass

    return 2  # default: d-type (2 OH groups)


def _build_sphingoid_base_smiles(n_carbon: int, n_oh: int, n_db: int) -> str:
    """Build the sphingoid base SMILES in methyl-first format.

    Returns a template string where:
    - {N_ACYL}  = what goes on C2 nitrogen (e.g. "NC(=O)CCCCC..." or just "N")
    - {C1_HEAD} = the C1 atom with its headgroup substituent (e.g. "CO" for free OH)

    Biological defaults used:
    - d-type (n_oh=2): C1-OH, C3-OH, C4=C5 (4E trans) for n_db=1
    - t-type (n_oh=3): C1-OH, C3-OH, C4-OH, no double bond
    - m-type (n_oh=1): C1-OH only
    """
    if n_oh == 3:
        # Phytosphingosine: C3-OH, C4-OH, no canonical double bond
        # Natural configuration: (2S,3S,4R) — PubChem CID 122121
        # Structure (methyl-first): Cn...C5-C4(OH)-C3(OH)-C2(N)-C1(OH/head)
        tail_n = n_carbon - 4
        tail = "C" * tail_n
        return f"{tail}[C@H](O)[C@@H](O)[C@H]({{N_ACYL}}){{C1_HEAD}}"

    elif n_oh == 2:
        # Standard d-type sphingosine or dihydrosphingosine
        if n_db == 0:
            # Dihydrosphingosine: no double bond
            # Natural configuration: (2S,3R).
            # In methyl-first SMILES both C3 and C2 are [C@@H] to give (2S,3R).
            # Verified against PubChem ceramide (CID 5283564) pattern.
            # Cn...C4-C3(OH)-C2(N)-C1(OH/head)
            tail_n = n_carbon - 3
            tail = "C" * tail_n
            return f"{tail}[C@@H](O)[C@@H]({{N_ACYL}}){{C1_HEAD}}"
        else:
            # Sphingosine: 4E double bond (standard biological default)
            # Natural configuration: (2S,3R).
            # In methyl-first SMILES both C3 and C2 are [C@@H] to give (2S,3R).
            # Verified: gives YDNKGFDKKRUKPY-TURZORIXSA-N for Cer d18:1/16:0.
            # Cn...C6-C5=C4-C3(OH)-C2(N)-C1(OH/head)   [4E = trans = /C=C/]
            tail_n = n_carbon - 5
            if tail_n < 0:
                raise StructureGenerationError(
                    f"Sphingoid base too short ({n_carbon}C) for d-type with double bond."
                )
            tail = "C" * tail_n
            return f"{tail}/C=C/[C@@H](O)[C@@H]({{N_ACYL}}){{C1_HEAD}}"

    elif n_oh == 1:
        # m-type: C1-OH only
        tail_n = n_carbon - 2
        tail = "C" * tail_n
        return f"{tail}[C@@H]({{N_ACYL}}){{C1_HEAD}}"

    else:
        raise UnsupportedLipidClassError(
            f"Unsupported sphingoid base: {n_carbon}:{n_db};O{n_oh}. "
            "Supported: O1 (m), O2 (d), O3 (t)."
        )


def _lcb_has_sugar(lcb_fa) -> bool:
    """Return True if the LCB has a sugar group ([X] FG), i.e. it's a glycosphingolipid."""
    raw = getattr(lcb_fa, "functional_groups", {}) or {}
    return any(str(name) == "[X]" for name in raw)


def _promote_unlocalized_mods(fa, mods: dict) -> dict:
    """Add unlocalized modifications (position <= 0) at the α-carbon (position 2).

    pygoslin reports e.g. ';OH' without a position as FG(position=-1). These would
    be silently skipped by _extract_modifications. Here we assign them to C2 (α-OH),
    the most common hydroxylation site in ceramide fatty acids.
    """
    raw = getattr(fa, "functional_groups", {}) or {}
    promoted = dict(mods)
    for name, fg_list in raw.items():
        name_upper = str(name).upper()
        if name_upper not in ("OH", "O"):
            continue
        for fg in fg_list:
            pos = getattr(fg, "position", -1)
            if pos is not None and int(pos) <= 0:
                # Unlocalized: default to α-carbon (C2)
                if 2 not in promoted:
                    promoted[2] = "OH"
    return promoted


class SphingolipidBuilder(AbstractLipidBuilder):
    """Build sphingolipid molecules from pygoslin LipidMolecule objects."""

    def build(self, lipid, lipid_name: str = "") -> Chem.Mol:
        headgroup = lipid.headgroup.headgroup
        fa_dict = lipid.fa

        lcb_fa, nacyl_fa = _extract_sp_chains(fa_dict)

        if lcb_fa is None:
            raise StructureGenerationError(
                f"Could not identify sphingoid base (LCB) in {headgroup}."
            )

        # Sphingoid base parameters
        lcb_n_carbon = int(lcb_fa.num_carbon)
        lcb_n_db = _extract_db_count(lcb_fa)
        lcb_n_oh = _get_lcb_n_oh(lcb_fa)

        # Build N-acyl part
        if nacyl_fa is not None and int(nacyl_fa.num_carbon) > 0:
            nacyl_n_carbon = int(nacyl_fa.num_carbon)
            nacyl_n_db = _extract_db_count(nacyl_fa)
            nacyl_db_pos = _extract_db_positions(nacyl_fa)
            nacyl_mods = _extract_modifications(nacyl_fa)
            # Promote unlocalized OH modifications to a default position (α-OH = pos 2)
            nacyl_mods = _promote_unlocalized_mods(nacyl_fa, nacyl_mods)

            if nacyl_n_db > 0 and not nacyl_db_pos:
                raise InsufficientStructuralDetailError(
                    f"N-acyl chain has {nacyl_n_db} double bond(s) but no positional info."
                )

            # Build acyl chain: methyl-first ending with C(=O)
            # e.g. for 16:0: "CCCCCCCCCCCCCCCC(=O)"
            acyl_fragment = build_acyl_chain(
                nacyl_n_carbon, nacyl_db_pos, nacyl_mods, terminus="ester"
            )
            # N-acyl SMILES: N bonded to C1(=O) which is at the end of acyl_fragment
            # acyl_fragment = "CCCC...CC(=O)"  (methyl→C1, C1=carbonyl)
            # We want: NC(=O)[C2..Cn]
            # The C(=O) is the last 5 chars; C2..Cn = acyl_fragment[:-5]
            if acyl_fragment.endswith("C(=O)"):
                chain_tail = acyl_fragment[:-5]  # C2..Cn (methyl end first)
                n_acyl_part = f"NC(=O){chain_tail}"
            else:
                n_acyl_part = f"NC(=O){acyl_fragment}"
        else:
            n_acyl_part = "N"  # Free sphingoid base

        # Resolve headgroup: pygoslin normalizes HexCer/GlcCer/GalCer/Hex2Cer to 'Cer'
        # with an '[X]' functional group on the LCB marking the sugar attachment.
        # We recover the original subtype from the user-supplied lipid_name because
        # pygoslin erases the distinction (GalCer vs GlcCer vs Hex2Cer all become 'Cer').
        resolved_hg = headgroup
        if headgroup == "Cer" and _lcb_has_sugar(lcb_fa):
            resolved_hg = "HexCer"  # default fallback
            if lipid_name:
                prefix = lipid_name.strip().split()[0]
                if prefix in SPHINGOLIPID_HEADGROUPS:
                    resolved_hg = prefix

        # Headgroup at C1-OH
        if resolved_hg not in SPHINGOLIPID_HEADGROUPS:
            raise UnsupportedLipidClassError(
                f"Sphingolipid headgroup {headgroup!r} is not yet supported. "
                f"Supported: {sorted(SPHINGOLIPID_HEADGROUPS.keys())}"
            )
        c1_substituent = SPHINGOLIPID_HEADGROUPS[resolved_hg]
        # C1 in the template is: C{substituent}
        # For Cer: substituent = "O"  → C1 = "CO" = primary alcohol
        # For SM:  substituent = "OP(=O)..." → C1 = "COP(=O)..."
        c1_smiles = f"C{c1_substituent}"

        # Build base template and fill placeholders
        base_template = _build_sphingoid_base_smiles(lcb_n_carbon, lcb_n_oh, lcb_n_db)
        smiles = base_template.format(N_ACYL=n_acyl_part, C1_HEAD=c1_smiles)

        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))


def _extract_sp_chains(fa_dict):
    """Extract LCB and N-acyl FA objects from pygoslin FA dict."""
    lcb = None
    nacyl = None

    for key, fa in fa_dict.items():
        key_upper = str(key).upper()
        if "LCB" in key_upper or "SPB" in key_upper:
            lcb = fa
        else:
            # "FA1" or any non-LCB key → N-acyl chain
            if nacyl is None:
                nacyl = fa

    # Fallback: if no LCB key found, first FA is LCB, second is N-acyl
    if lcb is None:
        values = list(fa_dict.values())
        if values:
            lcb = values[0]
        if len(values) > 1:
            nacyl = values[1]

    return lcb, nacyl
