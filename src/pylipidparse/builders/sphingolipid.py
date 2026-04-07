"""Sphingolipid (SP) builder.

Handles ceramides (Cer), sphingomyelins (SM), and glycosphingolipids
(HexCer, Hex2Cer, GlcCer, GalCer).

Sphingolipid anatomy
--------------------
A sphingolipid consists of:
1. A **sphingoid base** (long-chain amino alcohol, e.g. sphingosine d18:1)
2. An **N-acyl chain** (fatty acid attached via amide bond to C2-NH2)
3. A **headgroup** at C1-OH (phosphocholine for SM, sugars for glycolipids)

Sphingoid base notation
-----------------------
- ``d18:1`` = dihydrosphingosine: 18C, 1 double bond, 2 OH groups (at C1 and C3)
  with the double bond at C4 in the 4E configuration (trans)
  Stereochemistry: (2S,3R,4E) — natural d-erythro form
- ``d18:0`` = dihydrosphingosine: 18C, 0 double bonds, 2 OH groups
- ``t18:0`` = phytosphingosine: 18C, 0 double bonds, 3 OH groups (C1, C3, C4)

In LIPID MAPS notation:
- ``Cer 18:1;O2/16:0`` means sphingoid base 18:1 with 2 oxygens / N-acyl 16:0
- ``Cer d18:1/16:0`` uses the legacy 'd' prefix
"""
from rdkit import Chem

from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.fatty_acid import _extract_modifications
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    StructureGenerationError,
    UnsupportedLipidClassError,
)
from pylipidparse.scaffolds.headgroups import SPHINGOLIPID_HEADGROUPS
from pylipidparse.utils.chain import build_acyl_chain


def _build_sphingoid_base_template(n_carbon: int, n_oh: int, db_pos_dict: dict) -> str:
    """Build the sphingoid base SMILES portion (without N-acyl and headgroup).

    Returns a template string with ``{N_ACYL}`` placeholder at C2-N and
    ``{C1_OH}`` placeholder at C1-OH.

    Covers the most common cases:
    - d-type (2 OH): d18:1(4E), d18:0, d20:1(4E)
    - t-type (3 OH): t18:0

    Parameters
    ----------
    n_carbon : int
        Total carbons in the sphingoid base.
    n_oh : int
        Number of hydroxyl groups (2 = d-type, 3 = t-type).
    db_pos_dict : dict
        Double bond positions (same convention as chain builder).

    Returns
    -------
    str
        SMILES with {N_ACYL} and {C1_HEAD} placeholders.
    """
    # The sphingoid base is built from the methyl end:
    # Cn ... C5 - C4=C3/... - C3(OH) - C2(NH...) - C1(OH/head)
    #
    # Tail length: C5 to Cn = (n_carbon - 4) carbons for d-type with 4E double bond
    # For d-type with 4E double bond: C4=C5 (positions 4 and 5 from C1)
    # The "tail" is C5..Cn, then /C=C\ connects C5 to C4, etc.
    #
    # Standard d-erythro-sphingosine (d18:1(4E)):
    # SMILES (methyl-first): CCCCCCCCCCCCC/C=C/[C@H](O)[C@@H]({N_ACYL}){C1_HEAD}
    # Where: C18..C6 = CCCCCCCCCCCCC (13 carbons)
    #        /C=C/   = C5=C4 (4E = trans)
    #        [C@H](O) = C3-OH (3R)
    #        [C@@H]({N_ACYL}) = C2-N (2S)
    #        {C1_HEAD} = C1 attachment

    if n_oh == 3:
        # t-type (phytosphingosine): 3 OH at C1, C3, C4; no double bond
        # Structure: C1(OH) - C2(N) - C3(OH) - C4(OH) - C5 - ... - Cn
        tail_n = n_carbon - 4
        tail = "C" * tail_n
        return (
            f"{tail}[C@@H](O)[C@H](O)[C@@H]({{N_ACYL}}){{C1_HEAD}}"
        )

    elif n_oh == 2:
        # d-type: 2 OH at C1 and C3; double bond position from db_pos_dict
        # Most common: 4E double bond
        if 4 in db_pos_dict:
            # Standard sphingosine-type: 4E double bond
            # Tail: C6..Cn (n_carbon - 5 carbons)
            tail_n = n_carbon - 5
            tail = "C" * tail_n
            geom = db_pos_dict[4]
            # 4E = trans = /C=C/ in SMILES (see chain.py notes)
            bond_geom = "/C=C/" if geom == "E" else "/C=C\\"
            return (
                f"{tail}{bond_geom}[C@H](O)[C@@H]({{N_ACYL}}){{C1_HEAD}}"
            )
        elif not db_pos_dict:
            # Dihydrosphingosine (d18:0): no double bond
            tail_n = n_carbon - 3
            tail = "C" * tail_n
            return (
                f"{tail}[C@H](O)[C@@H]({{N_ACYL}}){{C1_HEAD}}"
            )
        else:
            # Other double bond position (unusual)
            # Fall through to generic handler
            tail_n = n_carbon - 3
            tail = "C" * tail_n
            # TODO: handle non-C4 double bonds in sphingoid bases
            return (
                f"{tail}[C@H](O)[C@@H]({{N_ACYL}}){{C1_HEAD}}"
            )

    elif n_oh == 1:
        # m-type (monohydro): 1 OH at C1 only
        tail_n = n_carbon - 2
        tail = "C" * tail_n
        return f"{tail}[C@@H]({{N_ACYL}}){{C1_HEAD}}"

    else:
        raise UnsupportedLipidClassError(
            f"Unsupported sphingoid base with {n_oh} hydroxyl groups. "
            "Supported: 1 (m-type), 2 (d-type), 3 (t-type)."
        )


class SphingolipidBuilder(AbstractLipidBuilder):
    """Build sphingolipid molecules from pygoslin parsed objects."""

    def build(self, lipid) -> Chem.Mol:
        """Build a sphingolipid molecule.

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
        fa_dict = lipid.fa

        # Extract sphingoid base (LCB) and N-acyl chain
        lcb_fa, nacyl_fa = _extract_sp_chains(fa_dict)

        if lcb_fa is None:
            raise StructureGenerationError(
                f"Could not identify sphingoid base (LCB) in {headgroup} lipid."
            )

        # Build sphingoid base SMILES
        lcb_n_carbon = int(lcb_fa.num_carbon)
        lcb_n_db = int(getattr(lcb_fa, "num_double_bonds", 0))
        lcb_db_pos = _get_lcb_db_positions(lcb_fa)
        lcb_n_oh = _get_lcb_n_oh(lcb_fa)

        if lcb_n_db > 0 and not lcb_db_pos:
            raise InsufficientStructuralDetailError(
                f"Sphingoid base has {lcb_n_db} double bond(s) but no positional info. "
                "Provide full structural notation, e.g. 'Cer 18:1;O2/16:0'."
            )

        # Build N-acyl chain fragment (or handle free base)
        if nacyl_fa is not None:
            nacyl_n_carbon = int(nacyl_fa.num_carbon)
            nacyl_n_db = int(getattr(nacyl_fa, "num_double_bonds", 0))
            nacyl_db_pos = _get_nacyl_db_positions(nacyl_fa)
            nacyl_mods = _extract_modifications(nacyl_fa)

            if nacyl_n_db > 0 and not nacyl_db_pos:
                raise InsufficientStructuralDetailError(
                    f"N-acyl chain has {nacyl_n_db} double bond(s) but no positional info."
                )

            # N-acyl chain fragment: methyl-first, ends with C(=O)
            nacyl_fragment = build_acyl_chain(
                nacyl_n_carbon, nacyl_db_pos, nacyl_mods, terminus="ester"
            )
            n_acyl_smiles = f"NC(=O){nacyl_fragment[:-5]}" if nacyl_fragment.endswith("C(=O)") else f"NC(=O){nacyl_fragment}"
            # Actually, build_acyl_chain terminus='ester' ends with C(=O)
            # Fragment looks like: CCCCC...CC(=O)
            # We want: N-C(=O)-CCCCC...C (amide bond)
            # For SMILES substitution into sphingoid base template:
            # {N_ACYL} → NC(=O)[nacyl_chain_without_C1]
            # build_acyl_chain returns methyl-first chain ending at C(=O) at C1
            # So nacyl_fragment = "CCCC...CC(=O)" where the last C(=O) is C1
            # In the sphingoid base: [C@@H]({N_ACYL}) needs {N_ACYL} = N-substituent
            # {N_ACYL} = NC(=O)CCCCC...C (N followed by amide followed by C2..Cn methyl-first)
            # Note: nacyl_fragment already has C(=O) at the end (which is C1 of the acyl chain)
            # We need: NC(=O)[C2..Cn] where C2..Cn is the chain minus C1
            # nacyl_fragment = CCCC..CC(=O)
            # The C2..Cn part = nacyl_fragment[:-5] = all but the last "(=O)"... hmm
            # Actually: nacyl_fragment = "CCCCCCCCCCCCCCCC(=O)" for 16:0
            # The "CCCCCCCCCCCCCCC" part is C2-C16 (15 carbons)
            # We want: NC(=O)CCCCCCCCCCCCCCC (amide to C2-C16 from C1)
            # = "N" + "C(=O)" + C2-C16
            # nacyl_fragment already is methyl-first: C16...C2-C1(=O)
            # So: {N_ACYL} = N-C(=O)[C16..C2] but the nacyl_fragment is C16..C2-C1(=O)
            # = we can use: {N_ACYL} = "N" + nacyl_fragment (since the C(=O) is at the end)
            # [C@@H](NC(=O)CCCC...C) — YES! This works because nacyl_fragment = CCCC...CC(=O)
            # becomes [C@@H](NCCCC...CC(=O)) = [C@@H](N-C16-C15-...-C2-C1(=O)) ← WRONG ORDER
            # The amide N connects to C1 (carbonyl C), not to Cn!
            # So we need: [C@@H](NC(=O)CCCCCC) = N is bonded to C(=O) which is C1
            # The nacyl_fragment in methyl-first is: CCCCCC(=O) = C16..C2-C1(=O)
            # In SMILES for N-acyl: N-C(=O)-C2-...-C16 = NC(=O)[C2..Cn methyl direction]
            # nacyl_fragment = "CCCCCC...CC(=O)" = Cn..C2-C1(=O)
            # For the amide: N-C1(=O)-C2-...-Cn = NC(=O)C2...Cn
            # But our fragment is reversed: C16..C2-C(=O) from methyl
            # We want: NC(=O)[same fragment without the terminal C(=O)]...

            # The cleanest fix: build the nacyl chain in COOH-first order for amides
            # OR: just concatenate as N + nacyl_fragment since SMILES reads left to right
            # NC(=O)CCCCCCCCCCCCCCCC means: N bonded to C(=O) bonded to C16..C1?
            # Actually NC(=O)CCCC = N-C(=O)-C-C-C = N-amide-C3-C2-... NO
            # In SMILES, NC(=O)CCCC: N → C(=O) → C → C → C → C
            # The C(=O) here is at position 2 in the string. After C(=O), we have CCCC.
            # So N-C(=O) is the amide, and CCCC is C2..C5 (from the acyl C1 carbonyl).
            # This means nacyl_fragment needs to be: C2..Cn (methyl end last)
            # But our nacyl_fragment from build_acyl_chain(terminus='ester') is: Cn..C2-C1(=O)
            # i.e., methyl end first, carboxyl end last.
            # For an amide: N-C1(=O)-C2-...-Cn
            # We need: NC(=O)[C2..Cn] where [C2..Cn] is the non-C1 part of the chain
            # nacyl_fragment = CCCCCCCCCCCCCCCC(=O) for 16:0
            # To get C2..C16 fragment: nacyl_fragment[:-5] = "CCCCCCCCCCCCCCC" (15 carbons) ← correct!
            # Because the last 5 chars of CCCCCCCCCCCCCCCC(=O) are "C(=O)" = C1
            # So C2..C16 = nacyl_fragment[:-5]

            if nacyl_fragment.endswith("C(=O)"):
                chain_without_c1 = nacyl_fragment[:-5]  # Remove terminal C(=O) = C1
                n_acyl_part = f"NC(=O){chain_without_c1}"
            else:
                # Fallback
                n_acyl_part = f"NC(=O){nacyl_fragment}"
        else:
            n_acyl_part = "N"  # Free amine (sphingoid base only)

        # Determine headgroup attachment at C1-OH
        hg_key = headgroup
        if hg_key not in SPHINGOLIPID_HEADGROUPS:
            raise UnsupportedLipidClassError(
                f"Sphingolipid headgroup {headgroup!r} is not yet supported. "
                f"Supported: {sorted(SPHINGOLIPID_HEADGROUPS.keys())}"
            )
        c1_head = SPHINGOLIPID_HEADGROUPS[hg_key]
        # c1_head is the group at C1-OH; for Cer it's just 'O' (free OH)
        # For SM it's 'OP(=O)([O-])OCC[N+](C)(C)C'
        # The C1 atom in the template is: C{C1_HEAD} = CO (for Cer) or C-OP... (for SM)
        c1_smiles = f"C{c1_head}"  # The C1 carbon with its substituent

        # Build the complete sphingoid base + modification
        base_template = _build_sphingoid_base_template(
            lcb_n_carbon, lcb_n_oh, lcb_db_pos
        )

        # Fill in the template
        smiles = base_template.format(N_ACYL=n_acyl_part, C1_HEAD=c1_smiles)

        mol = self._mol_from_smiles(smiles)
        return self._sanitize(Chem.RWMol(mol))


def _get_lcb_db_positions(lcb_fa) -> dict:
    """Extract double bond positions from the sphingoid base FA."""
    db_positions = {}
    raw_db = getattr(lcb_fa, "double_bond_positions", None)
    if raw_db:
        if hasattr(raw_db, "get_positions"):
            try:
                for pos, geom in raw_db.get_positions().items():
                    db_positions[int(pos)] = str(geom) if geom else "E"
            except AttributeError:
                pass
        elif hasattr(raw_db, "items"):
            for pos, geom in raw_db.items():
                db_positions[int(pos)] = str(geom) if geom else "E"
    return db_positions


def _get_nacyl_db_positions(nacyl_fa) -> dict:
    """Extract double bond positions from the N-acyl FA."""
    db_positions = {}
    raw_db = getattr(nacyl_fa, "double_bond_positions", None)
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
    return db_positions


def _get_lcb_n_oh(lcb_fa) -> int:
    """Determine the number of hydroxyl groups on the sphingoid base."""
    # pygoslin represents this via num_hydroxyl or through ;O2 notation
    n_oh = getattr(lcb_fa, "num_hydroxyl", None)
    if n_oh is not None and int(n_oh) > 0:
        return int(n_oh)

    # Try to infer from modifications
    mods = _extract_modifications(lcb_fa)
    n_oh_mods = sum(1 for m in mods.values() if m == "OH")
    if n_oh_mods > 0:
        return n_oh_mods

    # Default: d-type (2 OH) is standard for most ceramides
    return 2


def _extract_sp_chains(fa_dict):
    """Extract LCB (sphingoid base) and N-acyl FA from pygoslin FA dict.

    pygoslin uses "LCB" key for the sphingoid base and "FA1" for the N-acyl chain.
    """
    lcb = None
    nacyl = None

    for key, fa in fa_dict.items():
        key_upper = str(key).upper()
        if "LCB" in key_upper or key_upper == "SPB":
            lcb = fa
        elif "FA" in key_upper or key_upper in ("NACYL", "FA1"):
            nacyl = fa
        else:
            # Fallback assignment
            if lcb is None:
                lcb = fa
            else:
                nacyl = fa

    return lcb, nacyl
