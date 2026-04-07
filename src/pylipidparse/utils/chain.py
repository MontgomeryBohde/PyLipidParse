"""Fatty acid and alkyl chain SMILES generation.

This module is the mathematical core of PyLipidParse. It converts
parsed chain specifications (carbon count, double bond positions,
modifications, linkage type) into SMILES fragments.

Chain SMILES are always written in **methyl-first** order (Cn → C1),
which matches the convention used by most chemical databases.

Terminology
-----------
- C1: the carboxyl/attachment-end carbon
- Cn: the methyl-end carbon
- position p: 1-indexed from the carboxyl end (C1=1, Cn=n)
- double bond at position p: bond between C_p and C_{p+1}
- Z: cis configuration (same-side substituents)
- E: trans configuration (opposite-side substituents)
"""
import warnings
from typing import Dict, Optional


def build_acyl_chain(
    num_carbon: int,
    double_bond_positions: Optional[Dict[int, str]] = None,
    modifications: Optional[Dict[int, str]] = None,
    terminus: str = "ester",
) -> str:
    """Build an acyl (or alkyl) chain SMILES fragment in methyl-first format.

    The chain runs from C_n (methyl end) to C_1 (attachment end).
    The SMILES fragment does NOT include the ester/ether oxygen — that is
    provided by the scaffold. The fragment ends at the attachment carbon C_1.

    Parameters
    ----------
    num_carbon : int
        Total carbon count including C1 (e.g., 16 for palmitic acid).
    double_bond_positions : dict, optional
        Mapping of ``{position: geometry}`` where position is 1-indexed
        (the bond is between C_position and C_{position+1}) and geometry
        is ``'Z'`` (cis), ``'E'`` (trans), or ``''`` (unspecified → Z with warning).
    modifications : dict, optional
        Mapping of ``{position: mod_type}`` where position is 1-indexed.
        Supported mod_types: ``'OH'``, ``'oxo'``, ``'Me'``.
    terminus : str
        How to represent C1:
        - ``'ester'``:      C1 is ``C(=O)``  — for ester linkage to a scaffold
        - ``'free_acid'``:  C1 is ``C(=O)O`` — free carboxylic acid
        - ``'alkyl'``:      C1 is ``C``       — for ether linkage (no carbonyl)
        - ``'aldehyde'``:   C1 is ``C=O``     — fatty aldehyde
        - ``'amide'``:      C1 is ``C(=O)``   — for N-acyl bond (same as ester)

    Returns
    -------
    str
        SMILES fragment in methyl-first format.

    Examples
    --------
    >>> build_acyl_chain(16, terminus='free_acid')
    'CCCCCCCCCCCCCCCC(=O)O'
    >>> build_acyl_chain(18, {9: 'Z'}, terminus='free_acid')
    'CCCCCCCC/C=C\\\\CCCCCCCC(=O)O'
    """
    if num_carbon < 1:
        raise ValueError(f"num_carbon must be >= 1, got {num_carbon}")

    db_pos = double_bond_positions or {}
    mods = modifications or {}

    # Warn on unspecified geometry; default to Z (biologically most common)
    for pos, geom in list(db_pos.items()):
        if geom == "" or geom is None:
            warnings.warn(
                f"Double bond at position {pos} has unspecified geometry (Z or E). "
                "Defaulting to Z (cis). Specify explicitly to suppress this warning.",
                UserWarning,
                stacklevel=4,
            )
            db_pos[pos] = "Z"

    # Validate double bond positions
    for pos in db_pos:
        if not (1 <= pos <= num_carbon - 1):
            raise ValueError(
                f"Double bond position {pos} is out of range for a {num_carbon}-carbon chain. "
                f"Valid range: 1..{num_carbon - 1}"
            )

    # ------------------------------------------------------------------ #
    # Build atom tokens and bond tokens (1-indexed, C1=carboxyl end)      #
    # ------------------------------------------------------------------ #
    # atom_tokens[i]: SMILES token for C_i
    # bond_tokens[i]: bond between C_i and C_{i+1}  (for i in 1..n-1)
    atom_tokens: Dict[int, str] = {i: "C" for i in range(1, num_carbon + 1)}
    bond_tokens: Dict[int, str] = {i: "" for i in range(1, num_carbon)}

    # Apply modifications to internal atoms (C1 terminus overrides these)
    for pos, mod in mods.items():
        if not (1 <= pos <= num_carbon):
            continue
        if mod == "OH":
            atom_tokens[pos] = "C(O)"
        elif mod == "oxo":
            atom_tokens[pos] = "C(=O)"
        elif mod == "Me":
            atom_tokens[pos] = "C(C)"
        # Additional mods can be added here

    # Apply C1 terminus (overrides any modification at position 1)
    if terminus in ("ester", "amide"):
        atom_tokens[1] = "C(=O)"
    elif terminus == "free_acid":
        atom_tokens[1] = "C(=O)O"
    elif terminus == "alkyl":
        atom_tokens[1] = "C"
    elif terminus == "aldehyde":
        atom_tokens[1] = "C=O"
    else:
        raise ValueError(f"Unknown terminus type: {terminus!r}")

    # Apply double bonds to bond_tokens
    for pos, geom in db_pos.items():
        bond_tokens[pos] = "="

    # Apply stereochemistry directional bonds for Z/E.
    # Traversal direction: methyl-first (C_n → C_1).
    # For a Z double bond at position p (bond between C_p and C_{p+1}):
    #   - bond[p+1] (between C_{p+1} and C_{p+2}) gets '/' in methyl-first
    #   - bond[p-1] (between C_{p-1} and C_p) gets '\' for Z, '/' for E
    #
    # Verified against PubChem SMILES for oleic acid (9Z):
    #   CCCCCCCC/C=C\CCCCCCCC(=O)O  — the /C=C\ pattern = Z ✓
    for pos, geom in db_pos.items():
        if geom not in ("Z", "E"):
            continue

        # Bond entering C_{p+1} from C_{p+2} (in methyl-first traversal)
        if pos + 1 < num_carbon:  # C_{p+2} exists
            existing = bond_tokens.get(pos + 1, "")
            if existing and existing not in ("/", "\\"):
                # Another double bond's stereo already set this bond — conflict
                # (e.g., adjacent double bonds). Fall back to no stereo for safety.
                warnings.warn(
                    f"Conflicting directional bond at position {pos + 1}. "
                    "Stereochemistry for adjacent/overlapping double bonds may be incorrect.",
                    UserWarning,
                    stacklevel=4,
                )
            elif not existing:
                bond_tokens[pos + 1] = "/"

        # Bond leaving C_p toward C_{p-1} (in methyl-first traversal)
        if pos > 1:  # C_{p-1} exists
            existing = bond_tokens.get(pos - 1, "")
            after_dir = "\\" if geom == "Z" else "/"
            if existing and existing not in ("/", "\\"):
                warnings.warn(
                    f"Conflicting directional bond at position {pos - 1}. "
                    "Stereochemistry for adjacent/overlapping double bonds may be incorrect.",
                    UserWarning,
                    stacklevel=4,
                )
            elif not existing:
                bond_tokens[pos - 1] = after_dir

    # ------------------------------------------------------------------ #
    # Assemble SMILES in methyl-first order: C_n, bond[n-1], ..., C_1    #
    # ------------------------------------------------------------------ #
    parts = []
    for i in range(num_carbon, 0, -1):
        parts.append(atom_tokens[i])
        if i > 1:
            # Bond between C_{i-1} and C_i (appears in methyl-first between C_i and C_{i-1})
            parts.append(bond_tokens.get(i - 1, ""))

    return "".join(parts)


def build_alkyl_chain(
    num_carbon: int,
    double_bond_positions: Optional[Dict[int, str]] = None,
    modifications: Optional[Dict[int, str]] = None,
    plasmalogen: bool = False,
) -> str:
    """Build an alkyl/alkenyl chain fragment for ether linkages.

    Parameters
    ----------
    num_carbon : int
        Total carbon count (C1 is the ether carbon).
    double_bond_positions : dict, optional
        Double bond positions (same convention as :func:`build_acyl_chain`).
        Note: for plasmalogens the C1=C2 double bond is handled automatically
        via the ``plasmalogen`` flag; do not include it here.
    modifications : dict, optional
        Modifications (same convention as :func:`build_acyl_chain`).
    plasmalogen : bool
        If True, adds the C1=C2 alk-1-enyl double bond (Z configuration)
        characteristic of plasmalogen lipids. The chain fragment will
        start with ``/C=C\\`` (connecting to the scaffold ether oxygen).

    Returns
    -------
    str
        SMILES fragment. For ether: ends with ``C`` at C1.
        For plasmalogen: starts at C1 with the vinyl ether portion
        (to be used with a scaffold that has the ether O already).
    """
    if plasmalogen:
        # C1=C2 double bond (alk-1-enyl, always Z in biological plasmalogens)
        # The scaffold provides the ether oxygen; this fragment starts at C1.
        # Chain fragment: C3..Cn (methyl-first) + /C=C\ (C2=C1) ending at C1
        # The ether O in the scaffold connects to C1.
        if num_carbon < 2:
            raise ValueError("Plasmalogen chain must have at least 2 carbons.")

        # Build C3..Cn tail (n-2 carbons, no C1=C2 included in db_positions)
        tail_n = num_carbon - 2
        tail_db = {}
        tail_mods = {}
        for pos, geom in (double_bond_positions or {}).items():
            if pos >= 2:  # positions in original numbering, shift by -2 for tail
                tail_db[pos - 2] = geom
        for pos, mod in (modifications or {}).items():
            if pos >= 3:
                tail_mods[pos - 2] = mod

        if tail_n > 0:
            tail = build_acyl_chain(
                tail_n, tail_db, tail_mods, terminus="alkyl"
            )
        else:
            tail = ""

        # Plasmalogen SMILES: [tail_Cn..C3]/C=C\ where C2=C1 is the alk-1-enyl part
        # In methyl-first: Cn...C3 (tail) / C2=C1 \
        # The trailing \ is the bond from C1 to the scaffold ether O
        # When used in scaffold [SN]OC... the substitution gives: .../C=C\OC... ← vinyl ether
        # SMILES: tail + /C=C\ where the / bond is between C3 and C2, \ is between C1 and O
        return tail + "/C=C\\"

    else:
        return build_acyl_chain(
            num_carbon,
            double_bond_positions,
            modifications,
            terminus="alkyl",
        )


def smiles_to_mol_and_back(smiles: str) -> str:
    """Parse a SMILES string and return the canonical SMILES.

    Used for validation and normalization.

    Raises
    ------
    ValueError
        If the SMILES is invalid.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    return Chem.MolToSmiles(mol)
