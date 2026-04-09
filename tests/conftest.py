"""Shared pytest fixtures and helpers."""

import pytest
from rdkit import Chem


@pytest.fixture(scope="session")
def converter():
    """A shared LipidConverter instance for the test session."""
    from pylipidparse import LipidConverter

    return LipidConverter()


def mols_equal(smiles1: str, smiles2: str) -> bool:
    """Return True if two SMILES represent the same molecule (by canonical form)."""
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return False
    return Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)


def assert_smiles_equivalent(generated: str, expected: str, context: str = ""):
    """Assert two SMILES represent the same molecule.

    Uses RDKit canonical SMILES comparison to handle equivalent representations.
    """
    mol_gen = Chem.MolFromSmiles(generated)
    mol_exp = Chem.MolFromSmiles(expected)

    prefix = f"[{context}] " if context else ""

    assert mol_gen is not None, f"{prefix}Generated SMILES is invalid: {generated!r}"
    assert mol_exp is not None, f"{prefix}Expected SMILES is invalid: {expected!r}"

    gen_canonical = Chem.MolToSmiles(mol_gen)
    exp_canonical = Chem.MolToSmiles(mol_exp)

    assert gen_canonical == exp_canonical, (
        f"{prefix}SMILES mismatch:\n"
        f"  Generated: {gen_canonical}\n"
        f"  Expected:  {exp_canonical}\n"
        f"  (Input generated: {generated!r})\n"
        f"  (Input expected:  {expected!r})"
    )


def assert_formula(smiles: str, expected_formula: str, context: str = ""):
    """Assert that a SMILES has the expected molecular formula."""
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    prefix = f"[{context}] " if context else ""

    assert mol is not None, f"{prefix}Invalid SMILES: {smiles!r}"

    formula = rdMolDescriptors.CalcMolFormula(mol)
    assert formula == expected_formula, (
        f"{prefix}Formula mismatch:\n" f"  Got:      {formula}\n" f"  Expected: {expected_formula}"
    )


def assert_valid_inchikey(inchikey: str, context: str = ""):
    """Assert that a string is a valid InChIKey (27 chars, correct format)."""
    prefix = f"[{context}] " if context else ""
    assert inchikey is not None, f"{prefix}InChIKey is None"
    assert (
        len(inchikey) == 27
    ), f"{prefix}InChIKey wrong length: {len(inchikey)} (expected 27): {inchikey!r}"
    assert inchikey[14] == "-", f"{prefix}InChIKey missing hyphen at pos 14: {inchikey!r}"
    assert inchikey[25] == "-", f"{prefix}InChIKey missing hyphen at pos 25: {inchikey!r}"


def assert_inchikey_match(generated: str, expected: str, context: str = ""):
    """Assert InChIKey matches a reference value, with decomposed error on failure.

    On mismatch, splits the InChIKey into its three parts to distinguish:
    - CONNECTIVITY MISMATCH: the molecular graph (atom types + bonds) is wrong
    - STEREO MISMATCH: connectivity is correct but stereochemistry differs
    - CHARGE MISMATCH: connectivity and stereo match but protonation state differs
    """
    prefix = f"[{context}] " if context else ""
    assert generated is not None, f"{prefix}InChIKey is None"
    assert len(generated) == 27, f"{prefix}Invalid InChIKey length {len(generated)}: {generated!r}"

    if generated == expected:
        return

    gen_parts = generated.split("-")
    exp_parts = expected.split("-")

    if gen_parts[0] != exp_parts[0]:
        detail = "CONNECTIVITY MISMATCH — the molecular graph is wrong"
    elif gen_parts[1] != exp_parts[1]:
        detail = "STEREO MISMATCH — connectivity is correct but stereochemistry differs"
    else:
        detail = "CHARGE MISMATCH — connectivity and stereo match but charge/protonation differs"

    raise AssertionError(
        f"{prefix}InChIKey mismatch ({detail}):\n"
        f"  Generated: {generated}\n"
        f"  Expected:  {expected}"
    )
