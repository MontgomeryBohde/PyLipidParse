"""Integration tests for LipidConverter."""

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    LipidParseError,
)
from tests.conftest import assert_valid_inchikey


def test_import():
    from pylipidparse import LipidConverter

    assert LipidConverter is not None


def test_bad_input_raises():
    conv = LipidConverter()
    with pytest.raises(LipidParseError):
        conv.to_smiles("not a lipid string")


def test_empty_string_raises():
    conv = LipidConverter()
    with pytest.raises(LipidParseError):
        conv.to_smiles("")


def test_species_level_raises():
    """Sum-composition notation without chain details should raise."""
    conv = LipidConverter()
    with pytest.raises(InsufficientStructuralDetailError):
        conv.to_smiles("PC 34:1")


def test_to_mol_returns_mol():
    conv = LipidConverter()
    mol = conv.to_mol("FA 16:0")
    assert mol is not None
    assert isinstance(mol, Chem.Mol)
    assert mol.GetNumAtoms() > 0


def test_to_smiles_roundtrip():
    """Generated SMILES should be parseable back to a valid molecule."""
    conv = LipidConverter()
    for name in ["FA 16:0", "FA 18:1(9Z)", "PC 16:0/18:1(9Z)"]:
        smiles = conv.to_smiles(name)
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"SMILES roundtrip failed for {name}: {smiles!r}"


def test_to_inchi_returns_string():
    conv = LipidConverter()
    inchi = conv.to_inchi("FA 16:0")
    assert inchi is not None
    assert inchi.startswith("InChI=")


def test_to_inchikey_format():
    conv = LipidConverter()
    ik = conv.to_inchikey("FA 16:0")
    assert_valid_inchikey(ik, "FA 16:0")


def test_caching_returns_same_object():
    conv = LipidConverter(cache_size=10)
    mol1 = conv.to_mol("FA 16:0")
    mol2 = conv.to_mol("FA 16:0")
    # Same cached object
    assert mol1 is mol2


def test_cache_disabled():
    conv = LipidConverter(cache_size=0)
    mol1 = conv.to_mol("FA 16:0")
    mol2 = conv.to_mol("FA 16:0")
    # Not cached, but should still work
    assert Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)


def test_clear_cache():
    conv = LipidConverter(cache_size=10)
    conv.to_mol("FA 16:0")
    assert len(conv._cache) == 1
    conv.clear_cache()
    assert len(conv._cache) == 0


def test_to_mol_file(tmp_path):
    conv = LipidConverter()
    out = tmp_path / "test.mol"
    conv.to_mol_file("FA 18:1(9Z)", str(out))
    assert out.exists()
    assert out.stat().st_size > 0


def test_to_sdf_single(tmp_path):
    conv = LipidConverter()
    out = tmp_path / "single.sdf"
    conv.to_sdf("FA 16:0", str(out))
    content = out.read_text()
    assert "$$$$" in content


def test_to_sdf_batch(tmp_path):
    conv = LipidConverter()
    lipids = ["FA 16:0", "FA 18:1(9Z)", "FA 18:2(9Z,12Z)"]
    out = tmp_path / "batch.sdf"
    conv.to_sdf(lipids, str(out))
    content = out.read_text()
    assert content.count("$$$$") == 3


def test_sdf_has_properties(tmp_path):
    conv = LipidConverter()
    out = tmp_path / "props.sdf"
    conv.to_sdf("FA 16:0", str(out))
    content = out.read_text()
    assert "SMILES" in content
    assert "InChIKey" in content
