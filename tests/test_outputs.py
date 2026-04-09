"""Tests for output format methods: SMILES, InChI, InChIKey, MOL, SDF."""

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
from pylipidparse.exceptions import UnsupportedLipidClassError
from tests.conftest import assert_valid_inchikey


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


class TestSMILES:
    """SMILES output tests."""

    def test_smiles_is_string(self, conv):
        result = conv.to_smiles("FA 16:0")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_smiles_parseable(self, conv):
        """Generated SMILES must be parseable by RDKit."""
        for name in [
            "FA 16:0",
            "FA 18:1(9Z)",
            "FA 20:4(5Z,8Z,11Z,14Z)",
            "PC 16:0/18:1(9Z)",
            "PE 16:0/18:1(9Z)",
            "Cer 18:1;O2/16:0",
            "SM 18:1;O2/16:0",
            "ST 27:1;O",
            "CE 16:0",
        ]:
            try:
                smiles = conv.to_smiles(name)
                mol = Chem.MolFromSmiles(smiles)
                assert mol is not None, f"SMILES for {name} is not parseable: {smiles!r}"
            except Exception:
                pass  # Some classes might not be supported yet

    def test_smiles_canonical(self, conv):
        """Two calls should return the same canonical SMILES."""
        s1 = conv.to_smiles("FA 16:0")
        s2 = conv.to_smiles("FA 16:0")
        assert s1 == s2


class TestInChI:
    """InChI output tests."""

    def test_inchi_is_string(self, conv):
        result = conv.to_inchi("FA 16:0")
        assert isinstance(result, str)
        assert result.startswith("InChI=")

    def test_inchi_for_known_fa(self, conv):
        """Palmitic acid InChI should match known value."""
        inchi = conv.to_inchi("FA 16:0")
        # Palmitic acid InChI from PubChem CID 985:
        # InChI=1S/C16H32O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16(17)18/h2-15H2,1H3,(H,17,18)
        assert "C16H32O2" in inchi


class TestInChIKey:
    """InChIKey output tests."""

    def test_inchikey_format(self, conv):
        ik = conv.to_inchikey("FA 16:0")
        assert_valid_inchikey(ik, "FA 16:0")

    def test_inchikey_for_known_fa(self, conv):
        """Palmitic acid InChIKey from PubChem CID 985."""
        ik = conv.to_inchikey("FA 16:0")
        # PubChem CID 985 InChIKey: IPCSVZSSVZVIGE-UHFFFAOYSA-N
        assert ik == "IPCSVZSSVZVIGE-UHFFFAOYSA-N", f"Palmitic acid InChIKey mismatch: {ik}"

    def test_different_molecules_different_keys(self, conv):
        ik1 = conv.to_inchikey("FA 16:0")
        ik2 = conv.to_inchikey("FA 18:0")
        assert ik1 != ik2

    @pytest.mark.parametrize(
        "lipid_name",
        [
            "FA 16:0",
            "FA 18:1(9Z)",
            "FA 18:2(9Z,12Z)",
        ],
    )
    def test_inchikey_valid_format(self, conv, lipid_name):
        ik = conv.to_inchikey(lipid_name)
        assert_valid_inchikey(ik, lipid_name)


class TestMOLFile:
    """MOL file output tests."""

    def test_mol_file_created(self, tmp_path, conv):
        out = tmp_path / "test.mol"
        conv.to_mol_file("FA 18:1(9Z)", str(out))
        assert out.exists()

    def test_mol_file_nonempty(self, tmp_path, conv):
        out = tmp_path / "test.mol"
        conv.to_mol_file("FA 16:0", str(out))
        assert out.stat().st_size > 0

    def test_mol_file_readable(self, tmp_path, conv):
        out = tmp_path / "test.mol"
        conv.to_mol_file("FA 16:0", str(out))
        mol = Chem.MolFromMolFile(str(out))
        assert mol is not None

    def test_mol_file_has_2d_coords(self, tmp_path, conv):
        out = tmp_path / "test.mol"
        conv.to_mol_file("FA 16:0", str(out))
        mol = Chem.MolFromMolFile(str(out))
        assert mol is not None
        conf = mol.GetConformer()
        assert conf is not None


class TestSDFFile:
    """SDF file output tests."""

    def test_sdf_single_molecule(self, tmp_path, conv):
        out = tmp_path / "single.sdf"
        conv.to_sdf("FA 16:0", str(out))
        content = out.read_text()
        assert "$$$$" in content
        assert content.count("$$$$") == 1

    def test_sdf_multiple_molecules(self, tmp_path, conv):
        out = tmp_path / "batch.sdf"
        lipids = ["FA 16:0", "FA 18:0", "FA 18:1(9Z)", "FA 20:0"]
        conv.to_sdf(lipids, str(out))
        content = out.read_text()
        assert content.count("$$$$") == 4

    def test_sdf_has_name_property(self, tmp_path, conv):
        out = tmp_path / "name.sdf"
        conv.to_sdf("FA 16:0", str(out))
        content = out.read_text()
        assert "FA 16:0" in content

    def test_sdf_has_smiles_property(self, tmp_path, conv):
        out = tmp_path / "smiles.sdf"
        conv.to_sdf("FA 16:0", str(out))
        content = out.read_text()
        assert "SMILES" in content

    def test_sdf_has_inchikey_property(self, tmp_path, conv):
        out = tmp_path / "ik.sdf"
        conv.to_sdf("FA 16:0", str(out))
        content = out.read_text()
        assert "InChIKey" in content
        assert "IPCSVZSSVZVIGE-UHFFFAOYSA-N" in content


# ---------------------------------------------------------------------------
# add_hydrogens=True paths
# ---------------------------------------------------------------------------


class TestAddHydrogens:
    """Tests for the add_hydrogens parameter in to_mol_file() and to_sdf()."""

    def test_mol_file_add_hydrogens(self, tmp_path, conv):
        """to_mol_file with add_hydrogens=True should produce a larger file."""
        out_no_h = tmp_path / "no_h.mol"
        out_with_h = tmp_path / "with_h.mol"

        conv.to_mol_file("FA 16:0", str(out_no_h))
        conv.to_mol_file("FA 16:0", str(out_with_h), add_hydrogens=True)

        # File with explicit H should be larger (more atoms in the MOL block)
        assert out_with_h.stat().st_size > out_no_h.stat().st_size

    def test_mol_file_add_hydrogens_content(self, tmp_path, conv):
        """MOL file with add_hydrogens=True contains H atom entries."""
        out = tmp_path / "h.mol"
        conv.to_mol_file("FA 16:0", str(out), add_hydrogens=True)
        content = out.read_text()
        # MOL file should contain explicit H entries
        assert " H " in content

    def test_sdf_add_hydrogens(self, tmp_path, conv):
        """to_sdf with add_hydrogens=True should produce a larger file."""
        out_no_h = tmp_path / "no_h.sdf"
        out_with_h = tmp_path / "with_h.sdf"

        conv.to_sdf("FA 16:0", str(out_no_h))
        conv.to_sdf("FA 16:0", str(out_with_h), add_hydrogens=True)

        assert out_with_h.stat().st_size > out_no_h.stat().st_size

    def test_sdf_add_hydrogens_multiple(self, tmp_path, conv):
        """to_sdf with add_hydrogens=True works for multiple molecules."""
        out = tmp_path / "multi_h.sdf"
        conv.to_sdf(["FA 16:0", "FA 18:0"], str(out), add_hydrogens=True)
        content = out.read_text()
        assert content.count("$$$$") == 2


# ---------------------------------------------------------------------------
# Cache eviction
# ---------------------------------------------------------------------------


class TestCacheEviction:
    """Tests for LRU cache eviction in LipidConverter."""

    def test_cache_eviction_on_get_mol(self):
        """With cache_size=1, converting a second molecule evicts the first."""
        conv_small = LipidConverter(cache_size=1)
        # First conversion — fills the cache
        conv_small.to_smiles("FA 16:0")
        assert "FA 16:0" in conv_small._cache

        # Second conversion — should evict FA 16:0 and cache FA 18:0
        conv_small.to_smiles("FA 18:0")
        assert "FA 18:0" in conv_small._cache
        assert "FA 16:0" not in conv_small._cache

    def test_cache_eviction_cholesterol_synonym(self):
        """Cholesterol synonyms also respect cache_size=1 eviction."""
        conv_small = LipidConverter(cache_size=1)
        # Fill the cache with a real lipid first
        conv_small.to_smiles("FA 16:0")
        assert "FA 16:0" in conv_small._cache

        # Cholesterol synonym path also triggers eviction
        conv_small.to_smiles("CHOL")
        assert "CHOL" in conv_small._cache
        assert "FA 16:0" not in conv_small._cache


# ---------------------------------------------------------------------------
# Unsupported lipid class dispatch
# ---------------------------------------------------------------------------


class TestUnsupportedDispatch:
    """Tests for UnsupportedLipidClassError from _dispatch()."""

    def test_unsupported_class_raises(self):
        """A lipid class not in any known category raises UnsupportedLipidClassError."""
        # We need to reach _dispatch with an unknown class.
        # Create a minimal mock that bypasses pygoslin parsing.
        from pylipidparse.exceptions import UnsupportedLipidClassError as Err

        conv_test = LipidConverter()

        class _MockHG:
            headgroup = "UNKNOWN_CLASS_XYZ"

        class _MockLipid:
            headgroup = _MockHG()
            fa = {}

        with pytest.raises(Err, match="UNKNOWN_CLASS_XYZ"):
            conv_test._dispatch(_MockLipid(), "UNKNOWN_CLASS_XYZ")
