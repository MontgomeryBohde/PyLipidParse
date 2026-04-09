"""Tests for sterol (ST) structure generation."""

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
from pylipidparse.builders.sterol import SterolBuilder
from pylipidparse.exceptions import StructureGenerationError, UnsupportedLipidClassError
from tests.conftest import assert_formula


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


class TestCholesterol:
    """Tests for free cholesterol."""

    def test_cholesterol_smiles_parseable(self, conv):
        smiles = conv.to_smiles("ST 27:1;O")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_cholesterol_formula(self, conv):
        """Cholesterol: PubChem CID 5997, formula C27H46O."""
        smiles = conv.to_smiles("ST 27:1;O")
        assert_formula(smiles, "C27H46O", "cholesterol")

    def test_fc_label(self, conv):
        """Free cholesterol (FC label)."""
        smiles = conv.to_smiles("FC")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        assert_formula(smiles, "C27H46O", "FC")


class TestCholesterolEsters:
    """Tests for cholesterol esters."""

    @pytest.mark.parametrize(
        "lipid_name,expected_formula",
        [
            ("CE 16:0", "C43H76O2"),  # Cholesteryl palmitate
            ("CE 18:1(9Z)", "C45H78O2"),  # Cholesteryl oleate
            ("CE 18:0", "C45H80O2"),  # Cholesteryl stearate
            ("CE 20:4(5Z,8Z,11Z,14Z)", "C47H76O2"),  # Cholesteryl arachidonate
        ],
    )
    def test_ce_formula(self, conv, lipid_name, expected_formula):
        smiles = conv.to_smiles(lipid_name)
        assert_formula(smiles, expected_formula, lipid_name)

    def test_ce_smiles_parseable(self, conv):
        smiles = conv.to_smiles("CE 16:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None


# ---------------------------------------------------------------------------
# Mock helpers for direct builder tests
# ---------------------------------------------------------------------------


class _MockHeadgroup:
    def __init__(self, headgroup: str):
        self.headgroup = headgroup


class _MockLipid:
    def __init__(self, headgroup: str, fa: dict = None):
        self.headgroup = _MockHeadgroup(headgroup)
        self.fa = fa or {}


# ---------------------------------------------------------------------------
# Bile acids — tested via SterolBuilder directly (pygoslin normalizes to "BA 24:1")
# Reference: PubChem
# ---------------------------------------------------------------------------


class TestBileAcids:
    """Tests for bile acid sterol structures.

    Bile acid names (CA, DCA, etc.) are parsed by pygoslin as 'BA 24:1', so the
    SterolBuilder is tested directly with mock objects to exercise the alias lookup.
    """

    def _smiles_for(self, headgroup: str) -> str:
        lipid = _MockLipid(headgroup=headgroup, fa={})
        mol = SterolBuilder().build(lipid)
        return Chem.MolToSmiles(mol)

    def test_cholic_acid_formula(self):
        """Cholic acid (CA): PubChem CID 221493, C24H40O5."""
        smiles = self._smiles_for("CA")
        assert_formula(smiles, "C24H40O5", "CA")

    def test_dca_formula(self):
        """Deoxycholic acid (DCA): PubChem CID 222528, C24H40O4."""
        smiles = self._smiles_for("DCA")
        assert_formula(smiles, "C24H40O4", "DCA")

    def test_cdca_formula(self):
        """Chenodeoxycholic acid (CDCA): PubChem CID 10133, C24H40O4."""
        smiles = self._smiles_for("CDCA")
        assert_formula(smiles, "C24H40O4", "CDCA")

    def test_udca_formula(self):
        """Ursodeoxycholic acid (UDCA): PubChem CID 31401, C24H40O4."""
        smiles = self._smiles_for("UDCA")
        assert_formula(smiles, "C24H40O4", "UDCA")

    def test_lca_formula(self):
        """Lithocholic acid (LCA): PubChem CID 9903, C24H40O3."""
        smiles = self._smiles_for("LCA")
        assert_formula(smiles, "C24H40O3", "LCA")

    def test_bile_acid_alias_cholicacid(self):
        """Alias 'CHOLICACID' maps to cholic acid."""
        smiles = self._smiles_for("CHOLICACID")
        assert_formula(smiles, "C24H40O5", "CHOLICACID")

    def test_all_bile_acids_parseable(self):
        """All bile acid SMILES must be parseable by RDKit."""
        for name in ("CA", "DCA", "CDCA", "UDCA", "LCA"):
            smiles = self._smiles_for(name)
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None, f"{name} SMILES is not parseable: {smiles!r}"


# ---------------------------------------------------------------------------
# Sterol error paths
# ---------------------------------------------------------------------------


class TestSterolErrors:
    """Error paths in SterolBuilder."""

    def test_unsupported_sterol_raises(self):
        """Unrecognized sterol headgroup raises UnsupportedLipidClassError."""
        # Must not start with ST/SE/CE/CHE and not be a bile acid name
        lipid = _MockLipid(headgroup="ERGOSTEROL", fa={})
        with pytest.raises(UnsupportedLipidClassError):
            SterolBuilder().build(lipid)

    def test_ce_no_chains_raises(self):
        """CE with no FA chain raises StructureGenerationError."""
        lipid = _MockLipid(headgroup="CE", fa={})
        with pytest.raises(StructureGenerationError, match="No fatty acid chain"):
            SterolBuilder().build(lipid)
