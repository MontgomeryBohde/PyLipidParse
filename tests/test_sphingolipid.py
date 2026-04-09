"""Tests for sphingolipid (SP) structure generation."""

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
from tests.conftest import assert_formula


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


class TestCeramideFormulas:
    """Verify molecular formulas for ceramides."""

    @pytest.mark.parametrize(
        "lipid_name,expected_formula",
        [
            # Standard ceramide d18:1/16:0 (PubChem CID: 5283564)
            ("Cer 18:1;O2/16:0", "C34H67NO3"),
            # Ceramide d18:1/24:1(15Z)
            ("Cer 18:1;O2/24:1(15Z)", "C42H81NO3"),
            # Dihydroceramide d18:0/16:0
            ("Cer 18:0;O2/16:0", "C34H69NO3"),
            # Ceramide with hydroxylated FA
            ("Cer 18:1;O2/16:0;OH", "C34H67NO4"),
        ],
    )
    def test_formula(self, conv, lipid_name, expected_formula):
        smiles = conv.to_smiles(lipid_name)
        assert_formula(smiles, expected_formula, lipid_name)

    def test_ceramide_smiles_parseable(self, conv):
        smiles = conv.to_smiles("Cer 18:1;O2/16:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None


class TestSphingomyelin:
    """Tests for sphingomyelins."""

    def test_sm_formula(self, conv):
        """SM d18:1/16:0 formula check."""
        smiles = conv.to_smiles("SM 18:1;O2/16:0")
        # SM = Cer + phosphocholine head group
        # Cer d18:1/16:0 = C34H67NO3; SM adds choline phosphate = +C5H12NO4P - H2O
        # SM d18:1/16:0 = C39H79N2O6P
        assert_formula(smiles, "C39H79N2O6P", "SM d18:1/16:0")

    def test_sm_smiles_parseable(self, conv):
        smiles = conv.to_smiles("SM 18:1;O2/16:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None


class TestGlycosphingolipids:
    """Tests for hexosylceramides."""

    def test_hexcer_parseable(self, conv):
        smiles = conv.to_smiles("HexCer 18:1;O2/16:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_hexcer_formula(self, conv):
        """HexCer d18:1/16:0 formula check."""
        # Cer d18:1/16:0 = C34H67NO3; HexCer adds hexose = +C6H10O5
        # HexCer d18:1/16:0: Cer(C34H67NO3) + hexose(C6H12O6) - H2O = C40H77NO8
        smiles = conv.to_smiles("HexCer 18:1;O2/16:0")
        assert_formula(smiles, "C40H77NO8", "HexCer d18:1/16:0")

    def test_glccer_parseable(self, conv):
        """GlcCer d18:1/16:0 — explicit glucosylceramide."""
        smiles = conv.to_smiles("GlcCer 18:1;O2/16:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_glccer_formula(self, conv):
        """GlcCer d18:1/16:0 should have the same formula as HexCer."""
        smiles = conv.to_smiles("GlcCer 18:1;O2/16:0")
        assert_formula(smiles, "C40H77NO8", "GlcCer d18:1/16:0")

    def test_galcer_parseable(self, conv):
        """GalCer d18:1/16:0 — galactosylceramide."""
        smiles = conv.to_smiles("GalCer 18:1;O2/16:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_galcer_formula(self, conv):
        """GalCer d18:1/16:0 — same formula as GlcCer (isomer)."""
        smiles = conv.to_smiles("GalCer 18:1;O2/16:0")
        assert_formula(smiles, "C40H77NO8", "GalCer d18:1/16:0")

    def test_hex2cer_formula(self, conv):
        """Hex2Cer d18:1/16:0 (lactosylceramide) — C46H87NO13.

        Two hexose sugars attached: Cer(C34H67NO3) + 2*hexose - 2*H2O = C46H87NO13.
        """
        smiles = conv.to_smiles("Hex2Cer 18:1;O2/16:0")
        assert_formula(smiles, "C46H87NO13", "Hex2Cer d18:1/16:0")
