"""Tests for sterol (ST) structure generation."""

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
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
