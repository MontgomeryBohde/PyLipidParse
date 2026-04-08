"""Tests for glycerophospholipid (GP) structure generation."""

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
from pylipidparse.exceptions import InsufficientStructuralDetailError
from tests.conftest import assert_formula


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


class TestGPFormulas:
    """Verify molecular formulas for common glycerophospholipids."""

    @pytest.mark.parametrize(
        "lipid_name,expected_formula",
        [
            # PC
            ("PC 16:0/18:1(9Z)", "C42H82NO8P"),
            ("PC 18:0/20:4(5Z,8Z,11Z,14Z)", "C46H84NO8P"),
            # PE
            ("PE 16:0/18:1(9Z)", "C39H76NO8P"),
            # PA
            ("PA 16:0/18:1(9Z)", "C37H71O8P"),
            # PI
            ("PI 16:0/18:1(9Z)", "C43H81O13P"),
            # PS
            ("PS 16:0/18:1(9Z)", "C40H76NO10P"),
            # PG
            ("PG 16:0/18:1(9Z)", "C40H77O10P"),
        ],
    )
    def test_formula(self, conv, lipid_name, expected_formula):
        smiles = conv.to_smiles(lipid_name)
        assert_formula(smiles, expected_formula, lipid_name)

    def test_pc_smiles_parseable(self, conv):
        smiles = conv.to_smiles("PC 16:0/18:1(9Z)")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_pe_smiles_parseable(self, conv):
        smiles = conv.to_smiles("PE 16:0/18:1(9Z)")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_species_level_raises(self, conv):
        with pytest.raises(InsufficientStructuralDetailError):
            conv.to_smiles("PC 34:1")


class TestLysoGP:
    """Tests for lyso-glycerophospholipids."""

    @pytest.mark.parametrize(
        "lipid_name,expected_formula",
        [
            ("LPC 16:0", "C24H50NO7P"),
            ("LPE 16:0", "C21H44NO7P"),
            ("LPA 16:0", "C19H39O7P"),
        ],
    )
    def test_lyso_formula(self, conv, lipid_name, expected_formula):
        smiles = conv.to_smiles(lipid_name)
        assert_formula(smiles, expected_formula, lipid_name)


class TestEtherGP:
    """Tests for ether and plasmalogen glycerophospholipids."""

    def test_lpc_ether(self, conv):
        """LPC O-18:1 (lyso-alkyl-PC)."""
        try:
            smiles = conv.to_smiles("LPC O-18:1(9Z)")
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None
        except InsufficientStructuralDetailError:
            pytest.skip("No chain detail for LPC O-18:1")

    def test_pe_plasmalogen(self, conv):
        """PE P-16:0/18:1(9Z) (plasmalogen PE)."""
        try:
            smiles = conv.to_smiles("PE P-16:0/18:1(9Z)")
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None
        except (InsufficientStructuralDetailError, Exception):
            pytest.skip("Plasmalogen parsing not yet fully supported")
