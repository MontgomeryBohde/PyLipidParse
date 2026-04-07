"""Tests for glycerolipid (GL) structure generation: MAG, DAG, TAG."""
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from pylipidparse import LipidConverter
from pylipidparse.exceptions import InsufficientStructuralDetailError
from tests.conftest import assert_formula, assert_smiles_equivalent


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


class TestGlycerolipidFormulas:
    """Verify molecular formulas for common glycerolipids."""

    @pytest.mark.parametrize("lipid_name,expected_formula", [
        # MAG
        ("MG 16:0/0:0/0:0", "C19H38O4"),
        ("MG 18:1(9Z)/0:0/0:0", "C21H40O4"),
        # DAG
        ("DG 16:0/18:1(9Z)/0:0", "C37H70O5"),
        ("DG 16:0/0:0/18:1(9Z)", "C37H70O5"),
        # TAG
        ("TG 16:0/18:1(9Z)/18:2(9Z,12Z)", "C55H100O6"),
        ("TG 16:0/16:0/16:0", "C51H98O6"),
        ("TG 18:0/18:1(9Z)/18:2(9Z,12Z)", "C57H104O6"),
    ])
    def test_formula(self, conv, lipid_name, expected_formula):
        smiles = conv.to_smiles(lipid_name)
        assert_formula(smiles, expected_formula, lipid_name)

    def test_tag_smiles_parseable(self, conv):
        smiles = conv.to_smiles("TG 16:0/18:1(9Z)/18:2(9Z,12Z)")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_dag_smiles_parseable(self, conv):
        smiles = conv.to_smiles("DG 16:0/18:1(9Z)/0:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_mag_smiles_parseable(self, conv):
        smiles = conv.to_smiles("MG 16:0/0:0/0:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_species_level_raises(self, conv):
        """TAG with unknown positions should raise."""
        with pytest.raises(InsufficientStructuralDetailError):
            conv.to_smiles("TG 16:0_18:1_18:2")

    def test_ether_tag(self, conv):
        """Ether TAG (O- prefix at sn-1)."""
        try:
            smiles = conv.to_smiles("TG O-16:0/18:1(9Z)/18:2(9Z,12Z)")
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None
        except InsufficientStructuralDetailError:
            pytest.skip("pygoslin did not provide chain detail for this ether TAG")
