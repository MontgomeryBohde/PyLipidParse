"""Tests for fatty acid variants: aldehydes, alcohols, and extraction helpers.

All reference data sourced from PubChem.
"""

import pytest

from pylipidparse import LipidConverter
from pylipidparse.builders.fatty_acid import (
    FattyAcidBuilder,
    _extract_db_count,
    _extract_modifications,
)
from pylipidparse.exceptions import StructureGenerationError
from tests.conftest import assert_formula, assert_inchikey_match


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


# ---------------------------------------------------------------------------
# Mock objects to exercise internal extraction helpers
# ---------------------------------------------------------------------------


class _MockFG:
    """Mock pygoslin FunctionalGroup object."""

    def __init__(self, position: int, count: int = 1):
        self.position = position
        self.count = count


class _MockFA:
    """Mock pygoslin FattyAcid object."""

    def __init__(self, num_carbon: int = 16, double_bonds=0, functional_groups=None):
        self.num_carbon = num_carbon
        self.double_bonds = double_bonds
        self.functional_groups = functional_groups or {}

    def db_num(self):
        if isinstance(self.double_bonds, dict):
            return len(self.double_bonds)
        return self.double_bonds


class _MockHeadgroup:
    def __init__(self, headgroup: str = "FA"):
        self.headgroup = headgroup


class _MockLipid:
    def __init__(self, headgroup: str = "FA", fa: dict = None):
        self.headgroup = _MockHeadgroup(headgroup)
        self.fa = fa or {}


# ---------------------------------------------------------------------------
# Fatty aldehyde (FAL) — PubChem CID 984
# ---------------------------------------------------------------------------


class TestFattyAldehyde:
    """Tests for fatty aldehyde (FAL) headgroup."""

    def test_fal_16_0_formula(self, conv):
        """Palmitaldehyde (FAL 16:0): PubChem CID 984, C16H32O."""
        smiles = conv.to_smiles("FAL 16:0")
        assert_formula(smiles, "C16H32O", "FAL 16:0")

    def test_fal_16_0_inchikey(self, conv):
        """Palmitaldehyde InChIKey: PubChem CID 984."""
        ik = conv.to_inchikey("FAL 16:0")
        assert_inchikey_match(ik, "NIOYUNMRJMEDGI-UHFFFAOYSA-N", "FAL 16:0")

    def test_fal_smiles_parseable(self, conv):
        """FAL 16:0 SMILES must be parseable by RDKit."""
        from rdkit import Chem

        smiles = conv.to_smiles("FAL 16:0")
        assert Chem.MolFromSmiles(smiles) is not None


# ---------------------------------------------------------------------------
# Fatty alcohol (FOH) — PubChem CID 2682
# ---------------------------------------------------------------------------


class TestFattyAlcohol:
    """Tests for fatty alcohol (FOH) headgroup."""

    def test_foh_16_0_formula(self, conv):
        """Cetyl alcohol (FOH 16:0): PubChem CID 2682, C16H34O."""
        smiles = conv.to_smiles("FOH 16:0")
        assert_formula(smiles, "C16H34O", "FOH 16:0")

    def test_foh_16_0_inchikey(self, conv):
        """Cetyl alcohol InChIKey: PubChem CID 2682."""
        ik = conv.to_inchikey("FOH 16:0")
        assert_inchikey_match(ik, "BXWNKGSJHAJOGX-UHFFFAOYSA-N", "FOH 16:0")

    def test_foh_smiles_parseable(self, conv):
        """FOH 16:0 SMILES must be parseable by RDKit."""
        from rdkit import Chem

        smiles = conv.to_smiles("FOH 16:0")
        assert Chem.MolFromSmiles(smiles) is not None


# ---------------------------------------------------------------------------
# _extract_db_count() fallback paths
# ---------------------------------------------------------------------------


class TestExtractDbCount:
    """Tests for non-callable db_num fallback paths in _extract_db_count()."""

    def test_callable_db_num(self):
        """Normal path: db_num is callable."""
        fa = _MockFA(double_bonds=2)
        assert _extract_db_count(fa) == 2

    def test_int_double_bonds_not_callable(self):
        """Fallback: fa.db_num is not callable, fa.double_bonds is an int."""

        class _FAIntDb:
            num_carbon = 18
            double_bonds = 2
            db_num = 2  # not callable

        assert _extract_db_count(_FAIntDb()) == 2

    def test_dict_double_bonds_not_callable(self):
        """Fallback: fa.db_num is not callable, fa.double_bonds is a dict."""

        class _FADictDb:
            num_carbon = 18
            double_bonds = {9: "Z", 12: "Z"}
            db_num = {9: "Z", 12: "Z"}  # not callable (a dict)

        assert _extract_db_count(_FADictDb()) == 2

    def test_fallback_zero_when_no_info(self):
        """Final fallback: return 0 when double_bonds is neither int nor dict."""

        class _FANoDb:
            num_carbon = 18
            double_bonds = None
            db_num = None  # not callable

        assert _extract_db_count(_FANoDb()) == 0


# ---------------------------------------------------------------------------
# _extract_modifications() — oxo and Me paths
# ---------------------------------------------------------------------------


class TestExtractModifications:
    """Tests for oxo/keto and Me/methyl modification extraction."""

    def test_oxo_modification(self):
        """oxo at position 9 maps to 'oxo' mod type."""
        fa = _MockFA(functional_groups={"OXO": [_MockFG(position=9)]})
        mods = _extract_modifications(fa)
        assert mods == {9: "oxo"}

    def test_keto_modification(self):
        """'KETO' name also maps to 'oxo'."""
        fa = _MockFA(functional_groups={"KETO": [_MockFG(position=5)]})
        mods = _extract_modifications(fa)
        assert mods == {5: "oxo"}

    def test_ket_modification(self):
        """'KET' name also maps to 'oxo'."""
        fa = _MockFA(functional_groups={"KET": [_MockFG(position=3)]})
        mods = _extract_modifications(fa)
        assert mods == {3: "oxo"}

    def test_me_modification(self):
        """'Me' at position 5 maps to 'Me' mod type."""
        fa = _MockFA(functional_groups={"Me": [_MockFG(position=5)]})
        mods = _extract_modifications(fa)
        assert mods == {5: "Me"}

    def test_methyl_modification(self):
        """'METHYL' name also maps to 'Me'."""
        fa = _MockFA(functional_groups={"METHYL": [_MockFG(position=7)]})
        mods = _extract_modifications(fa)
        assert mods == {7: "Me"}

    def test_unlocalized_modification_skipped(self):
        """Modifications at position <= 0 are skipped."""
        fa = _MockFA(functional_groups={"OH": [_MockFG(position=-1)]})
        mods = _extract_modifications(fa)
        assert mods == {}


# ---------------------------------------------------------------------------
# FattyAcidBuilder error paths
# ---------------------------------------------------------------------------


class TestFattyAcidBuilderErrors:
    """Error paths in FattyAcidBuilder.build()."""

    def test_no_fa_chains_raises(self):
        """Empty fa dict raises StructureGenerationError."""
        lipid = _MockLipid(headgroup="FA", fa={})
        with pytest.raises(StructureGenerationError, match="No fatty acid chains"):
            FattyAcidBuilder().build(lipid)
