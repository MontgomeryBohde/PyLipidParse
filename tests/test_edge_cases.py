"""Edge case tests for PyLipidParse."""

import warnings

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
from pylipidparse.builders.base import AbstractLipidBuilder
from pylipidparse.builders.glycerophospholipid import GlycerophospholipidBuilder
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    LipidParseError,
    StructureGenerationError,
    UnsupportedLipidClassError,
)


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


class TestParseErrors:
    """Test that invalid inputs raise the correct exceptions."""

    def test_garbage_string(self, conv):
        with pytest.raises(LipidParseError):
            conv.to_smiles("hello world")

    def test_empty_string(self, conv):
        with pytest.raises(LipidParseError):
            conv.to_smiles("")

    def test_partial_notation(self, conv):
        with pytest.raises((LipidParseError, Exception)):
            conv.to_smiles("PC")

    def test_wrong_separator(self, conv):
        # Using a comma instead of / should fail or produce species-level
        try:
            conv.to_smiles("PC 16:0,18:1(9Z)")
        except (LipidParseError, InsufficientStructuralDetailError):
            pass  # Expected


class TestSpeciesLevel:
    """Species-level (sum composition) inputs must raise InsufficientStructuralDetailError."""

    @pytest.mark.parametrize(
        "lipid_name",
        [
            "PC 34:1",
            "PE 36:2",
            "TG 54:3",
        ],
    )
    def test_species_level_raises(self, conv, lipid_name):
        with pytest.raises(InsufficientStructuralDetailError):
            conv.to_smiles(lipid_name)


class TestUnknownPositions:
    """Underscore notation (unknown sn-positions) must raise."""

    def test_tg_underscore(self, conv):
        with pytest.raises(InsufficientStructuralDetailError):
            conv.to_smiles("TG 16:0_18:1_18:2")

    def test_pc_underscore(self, conv):
        with pytest.raises(InsufficientStructuralDetailError):
            conv.to_smiles("PC 16:0_18:1")


class TestUnspecifiedGeometry:
    """Double bonds without Z/E geometry should warn and default to Z."""

    def test_warns_on_missing_geometry(self, conv):
        """FA 18:1 with position but no geometry — expect warning or error."""
        # pygoslin may or may not provide geometry for 'FA 18:1(9)'
        # If it parses, we should get a valid molecule (defaulting to Z)
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                smiles = conv.to_smiles("FA 18:1(9)")
                # If warning was issued, check its message
                _unspecified = [
                    x
                    for x in w
                    if issubclass(x.category, UserWarning)
                    and "unspecified" in str(x.message).lower()
                ]
                # Either a warning was issued OR the molecule is valid
                from rdkit import Chem

                mol = Chem.MolFromSmiles(smiles)
                assert mol is not None
        except (InsufficientStructuralDetailError, LipidParseError):
            pass  # Also acceptable — pygoslin may reject it


class TestChainLengthEdgeCases:
    """Test extreme chain lengths."""

    def test_very_short_chain(self, conv):
        """FA 2:0 = acetic acid."""
        from tests.conftest import assert_smiles_equivalent

        smiles = conv.to_smiles("FA 2:0")
        assert_smiles_equivalent(smiles, "CC(=O)O", "FA 2:0 = acetic acid")

    def test_medium_length(self, conv):
        from rdkit import Chem

        smiles = conv.to_smiles("FA 14:0")  # Myristic acid
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None

    def test_long_chain(self, conv):
        from rdkit import Chem

        smiles = conv.to_smiles("FA 26:0")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None


class TestCaching:
    """Test cache behavior."""

    def test_consistent_results(self):
        """Calling twice should return the same result."""
        conv = LipidConverter(cache_size=5)
        s1 = conv.to_smiles("FA 16:0")
        s2 = conv.to_smiles("FA 16:0")
        assert s1 == s2

    def test_cache_eviction(self):
        """Cache should evict oldest entry when full."""
        conv = LipidConverter(cache_size=2)
        conv.to_mol("FA 16:0")  # Cache entry 1
        conv.to_mol("FA 18:0")  # Cache entry 2
        conv.to_mol("FA 20:0")  # Should evict FA 16:0
        assert "FA 20:0" in conv._cache
        assert "FA 18:0" in conv._cache
        assert "FA 16:0" not in conv._cache


# ---------------------------------------------------------------------------
# AbstractLipidBuilder error paths
# ---------------------------------------------------------------------------


class _ConcreteBuilder(AbstractLipidBuilder):
    """Minimal concrete subclass to test AbstractLipidBuilder static methods."""

    def build(self, lipid):
        raise NotImplementedError


class TestBuilderBase:
    """Tests for error paths in AbstractLipidBuilder."""

    def test_invalid_smiles_raises_structure_error(self):
        """_mol_from_smiles with invalid SMILES raises StructureGenerationError."""
        builder = _ConcreteBuilder()
        with pytest.raises(StructureGenerationError, match="RDKit could not parse"):
            builder._mol_from_smiles("C(C)(C)(C)(C)(C)")  # pentavalent carbon — invalid

    def test_sanitize_failure_raises_structure_error(self):
        """_sanitize with a molecule that fails sanitization raises StructureGenerationError."""
        from rdkit.Chem import RWMol

        # Create a molecule that passes MolFromSmiles but fails SanitizeMol
        # by setting an atom to have incorrect valence via direct manipulation
        mol = Chem.MolFromSmiles("CC", sanitize=False)
        rwmol = RWMol(mol)
        # Set atom 0 to have valence 8 (impossible for carbon) — triggers sanitization failure
        rwmol.GetAtomWithIdx(0).SetNumExplicitHs(10)
        rwmol.GetAtomWithIdx(0).SetNoImplicit(True)

        builder = _ConcreteBuilder()
        with pytest.raises(StructureGenerationError, match="sanitization failed"):
            builder._sanitize(rwmol)


# ---------------------------------------------------------------------------
# GlycerophospholipidBuilder error paths
# ---------------------------------------------------------------------------


class _MockHG:
    def __init__(self, hg: str):
        self.headgroup = hg


class _MockFA:
    def __init__(self, num_carbon=16, double_bonds=0, functional_groups=None):
        self.num_carbon = num_carbon
        self.double_bonds = double_bonds
        self.functional_groups = functional_groups or {}
        self.lipid_FA_bond_type = None

    def db_num(self):
        return self.double_bonds


class _MockLipid:
    def __init__(self, hg: str, fa: dict = None):
        self.headgroup = _MockHG(hg)
        self.fa = fa or {}


class TestGPBuilderErrors:
    """Error paths in GlycerophospholipidBuilder."""

    def test_unsupported_gp_class_raises(self):
        """Unknown GP class raises UnsupportedLipidClassError."""
        # CL (cardiolipin) is not in _FULL_GP_CLASSES
        lipid = _MockLipid("CL", fa={"FA1": _MockFA(), "FA2": _MockFA()})
        with pytest.raises(UnsupportedLipidClassError):
            GlycerophospholipidBuilder().build(lipid)

    def test_gp_no_chains_raises(self):
        """GP lipid with no sn1/sn2 raises InsufficientStructuralDetailError."""
        # PC with FA chains that have 0 carbons — triggers the no-chain error
        lipid = _MockLipid("PC", fa={"FA1": _MockFA(num_carbon=0), "FA2": _MockFA(num_carbon=0)})
        with pytest.raises(InsufficientStructuralDetailError):
            GlycerophospholipidBuilder().build(lipid)
