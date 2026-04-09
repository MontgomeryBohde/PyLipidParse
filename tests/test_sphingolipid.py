"""Tests for sphingolipid (SP) structure generation."""

import pytest
from rdkit import Chem

from pylipidparse import LipidConverter
from pylipidparse.builders.sphingolipid import (
    SphingolipidBuilder,
    _build_sphingoid_base_smiles,
    _extract_sp_chains,
)
from pylipidparse.exceptions import StructureGenerationError, UnsupportedLipidClassError
from tests.conftest import assert_formula, assert_inchikey_match


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


# ---------------------------------------------------------------------------
# Mock helpers
# ---------------------------------------------------------------------------


class _MockFG:
    def __init__(self, position: int, count: int = 1):
        self.position = position
        self.count = count


class _MockFA:
    def __init__(self, num_carbon: int = 18, double_bonds=0, functional_groups=None):
        self.num_carbon = num_carbon
        self.double_bonds = double_bonds
        self.functional_groups = functional_groups or {}

    def db_num(self):
        if isinstance(self.double_bonds, dict):
            return len(self.double_bonds)
        return self.double_bonds


class _MockHeadgroup:
    def __init__(self, headgroup: str = "Cer"):
        self.headgroup = headgroup


class _MockLipid:
    def __init__(self, headgroup: str = "Cer", fa: dict = None):
        self.headgroup = _MockHeadgroup(headgroup)
        self.fa = fa or {}


# ---------------------------------------------------------------------------
# Phytosphingolipids (t-type, n_oh=3) — PubChem verified
# ---------------------------------------------------------------------------


class TestPhytosphingolipids:
    """Tests for phytosphingolipids (t-type sphingoid base, 3 OH groups)."""

    def test_phytoceramide_formula(self, conv):
        """Phytoceramide Cer 18:0;O3/16:0: PubChem CID 10506988, C34H69NO4."""
        smiles = conv.to_smiles("Cer 18:0;O3/16:0")
        assert_formula(smiles, "C34H69NO4", "Cer 18:0;O3/16:0")

    def test_phytoceramide_parseable(self, conv):
        smiles = conv.to_smiles("Cer 18:0;O3/16:0")
        assert Chem.MolFromSmiles(smiles) is not None

    def test_build_sphingoid_base_t_type(self):
        """_build_sphingoid_base_smiles with n_oh=3 (t-type / phytosphingosine)."""
        template = _build_sphingoid_base_smiles(18, n_oh=3, n_db=0)
        # Fill template placeholders to get a standalone SMILES
        smiles = template.format(N_ACYL="N", C1_HEAD="CO")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        # Phytosphingosine C18 base: C18H39NO3
        from rdkit.Chem import rdMolDescriptors
        f = rdMolDescriptors.CalcMolFormula(mol)
        assert f == "C18H39NO3", f"Expected C18H39NO3, got {f}"


# ---------------------------------------------------------------------------
# Free sphingoid bases (SPB) — PubChem verified
# ---------------------------------------------------------------------------


class TestFreeSphingoidBases:
    """Tests for free sphingoid base (SPB) lipids."""

    def test_sphingosine_spb_formula(self, conv):
        """Sphingosine SPB 18:1;O2: PubChem CID 5280335, C18H37NO2."""
        smiles = conv.to_smiles("SPB 18:1;O2")
        assert_formula(smiles, "C18H37NO2", "SPB 18:1;O2")

    def test_sphingosine_spb_inchikey(self, conv):
        """Sphingosine InChIKey: PubChem CID 5280335."""
        ik = conv.to_inchikey("SPB 18:1;O2")
        assert_inchikey_match(ik, "WWUZIQQURGPMPG-KRWOKUGFSA-N", "SPB 18:1;O2")

    def test_sphinganine_spb_formula(self, conv):
        """Sphinganine SPB 18:0;O2: PubChem CID 91486, C18H39NO2."""
        smiles = conv.to_smiles("SPB 18:0;O2")
        assert_formula(smiles, "C18H39NO2", "SPB 18:0;O2")

    def test_phytosphingosine_spb_formula(self, conv):
        """Phytosphingosine SPB 18:0;O3: PubChem CID 122121, C18H39NO3."""
        smiles = conv.to_smiles("SPB 18:0;O3")
        assert_formula(smiles, "C18H39NO3", "SPB 18:0;O3")

    def test_spb_smiles_parseable(self, conv):
        for name in ("SPB 18:1;O2", "SPB 18:0;O2", "SPB 18:0;O3"):
            smiles = conv.to_smiles(name)
            assert Chem.MolFromSmiles(smiles) is not None, f"{name} SMILES is invalid"


# ---------------------------------------------------------------------------
# Direct unit tests for internal sphingolipid functions
# ---------------------------------------------------------------------------


class TestSphingoidBaseSmiles:
    """Unit tests for _build_sphingoid_base_smiles()."""

    def test_m_type_sphingoid_base(self):
        """n_oh=1 (m-type): C1-OH only, no C3-OH."""
        template = _build_sphingoid_base_smiles(18, n_oh=1, n_db=0)
        smiles = template.format(N_ACYL="N", C1_HEAD="CO")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        # m-type 18C base: no extra OHs beyond C1 → C18H39NO
        from rdkit.Chem import rdMolDescriptors
        f = rdMolDescriptors.CalcMolFormula(mol)
        assert f == "C18H39NO", f"Expected C18H39NO, got {f}"

    def test_unsupported_n_oh_raises(self):
        """n_oh=4 is not supported — should raise UnsupportedLipidClassError."""
        with pytest.raises(UnsupportedLipidClassError):
            _build_sphingoid_base_smiles(18, n_oh=4, n_db=0)

    def test_short_sphingoid_base_with_db_raises(self):
        """Very short d-type base with double bond raises StructureGenerationError."""
        with pytest.raises(StructureGenerationError, match="too short"):
            _build_sphingoid_base_smiles(4, n_oh=2, n_db=1)


class TestExtractSpChains:
    """Unit tests for _extract_sp_chains() fallback path."""

    def test_fallback_no_lcb_key(self):
        """When no LCB/SPB key exists, first FA is LCB and second is N-acyl."""
        fa1 = _MockFA(num_carbon=18)
        fa2 = _MockFA(num_carbon=16)
        lcb, nacyl = _extract_sp_chains({"FA1": fa1, "FA2": fa2})
        assert lcb is fa1
        assert nacyl is fa2

    def test_standard_lcb_key(self):
        """Standard LCB key is identified correctly."""
        lcb_fa = _MockFA(num_carbon=18)
        nacyl_fa = _MockFA(num_carbon=16)
        lcb, nacyl = _extract_sp_chains({"LCB": lcb_fa, "FA1": nacyl_fa})
        assert lcb is lcb_fa
        assert nacyl is nacyl_fa


class TestSphingolipidBuilderErrors:
    """Error paths in SphingolipidBuilder.build()."""

    def test_no_lcb_raises(self):
        """Empty fa dict raises StructureGenerationError."""
        lipid = _MockLipid(headgroup="Cer", fa={})
        with pytest.raises(StructureGenerationError, match="sphingoid base"):
            SphingolipidBuilder().build(lipid)

    def test_unsupported_headgroup_raises(self):
        """Unsupported headgroup raises UnsupportedLipidClassError."""
        lcb = _MockFA(
            num_carbon=18,
            double_bonds={4: "E"},
            functional_groups={"O": [_MockFG(position=1), _MockFG(position=3)]},
        )
        lipid = _MockLipid(headgroup="UNKNOWNSP", fa={"LCB": lcb})
        with pytest.raises(UnsupportedLipidClassError):
            SphingolipidBuilder().build(lipid)
