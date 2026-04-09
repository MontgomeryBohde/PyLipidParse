"""Unit tests for chain SMILES generation utilities.

Tests edge cases in build_acyl_chain(), build_alkyl_chain(), and
smiles_to_mol_and_back() that are not exercised by integration tests.
"""

import warnings

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from pylipidparse.utils.chain import build_acyl_chain, build_alkyl_chain, smiles_to_mol_and_back


def formula(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles!r}"
    return rdMolDescriptors.CalcMolFormula(mol)


class TestBuildAcylChainEdgeCases:
    """Edge cases in build_acyl_chain()."""

    def test_num_carbon_zero_raises(self):
        """Zero carbons should raise ValueError."""
        with pytest.raises(ValueError, match="num_carbon must be >= 1"):
            build_acyl_chain(0)

    def test_num_carbon_negative_raises(self):
        """Negative carbons should raise ValueError."""
        with pytest.raises(ValueError):
            build_acyl_chain(-1)

    def test_unspecified_geometry_warns_and_defaults_z(self):
        """Unspecified double bond geometry should warn and default to Z."""
        with pytest.warns(UserWarning, match="unspecified geometry"):
            smiles = build_acyl_chain(18, {9: ""}, terminus="free_acid")
        # Result should be the same as explicitly Z (oleic acid)
        oleic = build_acyl_chain(18, {9: "Z"}, terminus="free_acid")
        mol1 = Chem.MolFromSmiles(smiles)
        mol2 = Chem.MolFromSmiles(oleic)
        assert mol1 is not None
        assert Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)

    def test_unspecified_geometry_none_warns(self):
        """None geometry should also warn and default to Z."""
        with pytest.warns(UserWarning):
            smiles = build_acyl_chain(18, {9: None}, terminus="free_acid")
        assert Chem.MolFromSmiles(smiles) is not None

    def test_modification_out_of_range_skipped(self):
        """Modifications at position 0 or >n should be silently skipped."""
        # Position 0 is out of range (1-indexed)
        smiles_modified = build_acyl_chain(16, modifications={0: "OH"}, terminus="free_acid")
        smiles_plain = build_acyl_chain(16, terminus="free_acid")
        mol1 = Chem.MolFromSmiles(smiles_modified)
        mol2 = Chem.MolFromSmiles(smiles_plain)
        assert mol1 is not None
        # Out-of-range mod is ignored, so result equals unmodified chain
        assert Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)

    def test_oxo_modification(self):
        """oxo (keto) modification at position 9 on C18."""
        smiles = build_acyl_chain(18, modifications={9: "oxo"}, terminus="free_acid")
        assert Chem.MolFromSmiles(smiles) is not None
        # C18 free acid + oxo at C9: C18H34O3 (keto group adds =O, replaces CH2 with C(=O))
        f = formula(smiles)
        assert f == "C18H34O3", f"Expected C18H34O3, got {f}"

    def test_me_modification(self):
        """Me (methyl branch) modification at position 5 on C16."""
        smiles = build_acyl_chain(16, modifications={5: "Me"}, terminus="free_acid")
        assert Chem.MolFromSmiles(smiles) is not None
        # C16 free acid + methyl branch at C5: C17H34O2
        f = formula(smiles)
        assert f == "C17H34O2", f"Expected C17H34O2, got {f}"

    def test_aldehyde_terminus(self):
        """Aldehyde terminus (FAL) should give C16H32O for C16."""
        smiles = build_acyl_chain(16, terminus="aldehyde")
        assert Chem.MolFromSmiles(smiles) is not None
        # Palmitaldehyde: C16H32O (PubChem CID 984)
        f = formula(smiles)
        assert f == "C16H32O", f"Expected C16H32O, got {f}"

    def test_unknown_terminus_raises(self):
        """Unknown terminus type should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown terminus type"):
            build_acyl_chain(16, terminus="bogus")

    def test_adjacent_double_bonds_stereo_conflict_warns(self):
        """Adjacent double bonds with stereo can conflict; a warning should be issued."""
        # Positions 9 and 10 are adjacent — their directional bonds share an atom,
        # which triggers the conflict warning in the stereo assignment loop.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            smiles = build_acyl_chain(18, {9: "Z", 10: "Z"}, terminus="free_acid")
        # The SMILES should still be parseable even if stereo is partially dropped
        assert Chem.MolFromSmiles(smiles) is not None
        # At least one warning should have been issued about conflicting bonds
        warn_msgs = [str(w.message) for w in caught]
        assert any("conflict" in m.lower() or "directional" in m.lower() for m in warn_msgs), (
            f"Expected a stereo conflict warning, got: {warn_msgs}"
        )


class TestBuildAlkylChainPlasmalogen:
    """Edge cases in build_alkyl_chain() for plasmalogen chains."""

    def test_plasmalogen_1_carbon_raises(self):
        """Plasmalogen chains must have at least 2 carbons."""
        with pytest.raises(ValueError, match="at least 2 carbons"):
            build_alkyl_chain(1, plasmalogen=True)

    def test_plasmalogen_2_carbon_empty_tail(self):
        """2-carbon plasmalogen: no tail carbons beyond C1=C2 vinyl ether."""
        smiles = build_alkyl_chain(2, plasmalogen=True)
        # Should return the vinyl ether portion with empty tail: /C=C\
        assert smiles == "/C=C\\"

    def test_plasmalogen_shifted_positions(self):
        """Plasmalogen with additional double bonds past C2 — positions are shifted."""
        # 18-carbon plasmalogen with additional double bond at original position 9
        # (shifted to tail position 7 after removing C1=C2)
        # Fragment ends with /C=C\ (incomplete — attaches to scaffold ether oxygen)
        fragment = build_alkyl_chain(18, double_bond_positions={9: "Z"}, plasmalogen=True)
        assert fragment.endswith("/C=C\\")
        # Should contain one more double bond from the shifted position
        assert fragment.count("=") == 2  # one from tail, one from vinyl ether C1=C2

    def test_plasmalogen_modification_shifted(self):
        """Modifications at C3+ are shifted down by 2 for the tail."""
        # Fragment ends with /C=C\ — not standalone parseable (attaches to scaffold)
        fragment = build_alkyl_chain(18, modifications={5: "OH"}, plasmalogen=True)
        assert fragment.endswith("/C=C\\")
        assert "(O)" in fragment  # OH modification present


class TestNonStandardGeometry:
    """Tests for non-standard double bond geometry strings (not Z or E)."""

    def test_non_ze_geometry_skips_stereo(self):
        """A double bond with geometry 'X' (not Z or E) skips stereo assignment."""
        # This hits the 'continue' branch at line 149 of chain.py
        smiles = build_acyl_chain(18, {9: "X"}, terminus="free_acid")
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        # Should have a double bond but no directional bonds (stereo skipped)
        assert "=" in smiles
        assert "/" not in smiles
        assert "\\" not in smiles


class TestSmilesToMolAndBack:
    """Tests for smiles_to_mol_and_back()."""

    def test_valid_smiles_returns_canonical(self):
        """Valid SMILES returns canonical SMILES string."""
        result = smiles_to_mol_and_back("CCCC")
        assert isinstance(result, str)
        assert len(result) > 0
        # Result should be canonical (parseable by RDKit)
        assert Chem.MolFromSmiles(result) is not None

    def test_valid_smiles_canonicalizes(self):
        """Different representations of the same molecule return the same result."""
        r1 = smiles_to_mol_and_back("CCCC")
        r2 = smiles_to_mol_and_back("C(C)(CC)")
        # Both are butane variants — canonical form should match
        mol1 = Chem.MolFromSmiles(r1)
        mol2 = Chem.MolFromSmiles(r2)
        assert Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)

    def test_invalid_smiles_raises_value_error(self):
        """Invalid SMILES raises ValueError."""
        with pytest.raises(ValueError, match="Invalid SMILES"):
            smiles_to_mol_and_back("not_a_smiles_string")

    def test_complex_smiles_roundtrip(self):
        """Complex SMILES (oleic acid) should canonicalize correctly."""
        result = smiles_to_mol_and_back("CCCCCCCC/C=C\\CCCCCCCC(=O)O")
        assert isinstance(result, str)
        assert Chem.MolFromSmiles(result) is not None
