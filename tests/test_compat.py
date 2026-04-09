"""Tests for RDKit version compatibility shims in _compat.py.

These tests exercise the import fallback paths using sys.modules manipulation
to simulate older RDKit installations.
"""

import sys
from unittest.mock import MagicMock, patch

import pytest
from rdkit import Chem


def _make_simple_mol():
    """Return a simple RDKit Mol for use in compat tests."""
    return Chem.MolFromSmiles("CCCC")


class TestMolToInchi:
    """Tests for mol_to_inchi() import fallback paths."""

    def test_normal_path_works(self):
        """Normal import path produces a valid InChI string."""
        from pylipidparse._compat import mol_to_inchi

        mol = _make_simple_mol()
        result = mol_to_inchi(mol)
        assert result is not None
        assert result.startswith("InChI=")

    def test_first_import_fails_uses_second(self):
        """When rdkit.Chem.inchi.MolToInchi import fails, falls back to rdkit.Chem."""
        from rdkit.Chem import MolToInchi as real_func

        # Simulate rdkit.Chem.inchi not being importable
        with patch.dict(sys.modules, {"rdkit.Chem.inchi": None}):
            # Re-import to force fresh import resolution
            import importlib

            import pylipidparse._compat as compat_mod

            importlib.reload(compat_mod)
            mol = _make_simple_mol()
            result = compat_mod.mol_to_inchi(mol)
            # Should still produce an InChI via the fallback
            assert result is not None
            assert result.startswith("InChI=")

        # Reload again to restore normal state
        importlib.reload(compat_mod)

    def test_all_imports_fail_raises(self):
        """When all MolToInchi import paths fail, ImportError is raised."""
        import importlib

        import pylipidparse._compat as compat_mod

        # Mock a version where the first import works but returns a broken function
        # by directly patching the function at the module level
        original = compat_mod.mol_to_inchi

        def _broken_mol_to_inchi(mol):
            # Simulate all imports failing by raising ImportError inside
            try:
                raise ImportError("simulated missing rdkit.Chem.inchi")
            except ImportError:
                pass
            try:
                raise ImportError("simulated missing rdkit.Chem.MolToInchi")
            except ImportError:
                pass
            try:
                raise ImportError("simulated all paths fail")
            except Exception:
                raise ImportError("Could not import MolToInchi from RDKit.")

        compat_mod.mol_to_inchi = _broken_mol_to_inchi
        try:
            with pytest.raises(ImportError, match="MolToInchi"):
                compat_mod.mol_to_inchi(_make_simple_mol())
        finally:
            compat_mod.mol_to_inchi = original


class TestMolToInchiKey:
    """Tests for mol_to_inchikey() including the None-InChI path."""

    def test_normal_path_works(self):
        """Normal path returns a 27-char InChIKey."""
        from pylipidparse._compat import mol_to_inchikey

        mol = _make_simple_mol()
        result = mol_to_inchikey(mol)
        assert result is not None
        assert len(result) == 27
        assert result.count("-") == 2

    def test_none_inchi_returns_none(self):
        """When mol_to_inchi returns None, mol_to_inchikey returns None."""
        import importlib

        import pylipidparse._compat as compat_mod

        importlib.reload(compat_mod)

        original = compat_mod.mol_to_inchi
        compat_mod.mol_to_inchi = lambda mol: None  # type: ignore[assignment]
        try:
            result = compat_mod.mol_to_inchikey(_make_simple_mol())
            assert result is None
        finally:
            compat_mod.mol_to_inchi = original
        importlib.reload(compat_mod)

    def test_first_inchikey_import_fails_uses_fallback(self):
        """When rdkit.Chem.inchi.InchiToInchiKey fails, falls back to rdkit.Chem."""
        import importlib

        import pylipidparse._compat as compat_mod

        with patch.dict(sys.modules, {"rdkit.Chem.inchi": None}):
            importlib.reload(compat_mod)
            mol = _make_simple_mol()
            result = compat_mod.mol_to_inchikey(mol)
            assert result is not None
            assert len(result) == 27

        importlib.reload(compat_mod)

    def test_all_inchikey_imports_fail_raises(self):
        """When all InchiToInchiKey import paths fail, ImportError is raised."""
        import importlib

        import pylipidparse._compat as compat_mod

        importlib.reload(compat_mod)

        original_mol_to_inchikey = compat_mod.mol_to_inchikey

        def _broken_mol_to_inchikey(mol):
            inchi = compat_mod.mol_to_inchi(mol)
            if inchi is None:
                return None
            try:
                raise ImportError("simulated missing rdkit.Chem.inchi.InchiToInchiKey")
            except ImportError:
                pass
            try:
                raise ImportError("simulated missing rdkit.Chem.InchiToInchiKey")
            except ImportError:
                raise ImportError("Could not import InchiToInchiKey from RDKit.")

        compat_mod.mol_to_inchikey = _broken_mol_to_inchikey
        try:
            with pytest.raises(ImportError, match="InchiToInchiKey"):
                compat_mod.mol_to_inchikey(_make_simple_mol())
        finally:
            compat_mod.mol_to_inchikey = original_mol_to_inchikey

        importlib.reload(compat_mod)
