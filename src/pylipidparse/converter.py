"""Main LipidConverter class — the public entry point for PyLipidParse."""
import warnings
from typing import List, Optional, Union

from rdkit import Chem

from pylipidparse._compat import compute_2d_coords, mol_to_inchi, mol_to_inchikey
from pylipidparse.exceptions import (
    LipidParseError,
    UnsupportedLipidClassError,
)

# Lipid class dispatch tables
_FA_CLASSES = frozenset({
    "FA", "FAL", "FOH", "FAM", "FAO", "FAHFA",
    # Oxidized fatty acids
    "OXF", "OXFA",
})

_GL_CLASSES = frozenset({
    "MG", "DG", "TG",
    "MAG", "DAG", "TAG",  # alternative naming
    "MG_", "DG_", "TG_",
})

_GP_CLASSES = frozenset({
    "PC", "PE", "PA", "PI", "PS", "PG",
    "PIP", "PIP2", "PIP3",
    "LPC", "LPE", "LPA", "LPI", "LPS", "LPG",
    "LPIP", "LPIP2", "LPIP3",
    "LNPC", "LNPE",
    # Oxidized variants
    "OxPC", "OxPE", "OxPA", "OxPI",
    "OXPC", "OXPE", "OXPA", "OXPI",
})

_SP_CLASSES = frozenset({
    "Cer", "CerP", "Cer_P",
    "HexCer", "Hex1Cer", "Hex2Cer", "Hex3Cer",
    "SM", "LSM",
    "GlcCer", "GalCer", "LacCer",
    "SPB", "SPBP",
    "SHexCer", "SHex2Cer",
    "Cer1P",
    "CER",
    # Gangliosides — future support
})

_ST_CLASSES = frozenset({
    "ST", "FC", "CE", "ChE",
    "Cholesterol",
    "BA", "CA", "DCA", "CDCA", "UDCA", "LCA",
})


class LipidConverter:
    """Convert lipid shorthand notation to molecular structures.

    This is the main public class. It wraps a cached pygoslin parser and
    dispatches to the appropriate structure builder for each lipid class.

    Parameters
    ----------
    dialect : str
        Parsing dialect for pygoslin. Options: ``"LipidMaps"`` (default),
        ``"Goslin"``, ``"SwissLipids"``, ``"HMDB"``.
    cache_size : int
        Maximum number of parsed molecules to keep in the LRU cache.
        Set to 0 to disable caching (useful for memory-constrained environments).

    Examples
    --------
    >>> from pylipidparse import LipidConverter
    >>> conv = LipidConverter()
    >>> conv.to_smiles("FA 18:1(9Z)")
    'CCCCCCCC/C=C\\\\CCCCCCCC(=O)O'
    >>> conv.to_inchikey("PC 16:0/18:1(9Z)")
    'XXXXXXXXXXXXXX-XXXXXXXXXX-N'
    """

    def __init__(self, dialect: str = "LipidMaps", cache_size: int = 512):
        self._dialect = dialect
        self._cache_size = cache_size
        self._cache: dict = {}
        self._parser = None  # lazy init

    def _get_parser(self):
        """Lazy-initialize the pygoslin parser."""
        if self._parser is None:
            try:
                from pygoslin.parser.Parser import LipidParser
            except ImportError as exc:
                raise ImportError(
                    "pygoslin is required for PyLipidParse. "
                    "Install it with: pip install pygoslin"
                ) from exc
            self._parser = LipidParser()
        return self._parser

    def _parse(self, lipid_name: str):
        """Parse a lipid name with pygoslin."""
        parser = self._get_parser()
        lipid = parser.parse(lipid_name.strip())
        if lipid is None:
            raise LipidParseError(
                f"Could not parse lipid: {lipid_name!r}\n"
                "Ensure the name follows LIPID MAPS shorthand notation, "
                "e.g. 'PC 16:0/18:1(9Z)', 'FA 18:2(9Z,12Z)', 'Cer 18:1;O2/16:0'."
            )
        return lipid

    def _get_mol(self, lipid_name: str) -> Chem.Mol:
        """Get (or build and cache) the RDKit Mol for a lipid name."""
        if lipid_name in self._cache:
            return self._cache[lipid_name]

        lipid = self._parse(lipid_name)
        mol = self._dispatch(lipid)

        # LRU eviction: if cache is full, remove the oldest entry
        if self._cache_size > 0:
            if len(self._cache) >= self._cache_size:
                # Remove the oldest key (insertion order preserved in Python 3.7+)
                oldest_key = next(iter(self._cache))
                del self._cache[oldest_key]
            self._cache[lipid_name] = mol

        return mol

    def _dispatch(self, lipid) -> Chem.Mol:
        """Dispatch to the correct builder based on lipid class."""
        headgroup = lipid.headgroup.headgroup
        hg_upper = headgroup.upper()

        if hg_upper in {c.upper() for c in _FA_CLASSES}:
            from pylipidparse.builders.fatty_acid import FattyAcidBuilder
            return FattyAcidBuilder().build(lipid)

        if hg_upper in {c.upper() for c in _GL_CLASSES}:
            from pylipidparse.builders.glycerolipid import GlycerolipidBuilder
            return GlycerolipidBuilder().build(lipid)

        if hg_upper in {c.upper() for c in _GP_CLASSES}:
            from pylipidparse.builders.glycerophospholipid import GlycerophospholipidBuilder
            return GlycerophospholipidBuilder().build(lipid)

        if hg_upper in {c.upper() for c in _SP_CLASSES}:
            from pylipidparse.builders.sphingolipid import SphingolipidBuilder
            return SphingolipidBuilder().build(lipid)

        if hg_upper in {c.upper() for c in _ST_CLASSES}:
            from pylipidparse.builders.sterol import SterolBuilder
            return SterolBuilder().build(lipid)

        raise UnsupportedLipidClassError(
            f"Lipid class {headgroup!r} is not yet supported by PyLipidParse.\n"
            "Supported classes: Fatty Acids (FA), Glycerolipids (MG/DG/TG), "
            "Glycerophospholipids (PC/PE/PA/PI/PS/PG and lyso variants), "
            "Sphingolipids (Cer/SM/HexCer), Sterols (ST/CE/FC/bile acids).\n"
            "Please open an issue at https://github.com/MontgomeryBohde/PyLipidParse/issues"
        )

    def to_mol(self, lipid_name: str) -> Chem.Mol:
        """Convert a lipid name to an RDKit Mol object.

        The returned molecule is sanitized and has stereochemistry assigned,
        but does NOT have 2D coordinates. Use :meth:`to_mol_file` or
        :meth:`to_sdf` for a molecule with 2D layout.

        Parameters
        ----------
        lipid_name : str
            Lipid shorthand notation, e.g. ``'PC 16:0/18:1(9Z)'``.

        Returns
        -------
        Chem.Mol
            Sanitized RDKit molecule.

        Raises
        ------
        LipidParseError
            If the name cannot be parsed.
        UnsupportedLipidClassError
            If the lipid class is not yet supported.
        InsufficientStructuralDetailError
            If the input is species-level (sum composition).
        StructureGenerationError
            If molecule assembly fails.
        """
        return self._get_mol(lipid_name)

    def to_smiles(self, lipid_name: str) -> str:
        """Convert a lipid name to a canonical SMILES string.

        Parameters
        ----------
        lipid_name : str
            Lipid shorthand notation.

        Returns
        -------
        str
            Canonical SMILES string (RDKit canonical form).
        """
        return Chem.MolToSmiles(self.to_mol(lipid_name))

    def to_inchi(self, lipid_name: str) -> Optional[str]:
        """Convert a lipid name to an InChI string.

        Parameters
        ----------
        lipid_name : str
            Lipid shorthand notation.

        Returns
        -------
        str or None
            InChI string, or None if InChI generation failed.
        """
        return mol_to_inchi(self.to_mol(lipid_name))

    def to_inchikey(self, lipid_name: str) -> Optional[str]:
        """Convert a lipid name to an InChIKey string.

        Parameters
        ----------
        lipid_name : str
            Lipid shorthand notation.

        Returns
        -------
        str or None
            InChIKey string (27-character hash), or None if generation failed.
        """
        return mol_to_inchikey(self.to_mol(lipid_name))

    def to_mol_file(
        self,
        lipid_name: str,
        path: str,
        add_hydrogens: bool = False,
    ) -> None:
        """Write a lipid structure to a MOL file (with 2D coordinates).

        Parameters
        ----------
        lipid_name : str
            Lipid shorthand notation.
        path : str
            Output file path (should end in ``.mol``).
        add_hydrogens : bool
            If True, add explicit hydrogens before writing.
        """
        mol = Chem.RWMol(self.to_mol(lipid_name))
        if add_hydrogens:
            from rdkit.Chem import AddHs
            mol = Chem.RWMol(AddHs(mol))
        compute_2d_coords(mol)
        Chem.MolToMolFile(mol, path)

    def to_sdf(
        self,
        lipid_names: Union[str, List[str]],
        path: str,
        add_hydrogens: bool = False,
    ) -> None:
        """Write one or more lipid structures to an SDF file.

        Each molecule in the SDF has the following properties:
        - ``_Name``: lipid shorthand name
        - ``SMILES``: canonical SMILES
        - ``InChIKey``: InChIKey

        Parameters
        ----------
        lipid_names : str or list of str
            Lipid shorthand notation(s).
        path : str
            Output file path (should end in ``.sdf``).
        add_hydrogens : bool
            If True, add explicit hydrogens before writing.
        """
        from rdkit.Chem import SDWriter

        if isinstance(lipid_names, str):
            lipid_names = [lipid_names]

        writer = SDWriter(path)
        try:
            for name in lipid_names:
                mol = Chem.RWMol(self.to_mol(name))
                if add_hydrogens:
                    from rdkit.Chem import AddHs
                    mol = Chem.RWMol(AddHs(mol))
                compute_2d_coords(mol)
                mol.SetProp("_Name", name)
                mol.SetProp("SMILES", Chem.MolToSmiles(mol))
                ik = mol_to_inchikey(mol)
                if ik:
                    mol.SetProp("InChIKey", ik)
                writer.write(mol)
        finally:
            writer.close()

    def clear_cache(self) -> None:
        """Clear the molecule cache."""
        self._cache.clear()
