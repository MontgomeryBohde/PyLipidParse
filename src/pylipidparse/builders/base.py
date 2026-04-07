"""Abstract base class for all lipid builders."""
from abc import ABC, abstractmethod

from rdkit import Chem


class AbstractLipidBuilder(ABC):
    """Base class for all lipid class builders.

    Each subclass implements :meth:`build` to accept a parsed pygoslin
    lipid object and return a sanitized RDKit :class:`~rdkit.Chem.rdchem.Mol`.
    """

    @abstractmethod
    def build(self, lipid) -> Chem.Mol:
        """Build an RDKit Mol from a parsed pygoslin lipid object.

        Parameters
        ----------
        lipid :
            A parsed pygoslin ``LipidAdduct`` or ``LipidSpecies`` object.

        Returns
        -------
        Chem.Mol
            A sanitized, stereo-assigned RDKit molecule.

        Raises
        ------
        pylipidparse.exceptions.StructureGenerationError
            If the molecule cannot be assembled.
        pylipidparse.exceptions.InsufficientStructuralDetailError
            If the lipid is species-level and a unique structure cannot be generated.
        """
        ...

    @staticmethod
    def _mol_from_smiles(smiles: str) -> Chem.Mol:
        """Parse SMILES and return a sanitized Mol.

        Raises
        ------
        pylipidparse.exceptions.StructureGenerationError
            If the SMILES is invalid.
        """
        from pylipidparse.exceptions import StructureGenerationError

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise StructureGenerationError(
                f"RDKit could not parse the generated SMILES: {smiles!r}\n"
                "This is likely a bug in PyLipidParse. Please report it at "
                "https://github.com/MontgomeryBohde/PyLipidParse/issues"
            )
        return mol

    @staticmethod
    def _sanitize(mol: Chem.RWMol) -> Chem.Mol:
        """Sanitize and assign stereochemistry to a molecule.

        Raises
        ------
        pylipidparse.exceptions.StructureGenerationError
            If sanitization fails.
        """
        from rdkit.Chem import SanitizeMol, AssignStereochemistry
        from pylipidparse.exceptions import StructureGenerationError

        try:
            SanitizeMol(mol)
        except Exception as exc:
            raise StructureGenerationError(
                f"RDKit sanitization failed: {exc}\n"
                "This is likely a bug in PyLipidParse. Please report it at "
                "https://github.com/MontgomeryBohde/PyLipidParse/issues"
            ) from exc

        AssignStereochemistry(mol, cleanIt=True, force=True)
        return mol.GetMol() if hasattr(mol, "GetMol") else mol
