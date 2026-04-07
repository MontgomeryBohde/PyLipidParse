"""Custom exceptions for PyLipidParse."""


class PyLipidParseError(Exception):
    """Base exception for all PyLipidParse errors."""


class LipidParseError(PyLipidParseError):
    """Raised when pygoslin cannot parse the input lipid string."""


class StructureGenerationError(PyLipidParseError):
    """Raised when the molecular structure cannot be assembled (invalid chemistry)."""


class UnsupportedLipidClassError(PyLipidParseError):
    """Raised when the lipid class is recognized but not yet supported."""


class InsufficientStructuralDetailError(PyLipidParseError):
    """Raised when the input is species-level (sum composition) and a unique
    structure cannot be generated.

    Example: 'PC 34:1' without explicit chain breakdown cannot produce a unique SMILES.
    Use full structural notation like 'PC 16:0/18:1(9Z)' instead.
    """
