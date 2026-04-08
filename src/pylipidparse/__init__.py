"""PyLipidParse: Convert lipid shorthand notation to molecular structures.

Usage
-----
>>> from pylipidparse import LipidConverter
>>> conv = LipidConverter()
>>> conv.to_smiles("FA 18:1(9Z)")
'CCCCCCCC/C=C\\\\CCCCCCCC(=O)O'
>>> conv.to_inchikey("PC 16:0/18:1(9Z)")
'...'
"""

from pylipidparse.converter import LipidConverter
from pylipidparse.exceptions import (
    InsufficientStructuralDetailError,
    LipidParseError,
    PyLipidParseError,
    StructureGenerationError,
    UnsupportedLipidClassError,
)

__version__ = "0.1.0"
__all__ = [
    "LipidConverter",
    "PyLipidParseError",
    "LipidParseError",
    "StructureGenerationError",
    "UnsupportedLipidClassError",
    "InsufficientStructuralDetailError",
]
