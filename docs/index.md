# PyLipidParse

Convert standard lipid shorthand notation (LIPID MAPS, SwissLipids, HMDB) to molecular structures.

[![CI](https://github.com/MontgomeryBohde/PyLipidParse/actions/workflows/ci.yml/badge.svg)](https://github.com/MontgomeryBohde/PyLipidParse/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pylipidparse.svg)](https://pypi.org/project/pylipidparse/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/MontgomeryBohde/PyLipidParse/blob/main/LICENSE)

## What it does

No existing open-source tool converts lipid shorthand notation to SMILES, InChI, or RDKit molecules. PyLipidParse fills this gap by bridging [pygoslin](https://github.com/lifs-tools/pygoslin) (lipid name parsing) and [RDKit](https://www.rdkit.org/) (cheminformatics) to generate actual molecular structures.

```python
from pylipidparse import LipidConverter

conv = LipidConverter()

conv.to_smiles("PC 16:0/18:1(9Z)")
# -> canonical SMILES for POPC

conv.to_inchikey("FA 18:2(9Z,12Z)")
# -> 'OYHQOLUKZRVURQ-HZJYTTRNSA-N'

conv.to_mol("Cer 18:1;O2/16:0")
# -> <rdkit.Chem.rdchem.Mol object>
```

## Quick links

- [Installation](installation.md) -- pip, conda, or from source
- [Usage guide](usage.md) -- examples and workflows
- [Supported lipid classes](supported_classes.md) -- full list with examples
- [API Reference](api_reference.md) -- detailed method documentation
