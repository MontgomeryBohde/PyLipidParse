# Usage

## Basic conversion

```python
from pylipidparse import LipidConverter

conv = LipidConverter()

# SMILES
smiles = conv.to_smiles("FA 18:1(9Z)")        # Oleic acid
smiles = conv.to_smiles("PC 16:0/18:1(9Z)")   # POPC

# InChI / InChIKey
inchi = conv.to_inchi("FA 16:0")
ik = conv.to_inchikey("FA 16:0")   # "IPCSVZSSVZVIGE-UHFFFAOYSA-N"

# RDKit Mol (for downstream cheminformatics)
mol = conv.to_mol("Cer 18:1;O2/16:0")
```

## File output

```python
# Single molecule to MOL file (with 2D coordinates)
conv.to_mol_file("PC 16:0/18:1(9Z)", "popc.mol")

# Multiple molecules to SDF
conv.to_sdf(
    ["FA 16:0", "FA 18:1(9Z)", "PC 16:0/18:1(9Z)"],
    "lipids.sdf"
)
```

Each molecule in the SDF includes `_Name`, `SMILES`, and `InChIKey` properties.

## Notation requirements

PyLipidParse requires **full structural notation** with explicit chain positions.
Sum-composition notation (e.g., `PC 34:1`) is rejected because a unique structure
cannot be generated.

| Works | Fails |
|-------|-------|
| `PC 16:0/18:1(9Z)` | `PC 34:1` (no chain breakdown) |
| `TG 16:0/18:1(9Z)/18:2(9Z,12Z)` | `TG 16:0_18:1_18:2` (unknown positions) |
| `FA 18:1(9Z)` | `FA 18:1` (no double bond position) |

## Caching

`LipidConverter` caches parsed molecules (default: 512 entries). To control this:

```python
# Custom cache size
conv = LipidConverter(cache_size=1024)

# Disable caching
conv = LipidConverter(cache_size=0)

# Clear cache manually
conv.clear_cache()
```

## Error handling

```python
from pylipidparse import LipidConverter
from pylipidparse.exceptions import (
    LipidParseError,
    UnsupportedLipidClassError,
    InsufficientStructuralDetailError,
    StructureGenerationError,
)

conv = LipidConverter()

try:
    smiles = conv.to_smiles("PC 16:0/18:1(9Z)")
except LipidParseError:
    # pygoslin couldn't parse the input string
    pass
except UnsupportedLipidClassError:
    # Lipid class not yet implemented
    pass
except InsufficientStructuralDetailError:
    # Species-level input (e.g., "PC 34:1") -- no unique structure
    pass
except StructureGenerationError:
    # Molecule assembly failed (invalid chemistry)
    pass
```

## Integration with RDKit

The `to_mol()` method returns a standard RDKit `Chem.Mol` object, so you can use
any RDKit functionality downstream:

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = conv.to_mol("PC 16:0/18:1(9Z)")

# Molecular weight
mw = Descriptors.MolWt(mol)

# Fingerprints
from rdkit.Chem import AllChem
fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

# Substructure search
pattern = Chem.MolFromSmarts("[P](=O)([O-])(OCC[N+](C)(C)C)O")
has_pc_head = mol.HasSubstructMatch(pattern)
```
