# PyLipidParse

Convert standard lipid shorthand notation (LIPID MAPS, SwissLipids, HMDB) to molecular structures.

[![CI](https://github.com/MontgomeryBohde/PyLipidParse/actions/workflows/ci.yml/badge.svg)](https://github.com/MontgomeryBohde/PyLipidParse/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pylipidparse.svg)](https://pypi.org/project/pylipidparse/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## What it does

No existing open-source tool converts lipid shorthand notation to SMILES, InChI, or RDKit molecules. PyLipidParse fills this gap.

```python
from pylipidparse import LipidConverter

conv = LipidConverter()

conv.to_smiles("PC 16:0/18:1(9Z)")
# → 'CCCCCCCCCCCCCCCC(=O)OC[C@@H](OC(=O)CCCCCCCC/C=C\\CCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C'

conv.to_inchikey("FA 18:2(9Z,12Z)")
# → 'OYHQOLUKZRVURQ-HZJYTTRNSA-N'

conv.to_mol("Cer 18:1;O2/16:0")
# → <rdkit.Chem.rdchem.Mol object>

conv.to_sdf(["PC 16:0/18:1(9Z)", "PE 18:0/20:4(5Z,8Z,11Z,14Z)"], path="output.sdf")
```

## Supported lipid classes

| Class | Examples | Status |
|-------|----------|--------|
| Fatty Acids (FA) | FA 16:0, FA 18:1(9Z), FA 20:4(5Z,8Z,11Z,14Z) | ✅ |
| Glycerolipids (GL) | MG, DG, TG | ✅ |
| Glycerophospholipids (GP) | PC, PE, PA, PI, PS, PG, LPC, LPE, ... | ✅ |
| Sphingolipids (SP) | Cer, SM, HexCer, Hex2Cer | ✅ |
| Sterols (ST) | Cholesterol, CE, bile acids | ✅ |
| Ether/plasmalogen linkages | O-, P- prefix | ✅ |

## Installation

### pip

```bash
pip install pylipidparse
```

### uv

```bash
uv add pylipidparse
```

### MCP server (for Claude / AI assistants)

```bash
pip install "pylipidparse[mcp]"   # requires Python 3.10+
```

### From source

```bash
git clone https://github.com/MontgomeryBohde/PyLipidParse.git
cd PyLipidParse
pip install -e ".[dev]"
```

## Quick start

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

# Write to file
conv.to_mol_file("PC 16:0/18:1(9Z)", "popc.mol")
conv.to_sdf(["FA 16:0", "FA 18:1(9Z)", "PC 16:0/18:1(9Z)"], "lipids.sdf")
```

## Notation requirements

PyLipidParse requires **full structural notation** with explicit chain positions.
Sum-composition notation (e.g., `PC 34:1`) is rejected because a unique structure
cannot be generated from it.

| Works | Fails |
|-------|-------|
| `PC 16:0/18:1(9Z)` | `PC 34:1` (no chain breakdown) |
| `TG 16:0/18:1(9Z)/18:2(9Z,12Z)` | `TG 16:0_18:1_18:2` (unknown positions) |
| `FA 18:1(9Z)` | `FA 18:1` (no double bond position) |

## Use with Claude (MCP server)

PyLipidParse includes an MCP server so you can convert lipid names directly inside Claude Code
or Claude Desktop — no Python required.

**Claude Code plugin** — install from [claude.com/plugins](https://claude.com/plugins) (search
*PyLipidParse*) or load locally:

```bash
claude --plugin-dir ./claude-code-plugin
```

**Manual MCP config** — add to `.mcp.json` or run:

```bash
claude mcp add pylipidparse -- uvx --from "pylipidparse[mcp]" pylipidparse-mcp
```

Then ask Claude naturally:

```
What is the SMILES for PC 16:0/18:1(9Z)?
Convert these ceramides to InChIKey: Cer 18:1;O2/16:0, Cer 18:1;O2/24:1(15Z)
Export these 20 lipids as an SDF file.
```

Six tools are available: `lipid_to_smiles`, `lipid_to_inchi`, `lipid_to_inchikey`,
`batch_convert_lipids`, `lipid_to_mol_file`, `lipids_to_sdf`.
See the [MCP documentation](https://montgomerybohde.github.io/PyLipidParse/mcp/) for full details.

## API reference

### `LipidConverter(dialect="LipidMaps", cache_size=512)`

| Method | Returns | Description |
|--------|---------|-------------|
| `to_mol(name)` | `Chem.Mol` | RDKit molecule (no 2D coords) |
| `to_smiles(name)` | `str` | Canonical SMILES |
| `to_inchi(name)` | `str` | InChI string |
| `to_inchikey(name)` | `str` | InChIKey (27-char hash) |
| `to_mol_file(name, path)` | — | Write `.mol` file with 2D coords |
| `to_sdf(names, path)` | — | Write `.sdf` file (supports batch) |

### Exceptions

```python
from pylipidparse.exceptions import (
    LipidParseError,                  # pygoslin couldn't parse the input
    UnsupportedLipidClassError,       # lipid class not yet implemented
    InsufficientStructuralDetailError, # species-level / unknown positions
    StructureGenerationError,         # molecule assembly failed
)
```

## Contributing

Issues and pull requests are welcome at [github.com/MontgomeryBohde/PyLipidParse](https://github.com/MontgomeryBohde/PyLipidParse/issues).

## Citation

If you use PyLipidParse in published work, please cite this repository and the
underlying tools:

- **pygoslin:** Kopczynski et al., *Analytical Chemistry* (2022). DOI: 10.1021/acs.analchem.1c05430
- **RDKit:** Landrum, G. *RDKit: Open-source cheminformatics.* [rdkit.org](https://www.rdkit.org)

## License

MIT License — see [LICENSE](LICENSE).
