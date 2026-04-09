# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PyLipidParse converts lipid shorthand notation (LIPID MAPS, SwissLipids, HMDB) into SMILES, InChI, InChIKey, and RDKit Mol objects. It bridges pygoslin (which only parses notation) and RDKit (which handles cheminformatics) to generate actual molecular structures.

## Build & Development Commands

```bash
# Install in development mode (use the appropriate conda env)
conda run -n pylipidparse pip install -e ".[dev]"

# Run all tests (exclude benchmarks)
conda run -n pylipidparse pytest tests/ --ignore=tests/test_performance.py -v --override-ini="addopts="

# Run a single test file
conda run -n pylipidparse pytest tests/test_inchikey_validation.py -v --override-ini="addopts="

# Run a single test by name
conda run -n pylipidparse pytest tests/test_inchikey_validation.py -v -k "test_ceramide_inchikey" --override-ini="addopts="

# Lint and format checks (what CI runs)
ruff check src/ tests/
black --check src/ tests/

# Run benchmarks only
conda run -n pylipidparse pytest tests/test_performance.py --benchmark-only
```

Note: `--override-ini="addopts="` is needed because pyproject.toml sets `addopts` with `--cov` flags that require pytest-cov to be installed.

NEVER try to run code using the base Python installation. You must always run code using the correct conda env.

## Architecture

### Data Flow

```
User string → pygoslin.LipidParser.parse() → LipidMolecule object
  → LipidConverter._dispatch() → {FA,GL,GP,SP,ST}Builder.build()
  → RDKit Mol → SMILES/InChI/InChIKey/MOL/SDF
```

### Key Design Patterns

**SMILES string substitution for molecule assembly.** Scaffolds in `scaffolds/headgroups.py` use `{sn1}`, `{sn2}`, `{sn3}` placeholders. Chain fragments from `utils/chain.py` are substituted in, then the full SMILES is parsed by RDKit. No RWMol atom-by-atom assembly.

**Chain direction matters.** `build_acyl_chain()` in `utils/chain.py` has a `c1_first` parameter:
- `c1_first=False` (default): methyl-first order (Cn→C1), chain ends with `C(=O)`. Used for **sn-1** where scaffold has `{sn1}O` (chain LEFT of oxygen).
- `c1_first=True`: carboxyl-first order (C1→Cn), chain starts with `C(=O)`. Used for **sn-2** and **sn-3** where scaffold has `O{sn2}` (chain RIGHT of oxygen).
- Directional stereo bonds (`/`, `\`) are automatically flipped when `c1_first=True`.

**Builders extract pygoslin data via helpers in `builders/fatty_acid.py`:** `_extract_db_count()`, `_extract_db_positions()`, `_extract_modifications()` handle pygoslin's inconsistent API (e.g., `fa.double_bonds` can be `int` or `dict`).

### Module Roles

| Module | Role |
|---|---|
| `converter.py` | Public API (`LipidConverter`), pygoslin parsing, dispatch, caching |
| `builders/fatty_acid.py` | FA builder + pygoslin extraction helpers shared by all builders |
| `builders/glycerolipid.py` | GL builder (MAG/DAG/TAG); also exports `_build_chain_fragment()` reused by GP builder |
| `builders/glycerophospholipid.py` | GP builder (PC/PE/PA/PI/PS/PG + lyso); imports chain-building from GL |
| `builders/sphingolipid.py` | SP builder; builds sphingoid base SMILES with template placeholders `{N_ACYL}` and `{C1_HEAD}` |
| `builders/sterol.py` | ST builder; hardcoded ring SMILES for cholesterol, CE scaffold, bile acids |
| `scaffolds/headgroups.py` | All scaffold SMILES: glycerolipid, glycerophospholipid, and sphingolipid headgroups |
| `utils/chain.py` | Core chain SMILES generation: `build_acyl_chain()` and `build_alkyl_chain()` |

### Stereochemistry Notes

- Glycerol sn-2: always `[C@@H]` for natural (R)-configuration
- Sphingosine (2S,3R): in methyl-first SMILES, both C3 and C2 use `[C@@H]`
- PI inositol ring: C4 has no stereo annotation (plane of symmetry per PubChem)
- PS serine: `[C@H](N)` encodes L-serine (S) in the scaffold's SMILES context
- PG headgroup glycerol: no stereo (matches PubChem)

### Test Helpers (tests/conftest.py)

- `assert_smiles_equivalent()`: compare molecules by canonical SMILES, not string equality
- `assert_inchikey_match()`: decomposes failures into CONNECTIVITY / STEREO / CHARGE mismatch
- `assert_formula()`: check molecular formula via RDKit

### InChIKey Validation Tests

`tests/test_inchikey_validation.py` contains 49 tests validated against PubChem CIDs and LIPID MAPS. These are the gold-standard correctness tests. All reference InChIKeys are documented with their PubChem CID source.

## CI

GitHub Actions runs two jobs on push to main:
1. **Lint & type check**: `ruff check` + `black --check` (Python 3.11)
2. **Test matrix**: Python 3.10/3.11/3.12 via pip, plus Python 3.8/3.9/3.10 via conda with pinned RDKit versions

Line length limit is **100** (configured in pyproject.toml for both ruff and black).

## Git

You should make frequent commits to Github. Every time you complete a new standalone feature, you should make a new commit. 