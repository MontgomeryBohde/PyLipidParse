# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PyLipidParse converts lipid shorthand notation (LIPID MAPS, SwissLipids, HMDB, Goslin) into SMILES, InChI, InChIKey, and RDKit Mol objects. It bridges pygoslin (which only parses notation) and RDKit (which handles cheminformatics) to generate actual molecular structures.

Published on PyPI as `pylipidparse`. Docs at https://montgomerybohde.github.io/PyLipidParse.

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

# Build and serve docs locally
conda run -n pylipidparse pip install ".[docs]"
conda run -n pylipidparse mkdocs serve        # dev server at http://localhost:8000
conda run -n pylipidparse mkdocs build --strict  # production build (outputs to site/)

# Run MCP server (requires Python 3.10+)
conda run -n pylipidparse pip install ".[mcp]"
conda run -n pylipidparse pylipidparse-mcp
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

**Abstract builder base class** in `builders/base.py`: all builders inherit from `AbstractLipidBuilder`, which provides `_mol_from_smiles()` and `_sanitize()`.

**RDKit version compatibility** is isolated in `_compat.py`: `mol_to_inchi()`, `mol_to_inchikey()`, and `compute_2d_coords()` abstract over API differences across RDKit 2021.03+.

### Module Roles

| Module | Role |
|---|---|
| `converter.py` | Public API (`LipidConverter`), pygoslin parsing, dispatch, caching |
| `builders/base.py` | `AbstractLipidBuilder` base class with `_mol_from_smiles()` and `_sanitize()` |
| `builders/fatty_acid.py` | FA builder + pygoslin extraction helpers shared by all builders |
| `builders/glycerolipid.py` | GL builder (MAG/DAG/TAG); also exports `_build_chain_fragment()` reused by GP builder |
| `builders/glycerophospholipid.py` | GP builder (PC/PE/PA/PI/PS/PG + lyso); imports chain-building from GL |
| `builders/sphingolipid.py` | SP builder; builds sphingoid base SMILES with template placeholders `{N_ACYL}` and `{C1_HEAD}` |
| `builders/sterol.py` | ST builder; hardcoded ring SMILES for cholesterol, CE scaffold, bile acids |
| `scaffolds/headgroups.py` | All scaffold SMILES: glycerolipid, glycerophospholipid, and sphingolipid headgroups |
| `utils/chain.py` | Core chain SMILES generation: `build_acyl_chain()` and `build_alkyl_chain()` |
| `_compat.py` | RDKit version compatibility shims (InChI, InChIKey, 2D coords) |
| `mcp_server.py` | MCP server exposing converter tools over stdio (Python 3.10+) |

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

## MCP Server

`mcp_server.py` exposes PyLipidParse as an MCP (Model Context Protocol) server over stdio, allowing Claude and other MCP-compatible tools to convert lipid names without writing Python.

**6 tools:** `lipid_to_smiles`, `lipid_to_inchi`, `lipid_to_inchikey`, `batch_convert_lipids`, `lipid_to_mol_file`, `lipids_to_sdf`.

**Dialects:** LipidMaps (default), Goslin, SwissLipids, HMDB.

**Error handling:** Tools return structured error dicts (never raise). Error types: `parse_error`, `insufficient_detail`, `unsupported_class`, `structure_error`, `invalid_argument`, `io_error`.

**Entry point:** `pylipidparse-mcp` (defined in pyproject.toml scripts).

**Claude Code plugin** in `claude-code-plugin/` — wires up the MCP server via uvx so no pre-installation is needed. Can be installed from the Claude plugins directory or locally with `claude mcp add`.

**MCP test files** require Python 3.10+ and `pytest-asyncio`.

## Documentation

Built with **MkDocs + Material theme** + `mkdocstrings[python]` for auto-generated API docs.

- Config: `mkdocs.yml`
- Pages: `docs/` — index, installation, usage, mcp, supported_classes, api_reference
- API reference is auto-generated from docstrings in `src/pylipidparse`
- Deployed to GitHub Pages on push to main via `.github/workflows/docs.yml`

## CI/CD

Three GitHub Actions workflows in `.github/workflows/`:

### `ci.yml` — Lint + Tests (on push to main & PRs)

1. **`lint`**: ruff + black (Python 3.11)
2. **`test-pip`**: Python 3.10/3.11/3.12 with latest RDKit via pip. Uploads coverage to Codecov.
3. **`test-conda`**: Older Python + pinned RDKit versions via conda-forge:
   - Python 3.8 + RDKit 2021.09 (no MCP tests)
   - Python 3.9 + RDKit 2022.09 (no MCP tests)
   - Python 3.10 + RDKit 2023.03 (includes MCP + pytest-asyncio)

### `docs.yml` — Deploy Documentation (on push to main)

Builds with `mkdocs build --strict` and deploys to GitHub Pages.

### `publish.yml` — Publish to PyPI (on `v*` tags)

Uses trusted publishing (OIDC) with hatchling build backend. No API token needed.

Line length limit is **100** (configured in pyproject.toml for both ruff and black).

## Git

You should make frequent commits to Github. Every time you complete a new standalone feature, you should make a new commit.
