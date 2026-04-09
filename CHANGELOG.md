# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2026-04-09

### Added
- MCP server (`pylipidparse-mcp`) exposing the full `LipidConverter` API as 6 MCP tools:
  `lipid_to_smiles`, `lipid_to_inchi`, `lipid_to_inchikey`, `batch_convert_lipids`,
  `lipid_to_mol_file`, `lipids_to_sdf`
- `[mcp]` optional dependency: `pip install "pylipidparse[mcp]"` (requires Python 3.10+)
- Claude Code plugin (`claude-code-plugin/`) for one-click MCP integration via `uvx`
- Smithery registry config (`smithery.yaml`)
- 28 new tests for the MCP server (319 total); MCP tests auto-skip on Python < 3.10

## [1.0.0] - 2026-04-09

### Added
- `LipidConverter` class with `to_mol`, `to_smiles`, `to_inchi`, `to_inchikey`, `to_mol_file`, `to_sdf` methods
- Fatty acid (FA) support: saturated, mono/polyunsaturated with positional double bonds
- Glycerolipid (GL) support: MG, DG, TG with sn-position assignment
- Glycerophospholipid (GP) support: PC, PE, PA, PI, PS, PG and all lyso variants
- Sphingolipid (SP) support: Cer, SM, HexCer, GlcCer, GalCer, Hex2Cer, Cer1P
- Sterol (ST) support: cholesterol, cholesterol esters (CE), bile acids (CA, DCA, CDCA, UDCA, LCA)
- Ether (O-) and plasmalogen (P-) linkage support
- LRU caching for fast repeated conversions (configurable size, default 512)
- 49 InChIKey validation tests verified against PubChem CIDs
- mkdocs-material documentation site with GitHub Pages deployment
- CI: lint/type-check + test matrix across Python 3.8-3.12 (pip and conda)
- PyPI publish workflow via trusted publishing (OIDC)
- Python 3.8+ and RDKit 2021.03+ compatibility
