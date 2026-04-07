# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2024-01-01

### Added
- Initial release
- `LipidConverter` class with `to_mol`, `to_smiles`, `to_inchi`, `to_inchikey`, `to_mol_file`, `to_sdf` methods
- Support for Fatty Acids (FA)
- Support for Glycerolipids (MG, DG, TG)
- Support for Glycerophospholipids (PC, PE, PA, PI, PS, PG and lyso variants)
- Support for Sphingolipids (Cer, SM, HexCer, Hex2Cer)
- Support for Sterols (cholesterol, cholesterol esters)
- Ether and plasmalogen linkage support
- LRU caching for fast repeated conversions
- Python 3.8+ compatibility
- RDKit 2021.03+ compatibility
