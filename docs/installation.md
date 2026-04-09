# Installation

## pip (Python 3.10+)

```bash
pip install pylipidparse
```

## uv

```bash
uv add pylipidparse
```

## conda (Python 3.8+, recommended for older RDKit)

```bash
conda install -c conda-forge pylipidparse
```

## From source

```bash
git clone https://github.com/MontgomeryBohde/PyLipidParse.git
cd PyLipidParse
pip install -e ".[dev]"
```

!!! note "RDKit on older Python"
    For Python 3.8/3.9, pip wheels for RDKit may not be available.
    Use conda instead: `conda install -c conda-forge rdkit pygoslin`

## Dependencies

PyLipidParse requires:

- **pygoslin** (>=2.0.0) -- parses lipid shorthand notation into structured objects
- **RDKit** (>=2021.03) -- cheminformatics toolkit for molecule generation

Both are installed automatically via pip on Python 3.10+.
