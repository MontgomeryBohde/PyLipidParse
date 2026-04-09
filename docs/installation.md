# Installation

## pip

```bash
pip install pylipidparse
```

## uv

```bash
uv add pylipidparse
```

## From source

```bash
git clone https://github.com/MontgomeryBohde/PyLipidParse.git
cd PyLipidParse
pip install -e ".[dev]"
```

## Dependencies

PyLipidParse requires:

- **pygoslin** (>=2.0.0) -- parses lipid shorthand notation into structured objects
- **RDKit** (>=2021.03) -- cheminformatics toolkit for molecule generation

Both are installed automatically via pip.
