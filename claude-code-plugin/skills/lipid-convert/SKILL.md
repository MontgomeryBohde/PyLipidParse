---
name: lipid-convert
description: >
  Convert lipid shorthand notation to molecular structures (SMILES, InChI, InChIKey,
  MOL files, SDF files). Use when the user asks about lipid structures, needs to look up
  a lipid by name, wants SMILES or InChIKey for database searches, or needs structure
  files for visualization or cheminformatics tools.
---

# Lipid Structure Conversion

Use the `pylipidparse` MCP tools to convert lipid shorthand notation to molecular structures.

## Available Tools

| Tool | Use for |
|------|---------|
| `lipid_to_smiles` | Single lipid → canonical SMILES string |
| `lipid_to_inchi` | Single lipid → InChI string |
| `lipid_to_inchikey` | Single lipid → InChIKey (27-char hash for database lookups) |
| `batch_convert_lipids` | Multiple lipids → SMILES/InChI/InChIKey in one call |
| `lipid_to_mol_file` | Single lipid → MDL MOL file content (2D coords) |
| `lipids_to_sdf` | Multiple lipids → SDF file content (2D coords + metadata) |

## Input Format Requirements

Lipid names must use **full structural notation** with explicit chain positions.

| Valid (full structural) | Invalid (species-level) |
|------------------------|------------------------|
| `PC 16:0/18:1(9Z)` | `PC 34:1` |
| `FA 18:2(9Z,12Z)` | `FA 18:2` |
| `Cer 18:1;O2/16:0` | `Cer 18:1/16:0` |
| `TAG 16:0/18:1(9Z)/18:2(9Z,12Z)` | `TG 52:2` |

Special cases that work without chain notation: `Cholesterol`, `FC`, `CHOL`.

## Supported Dialects

- `LipidMaps` (default) — LIPID MAPS shorthand notation
- `Goslin` — Goslin shorthand
- `SwissLipids` — SwissLipids notation
- `HMDB` — Human Metabolome Database notation

## Example Usage

**Single conversion:**
```
Convert PC 16:0/18:1(9Z) to SMILES
→ Use lipid_to_smiles with lipid_name="PC 16:0/18:1(9Z)"
```

**Batch conversion with multiple formats:**
```
Get SMILES and InChIKey for a list of lipids
→ Use batch_convert_lipids with output_formats=["smiles", "inchikey"]
```

**Structure file for visualization:**
```
Give me a MOL file for FA 18:1(9Z)
→ Use lipid_to_mol_file; the mol_content field contains the file text
```

**Export multiple lipids for a cheminformatics tool:**
```
Generate an SDF file for these 5 lipids
→ Use lipids_to_sdf; the sdf_content field contains the full SDF text
```

## Error Handling

If a conversion fails, the tool returns an `error` field instead of the structure.
Common errors:
- `parse_error` — name not recognized; check notation and dialect
- `insufficient_detail` — species-level notation; use full structural notation
- `unsupported_class` — lipid class not yet supported by PyLipidParse
- `structure_error` — parsed but could not generate a valid 3D structure
