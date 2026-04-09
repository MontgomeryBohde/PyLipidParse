# PyLipidParse Claude Code Plugin

This plugin exposes [PyLipidParse](https://github.com/MontgomeryBohde/PyLipidParse) as an MCP server inside Claude Code, letting you convert lipid shorthand notation to molecular structures directly in your AI assistant.

## What it does

Converts lipid names like `PC 16:0/18:1(9Z)` or `Cer 18:1;O2/16:0` into:
- Canonical **SMILES** strings
- **InChI** identifiers
- **InChIKey** hashes (for PubChem / LIPID MAPS database lookups)
- **MOL files** with 2D coordinates
- **SDF files** for multiple lipids at once

Supports LIPID MAPS, Goslin, SwissLipids, and HMDB notation dialects.

## Requirements

- [uv](https://docs.astral.sh/uv/) must be installed (used to run the server without pre-installation)

## Installation

Install this plugin from the Claude Code plugins directory, or manually:

```bash
claude --plugin-dir /path/to/PyLipidParse/claude-code-plugin
```

The MCP server runs via `uvx` and installs `pylipidparse[mcp]` automatically on first use. Subsequent runs use the cached environment.

## Manual MCP configuration (without the plugin)

Add to your `.mcp.json` or Claude Code MCP config:

```json
{
  "mcpServers": {
    "pylipidparse": {
      "command": "uvx",
      "args": ["--from", "pylipidparse[mcp]", "pylipidparse-mcp"]
    }
  }
}
```

Or if you have `pylipidparse[mcp]` already installed:

```json
{
  "mcpServers": {
    "pylipidparse": {
      "command": "pylipidparse-mcp"
    }
  }
}
```

## Example prompts

- *"Convert PC 16:0/18:1(9Z) to SMILES"*
- *"What is the InChIKey for oleic acid (FA 18:1(9Z))?"*
- *"Give me SMILES and InChIKey for these 10 lipids: ..."*
- *"Generate an SDF file for this list of ceramides"*

## Supported lipid classes

Fatty acids (FA), glycerolipids (MAG/DAG/TAG), glycerophospholipids (PC/PE/PA/PI/PS/PG and lyso variants), sphingolipids (Cer/SM/HexCer/LacCer/SPB), and sterols (Cholesterol/CE/bile acids).

## Notes

- **Full structural notation required** — use `PC 16:0/18:1(9Z)`, not `PC 34:1`
- First run may take 30–60 seconds while `uvx` installs dependencies
- Source: [github.com/MontgomeryBohde/PyLipidParse](https://github.com/MontgomeryBohde/PyLipidParse)
