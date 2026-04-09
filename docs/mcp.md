# MCP Server & Claude Plugin

PyLipidParse can run as an [MCP (Model Context Protocol)](https://modelcontextprotocol.io) server,
letting Claude Code, Claude Desktop, and any MCP-compatible AI assistant convert lipid names to
molecular structures directly in conversation — no Python code required.

## Quick start

### Option 1: One-line setup (no pre-install)

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

[uv](https://docs.astral.sh/uv/) must be installed. The first run downloads and caches
`pylipidparse[mcp]` automatically (~30–60 s due to RDKit).

### Option 2: Install then run

```bash
pip install "pylipidparse[mcp]"   # requires Python 3.10+
pylipidparse-mcp                  # starts the stdio MCP server
```

Then point Claude Code at the installed command:

```json
{
  "mcpServers": {
    "pylipidparse": {
      "command": "pylipidparse-mcp"
    }
  }
}
```

Or register it with the Claude Code CLI directly:

```bash
claude mcp add pylipidparse -- uvx --from "pylipidparse[mcp]" pylipidparse-mcp
```

## Claude Code plugin

The easiest way to use PyLipidParse in Claude Code is via the plugin, which wires up the MCP
server automatically.

**Install from the Claude Code plugins directory** at [claude.com/plugins](https://claude.com/plugins)
— search for *PyLipidParse*.

**Or load locally** from this repository:

```bash
git clone https://github.com/MontgomeryBohde/PyLipidParse.git
claude --plugin-dir ./PyLipidParse/claude-code-plugin
```

Once the plugin is active, you can ask Claude naturally:

```
What is the SMILES for PC 16:0/18:1(9Z)?
Give me the InChIKey for oleic acid (FA 18:1(9Z)).
Convert these 10 ceramides to SMILES and InChIKey: ...
Export these lipids as an SDF file.
```

## Available tools

| Tool | Description |
|------|-------------|
| `lipid_to_smiles` | Single lipid → canonical SMILES |
| `lipid_to_inchi` | Single lipid → InChI string |
| `lipid_to_inchikey` | Single lipid → InChIKey (27-char hash) |
| `batch_convert_lipids` | Multiple lipids → SMILES / InChI / InChIKey in one call |
| `lipid_to_mol_file` | Single lipid → MDL MOL file content (2D coordinates) |
| `lipids_to_sdf` | Multiple lipids → SDF file content (2D coordinates + metadata) |

## Input format

Lipid names must use **full structural notation** with explicit chain positions.
Species-level notation is rejected because it does not define a unique structure.

| Accepted | Rejected |
|----------|----------|
| `PC 16:0/18:1(9Z)` | `PC 34:1` |
| `FA 18:2(9Z,12Z)` | `FA 18:2` |
| `Cer 18:1;O2/16:0` | `Cer d18:1/16:0` |
| `TAG 16:0/18:1(9Z)/18:2(9Z,12Z)` | `TG 52:2` |
| `Cholesterol` | — |

## Supported dialects

Pass `dialect` to any tool to switch notation standard:

| Dialect | Standard |
|---------|----------|
| `LipidMaps` (default) | LIPID MAPS shorthand |
| `Goslin` | Goslin shorthand |
| `SwissLipids` | SwissLipids notation |
| `HMDB` | Human Metabolome Database |

## Example: batch conversion

```
Convert these lipids to SMILES and InChIKey:
PC 16:0/18:1(9Z), PE 18:0/20:4(5Z,8Z,11Z,14Z), Cer 18:1;O2/16:0
```

Claude will call `batch_convert_lipids` with `output_formats=["smiles","inchikey"]` and return
a table of results.

## Error responses

If a conversion fails, the tool returns a structured error object rather than raising an
exception, so batch jobs can continue past individual failures:

```json
{
  "lipid_name": "PC 34:1",
  "error": "insufficient_detail",
  "message": "...",
  "hint": "Use full structural notation such as 'PC 16:0/18:1(9Z)'."
}
```

| Error type | Cause |
|------------|-------|
| `parse_error` | Name not recognized by the selected dialect |
| `insufficient_detail` | Species-level notation — chain positions not specified |
| `unsupported_class` | Lipid class not yet supported by PyLipidParse |
| `structure_error` | Parsed successfully but molecule assembly failed |
