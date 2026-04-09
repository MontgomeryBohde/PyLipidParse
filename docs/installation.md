# Installation

## pip

```bash
pip install pylipidparse
```

## uv

```bash
uv add pylipidparse
```

## MCP server (for Claude / AI assistants)

To use PyLipidParse as an MCP server with Claude Code or Claude Desktop, install the `[mcp]`
extra. This requires **Python 3.10+**.

```bash
pip install "pylipidparse[mcp]"
```

Or run it directly without pre-installing via [uvx](https://docs.astral.sh/uv/):

```bash
uvx --from "pylipidparse[mcp]" pylipidparse-mcp
```

See the [MCP Server & Claude Plugin](mcp.md) page for full setup instructions.

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

Both are installed automatically via pip. The optional `[mcp]` extra additionally requires
**mcp** (>=1.0.0) and **Python 3.10+**.
