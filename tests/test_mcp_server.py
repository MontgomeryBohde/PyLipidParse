"""Tests for the PyLipidParse MCP server.

Three levels:
  - Unit: call tool functions directly, check JSON output shape and values
  - Error handling: verify each error type returns a structured error dict
  - Integration: use the MCP in-process test client to verify tool discovery
"""

import sys

import pytest

if sys.version_info < (3, 10):
    pytest.skip("MCP server requires Python 3.10+", allow_module_level=True)

import json  # noqa: E402

from pylipidparse.mcp_server import (
    batch_convert_lipids,
    lipid_to_inchi,
    lipid_to_inchikey,
    lipid_to_mol_file,
    lipid_to_smiles,
    lipids_to_sdf,
    mcp,
)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def parse(result: str) -> dict | list:
    return json.loads(result)


# ─────────────────────────────────────────────────────────────────────────────
# lipid_to_smiles
# ─────────────────────────────────────────────────────────────────────────────


class TestLipidToSmiles:
    def test_valid_fa(self):
        result = parse(lipid_to_smiles("FA 18:1(9Z)"))
        assert result["lipid_name"] == "FA 18:1(9Z)"
        assert "smiles" in result
        assert isinstance(result["smiles"], str)
        assert len(result["smiles"]) > 0

    def test_valid_pc(self):
        result = parse(lipid_to_smiles("PC 16:0/18:1(9Z)"))
        assert "smiles" in result
        assert "error" not in result

    def test_valid_cholesterol(self):
        result = parse(lipid_to_smiles("Cholesterol"))
        assert "smiles" in result
        assert "error" not in result

    def test_invalid_lipid_name(self):
        result = parse(lipid_to_smiles("NOTVALID XYZ"))
        assert "error" in result
        assert result["error"] == "parse_error"
        assert "hint" in result

    def test_species_level_notation(self):
        result = parse(lipid_to_smiles("PC 34:1"))
        assert "error" in result
        assert result["error"] == "insufficient_detail"

    def test_invalid_dialect(self):
        result = parse(lipid_to_smiles("FA 18:1(9Z)", dialect="BadDialect"))
        assert "error" in result

    def test_valid_cer(self):
        result = parse(lipid_to_smiles("Cer 18:1;O2/16:0"))
        assert "smiles" in result
        assert "error" not in result


# ─────────────────────────────────────────────────────────────────────────────
# lipid_to_inchi
# ─────────────────────────────────────────────────────────────────────────────


class TestLipidToInchi:
    def test_valid(self):
        result = parse(lipid_to_inchi("FA 18:1(9Z)"))
        assert "inchi" in result
        assert result["inchi"].startswith("InChI=")

    def test_invalid(self):
        result = parse(lipid_to_inchi("NOTVALID XYZ"))
        assert "error" in result


# ─────────────────────────────────────────────────────────────────────────────
# lipid_to_inchikey
# ─────────────────────────────────────────────────────────────────────────────


class TestLipidToInchikey:
    def test_valid(self):
        result = parse(lipid_to_inchikey("FA 18:1(9Z)"))
        assert "inchikey" in result
        key = result["inchikey"]
        assert isinstance(key, str)
        assert len(key) == 27
        assert key[14] == "-"

    def test_invalid(self):
        result = parse(lipid_to_inchikey("NOTVALID XYZ"))
        assert "error" in result


# ─────────────────────────────────────────────────────────────────────────────
# batch_convert_lipids
# ─────────────────────────────────────────────────────────────────────────────


class TestBatchConvertLipids:
    def test_single_smiles(self):
        result = parse(batch_convert_lipids(["FA 18:1(9Z)"]))
        assert isinstance(result, list)
        assert len(result) == 1
        assert "smiles" in result[0]

    def test_multiple_formats(self):
        result = parse(
            batch_convert_lipids(
                ["FA 18:1(9Z)", "PC 16:0/18:1(9Z)"],
                output_formats=["smiles", "inchi", "inchikey"],
            )
        )
        assert len(result) == 2
        for entry in result:
            assert "smiles" in entry
            assert "inchi" in entry
            assert "inchikey" in entry

    def test_partial_failure(self):
        result = parse(batch_convert_lipids(["FA 18:1(9Z)", "NOTVALID XYZ", "Cholesterol"]))
        assert len(result) == 3
        assert "smiles" in result[0]
        assert "error" in result[1]
        assert "smiles" in result[2]

    def test_invalid_format(self):
        result = parse(batch_convert_lipids(["FA 18:1(9Z)"], output_formats=["badformat"]))
        assert "error" in result
        assert result["error"] == "invalid_argument"

    def test_default_output_format_is_smiles(self):
        result = parse(batch_convert_lipids(["FA 18:1(9Z)"]))
        assert "smiles" in result[0]
        assert "inchi" not in result[0]

    def test_empty_list(self):
        result = parse(batch_convert_lipids([]))
        assert result == []


# ─────────────────────────────────────────────────────────────────────────────
# lipid_to_mol_file
# ─────────────────────────────────────────────────────────────────────────────


class TestLipidToMolFile:
    def test_valid(self):
        result = parse(lipid_to_mol_file("FA 18:1(9Z)"))
        assert "mol_content" in result
        content = result["mol_content"]
        assert isinstance(content, str)
        assert "M  END" in content  # All V2000 MOL files end with M  END

    def test_invalid(self):
        result = parse(lipid_to_mol_file("NOTVALID XYZ"))
        assert "error" in result

    def test_add_hydrogens(self):
        result_no_h = parse(lipid_to_mol_file("FA 18:1(9Z)", add_hydrogens=False))
        result_with_h = parse(lipid_to_mol_file("FA 18:1(9Z)", add_hydrogens=True))
        # With hydrogens the MOL file should be longer (more atom/bond lines)
        assert len(result_with_h["mol_content"]) > len(result_no_h["mol_content"])


# ─────────────────────────────────────────────────────────────────────────────
# lipids_to_sdf
# ─────────────────────────────────────────────────────────────────────────────


class TestLipidsToSdf:
    def test_valid_single(self):
        result = parse(lipids_to_sdf(["FA 18:1(9Z)"]))
        assert "sdf_content" in result
        assert result["failed"] == []
        assert "$$$$" in result["sdf_content"]  # SDF record terminator

    def test_valid_multiple(self):
        result = parse(lipids_to_sdf(["FA 18:1(9Z)", "PC 16:0/18:1(9Z)"]))
        assert result["sdf_content"].count("$$$$") == 2

    def test_partial_failure(self):
        result = parse(lipids_to_sdf(["FA 18:1(9Z)", "NOTVALID XYZ"]))
        assert result["sdf_content"].count("$$$$") == 1
        assert len(result["failed"]) == 1
        assert result["failed"][0]["lipid_name"] == "NOTVALID XYZ"

    def test_all_fail(self):
        result = parse(lipids_to_sdf(["NOTVALID1", "NOTVALID2"]))
        assert result["sdf_content"] is None
        assert len(result["failed"]) == 2

    def test_invalid_dialect(self):
        result = parse(lipids_to_sdf(["FA 18:1(9Z)"], dialect="BadDialect"))
        assert "error" in result


# ─────────────────────────────────────────────────────────────────────────────
# MCP integration: tool discovery via FastMCP server introspection
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.asyncio
async def test_mcp_tool_list():
    """Verify the MCP server exposes all 6 expected tools."""
    tools = await mcp.list_tools()
    tool_names = {t.name for t in tools}

    expected = {
        "lipid_to_smiles",
        "lipid_to_inchi",
        "lipid_to_inchikey",
        "batch_convert_lipids",
        "lipid_to_mol_file",
        "lipids_to_sdf",
    }
    assert expected == tool_names


@pytest.mark.asyncio
async def test_mcp_call_lipid_to_smiles():
    """End-to-end MCP server call_tool for lipid_to_smiles."""
    # call_tool returns (content_list, metadata_dict)
    result_content, _ = await mcp.call_tool("lipid_to_smiles", {"lipid_name": "FA 18:1(9Z)"})

    assert len(result_content) == 1
    result = json.loads(result_content[0].text)
    assert "smiles" in result
    assert "error" not in result


@pytest.mark.asyncio
async def test_mcp_tool_descriptions_not_empty():
    """Each tool must have a non-empty description (used by LLMs for tool selection)."""
    tools = await mcp.list_tools()

    for tool in tools:
        assert tool.description, f"Tool '{tool.name}' has no description"
        assert len(tool.description) > 20, f"Tool '{tool.name}' description is too short"
