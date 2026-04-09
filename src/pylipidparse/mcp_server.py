"""MCP server for PyLipidParse: convert lipid shorthand notation to molecular structures.

Exposes the LipidConverter API as MCP tools over stdio transport.
Install with: pip install "pylipidparse[mcp]"  (requires Python >= 3.10)
Run with: pylipidparse-mcp
"""

from __future__ import annotations

import json
import tempfile

from mcp.server.fastmcp import FastMCP

from pylipidparse import (
    InsufficientStructuralDetailError,
    LipidConverter,
    LipidParseError,
    PyLipidParseError,
    StructureGenerationError,
    UnsupportedLipidClassError,
)

mcp = FastMCP(
    "pylipidparse",
    instructions=(
        "Convert lipid shorthand notation (LIPID MAPS, SwissLipids, HMDB, Goslin) "
        "to molecular structures: SMILES, InChI, InChIKey, MOL files, and SDF files."
    ),
)

# Global converter instance — lazy-initialized, one per dialect.
_converters: dict[str, LipidConverter] = {}

VALID_DIALECTS = ("LipidMaps", "Goslin", "SwissLipids", "HMDB")
VALID_FORMATS = ("smiles", "inchi", "inchikey")


def _get_converter(dialect: str) -> LipidConverter:
    """Return a cached LipidConverter for the given dialect."""
    if dialect not in VALID_DIALECTS:
        raise ValueError(f"Unknown dialect '{dialect}'. Valid options: {', '.join(VALID_DIALECTS)}")
    if dialect not in _converters:
        _converters[dialect] = LipidConverter(dialect=dialect)
    return _converters[dialect]


def _error_response(lipid_name: str, exc: Exception) -> dict:
    """Build a structured error dict from a PyLipidParse exception."""
    if isinstance(exc, LipidParseError):
        error_type = "parse_error"
        hint = "Check that the lipid name uses valid shorthand notation for the selected dialect."
    elif isinstance(exc, UnsupportedLipidClassError):
        error_type = "unsupported_class"
        hint = "This lipid class is recognized but not yet supported by PyLipidParse."
    elif isinstance(exc, InsufficientStructuralDetailError):
        error_type = "insufficient_detail"
        hint = (
            "Species-level notation (e.g. 'PC 34:1') does not specify individual chains. "
            "Use full structural notation such as 'PC 16:0/18:1(9Z)'."
        )
    elif isinstance(exc, StructureGenerationError):
        error_type = "structure_error"
        hint = "The parsed lipid could not be converted to a valid molecular structure."
    elif isinstance(exc, ValueError):
        error_type = "invalid_argument"
        hint = str(exc)
    else:
        error_type = "error"
        hint = str(exc)

    return {"lipid_name": lipid_name, "error": error_type, "message": str(exc), "hint": hint}


# ─────────────────────────────────────────────────────────────────────────────
# Tools
# ─────────────────────────────────────────────────────────────────────────────


@mcp.tool()
def lipid_to_smiles(lipid_name: str, dialect: str = "LipidMaps") -> str:
    """Convert a lipid shorthand name to a canonical SMILES string.

    Supports full structural notation where fatty acid chain positions are explicit.
    Species-level notation (e.g. 'PC 34:1') will fail because it is ambiguous.

    Args:
        lipid_name: Lipid in shorthand notation.
            Examples: 'PC 16:0/18:1(9Z)', 'FA 18:2(9Z,12Z)', 'Cer 18:1;O2/16:0',
            'TAG 16:0/18:1(9Z)/18:2(9Z,12Z)', 'Cholesterol'.
        dialect: Parsing dialect. Options: LipidMaps (default), Goslin, SwissLipids, HMDB.

    Returns:
        JSON object with 'lipid_name' and 'smiles', or an error object if conversion fails.
    """
    try:
        conv = _get_converter(dialect)
        smiles = conv.to_smiles(lipid_name)
        return json.dumps({"lipid_name": lipid_name, "smiles": smiles})
    except (PyLipidParseError, ValueError) as exc:
        return json.dumps(_error_response(lipid_name, exc))


@mcp.tool()
def lipid_to_inchi(lipid_name: str, dialect: str = "LipidMaps") -> str:
    """Convert a lipid shorthand name to an InChI string.

    Args:
        lipid_name: Lipid in shorthand notation.
            Examples: 'PC 16:0/18:1(9Z)', 'FA 18:2(9Z,12Z)', 'Cer 18:1;O2/16:0'.
        dialect: Parsing dialect. Options: LipidMaps (default), Goslin, SwissLipids, HMDB.

    Returns:
        JSON object with 'lipid_name' and 'inchi', or an error object if conversion fails.
    """
    try:
        conv = _get_converter(dialect)
        inchi = conv.to_inchi(lipid_name)
        return json.dumps({"lipid_name": lipid_name, "inchi": inchi})
    except (PyLipidParseError, ValueError) as exc:
        return json.dumps(_error_response(lipid_name, exc))


@mcp.tool()
def lipid_to_inchikey(lipid_name: str, dialect: str = "LipidMaps") -> str:
    """Convert a lipid shorthand name to an InChIKey (27-character hash).

    InChIKeys are useful for database lookups (PubChem, LIPID MAPS, ChEMBL).
    The first 14 characters encode connectivity; characters 15-25 encode stereochemistry.

    Args:
        lipid_name: Lipid in shorthand notation.
            Examples: 'PC 16:0/18:1(9Z)', 'FA 18:2(9Z,12Z)', 'Cer 18:1;O2/16:0'.
        dialect: Parsing dialect. Options: LipidMaps (default), Goslin, SwissLipids, HMDB.

    Returns:
        JSON object with 'lipid_name' and 'inchikey', or an error object if conversion fails.
    """
    try:
        conv = _get_converter(dialect)
        inchikey = conv.to_inchikey(lipid_name)
        return json.dumps({"lipid_name": lipid_name, "inchikey": inchikey})
    except (PyLipidParseError, ValueError) as exc:
        return json.dumps(_error_response(lipid_name, exc))


@mcp.tool()
def batch_convert_lipids(
    lipid_names: list[str],
    output_formats: list[str] | None = None,
    dialect: str = "LipidMaps",
) -> str:
    """Convert multiple lipid names to molecular identifiers in a single call.

    Errors for individual lipids are reported per-entry and do not abort the batch.

    Args:
        lipid_names: List of lipid names in shorthand notation.
            Example: ['PC 16:0/18:1(9Z)', 'FA 18:2(9Z,12Z)', 'Cer 18:1;O2/16:0']
        output_formats: Which identifiers to compute. Any combination of:
            'smiles', 'inchi', 'inchikey'. Defaults to ['smiles'].
        dialect: Parsing dialect. Options: LipidMaps (default), Goslin, SwissLipids, HMDB.

    Returns:
        JSON array of result objects, one per input lipid. Each object contains
        'lipid_name' plus the requested format fields, or an 'error' field on failure.
    """
    if output_formats is None:
        output_formats = ["smiles"]

    invalid_formats = [f for f in output_formats if f not in VALID_FORMATS]
    if invalid_formats:
        return json.dumps(
            {
                "error": "invalid_argument",
                "message": f"Unknown format(s): {invalid_formats}. Valid: {list(VALID_FORMATS)}",
            }
        )

    try:
        conv = _get_converter(dialect)
    except ValueError as exc:
        return json.dumps({"error": "invalid_argument", "message": str(exc)})

    results = []
    for name in lipid_names:
        entry: dict = {"lipid_name": name}
        try:
            if "smiles" in output_formats:
                entry["smiles"] = conv.to_smiles(name)
            if "inchi" in output_formats:
                entry["inchi"] = conv.to_inchi(name)
            if "inchikey" in output_formats:
                entry["inchikey"] = conv.to_inchikey(name)
        except (PyLipidParseError, ValueError) as exc:
            err = _error_response(name, exc)
            entry.update({k: v for k, v in err.items() if k != "lipid_name"})
        results.append(entry)

    return json.dumps(results)


@mcp.tool()
def lipid_to_mol_file(
    lipid_name: str,
    dialect: str = "LipidMaps",
    add_hydrogens: bool = False,
) -> str:
    """Convert a lipid shorthand name to MDL MOL file format (V2000).

    MOL files contain 2D atomic coordinates and bond tables, suitable for
    visualization tools (ChemDraw, MarvinJS, RDKit depiction) and cheminformatics
    workflows. The file content is returned as plain text.

    Args:
        lipid_name: Lipid in shorthand notation.
            Example: 'PC 16:0/18:1(9Z)', 'FA 18:1(9Z)'.
        dialect: Parsing dialect. Options: LipidMaps (default), Goslin, SwissLipids, HMDB.
        add_hydrogens: If True, add explicit hydrogen atoms to the structure.

    Returns:
        JSON object with 'lipid_name' and 'mol_content' (the MOL file as a string),
        or an error object if conversion fails.
    """
    try:
        conv = _get_converter(dialect)
        with tempfile.NamedTemporaryFile(suffix=".mol", mode="r", delete=False) as f:
            tmp_path = f.name
        conv.to_mol_file(lipid_name, tmp_path, add_hydrogens=add_hydrogens)
        with open(tmp_path, "r") as f:
            mol_content = f.read()
        return json.dumps({"lipid_name": lipid_name, "mol_content": mol_content})
    except (PyLipidParseError, ValueError) as exc:
        return json.dumps(_error_response(lipid_name, exc))
    except OSError as exc:
        return json.dumps(
            {"lipid_name": lipid_name, "error": "io_error", "message": str(exc), "hint": ""}
        )


@mcp.tool()
def lipids_to_sdf(
    lipid_names: list[str],
    dialect: str = "LipidMaps",
    add_hydrogens: bool = False,
) -> str:
    """Convert multiple lipids to SDF (Structure Data File) format.

    SDF files contain 2D atomic coordinates, bond tables, and embedded metadata
    (lipid name, SMILES, InChIKey) for each molecule. They are the standard exchange
    format for chemical structure databases and can be opened by any cheminformatics tool.

    Note: If any individual lipid fails, only the successful ones are included in the
    SDF. The 'failed' field lists any names that could not be converted.

    Args:
        lipid_names: List of lipid names in shorthand notation.
            Example: ['PC 16:0/18:1(9Z)', 'FA 18:2(9Z,12Z)']
        dialect: Parsing dialect. Options: LipidMaps (default), Goslin, SwissLipids, HMDB.
        add_hydrogens: If True, add explicit hydrogen atoms to all structures.

    Returns:
        JSON object with 'sdf_content' (the SDF file as a string) and 'failed' (list
        of lipid names that could not be converted), or an error object on total failure.
    """
    try:
        conv = _get_converter(dialect)
    except ValueError as exc:
        return json.dumps({"error": "invalid_argument", "message": str(exc)})

    # Filter out lipids that fail conversion, collecting failures separately.
    successful: list[str] = []
    failed: list[dict] = []
    for name in lipid_names:
        try:
            conv.to_smiles(name)  # Quick validation before writing SDF
            successful.append(name)
        except (PyLipidParseError, ValueError) as exc:
            failed.append(_error_response(name, exc))

    if not successful:
        return json.dumps(
            {
                "sdf_content": None,
                "failed": failed,
                "message": "No lipids could be converted successfully.",
            }
        )

    try:
        with tempfile.NamedTemporaryFile(suffix=".sdf", mode="r", delete=False) as f:
            tmp_path = f.name
        conv.to_sdf(successful, tmp_path, add_hydrogens=add_hydrogens)
        with open(tmp_path, "r") as f:
            sdf_content = f.read()
        return json.dumps({"sdf_content": sdf_content, "failed": failed})
    except OSError as exc:
        return json.dumps({"error": "io_error", "message": str(exc), "hint": ""})


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────


def main() -> None:
    """Start the PyLipidParse MCP server over stdio."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
