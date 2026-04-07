"""Headgroup scaffold SMILES templates for glycerolipids and glycerophospholipids.

Conventions
-----------
- ``{sn1}`` is replaced by the sn-1 chain fragment (methyl-first, ends at C1)
- ``{sn2}`` is replaced by the sn-2 chain fragment
- ``{sn3}`` is replaced by the sn-3 chain fragment

For **ester** linkages, chain fragments end with ``C(=O)``:
    The scaffold has ``-O{sn}`` which becomes ``-OC(=O)[chain]`` = ester bond.

For **ether** (O-) linkages, chain fragments end with ``C``:
    The scaffold has ``-O{sn}`` which becomes ``-OC[chain]`` = ether bond.

For **plasmalogen** (P-) linkages, chain fragments start with ``/C=C\\``:
    The scaffold has ``-O{sn}`` which becomes ``-O/C=C\\[chain]`` = vinyl ether.

Stereochemistry
---------------
The glycerol sn-2 carbon is ``[C@@H]`` for the natural (R)-configuration
(sn-glycero-3-phospho... configuration in glycerophospholipids).

IMPORTANT: All scaffold SMILES have been written based on known reference
structures. They should be validated against PubChem or LIPID MAPS references
before production use. Run the scaffold validation script in Appendix B of
ImplementationPlan.md.

References
----------
PC 16:0/18:1 — PubChem CID: 5497103
PE 16:0/18:1 — PubChem CID: 5497103 (verify)
Cholesterol   — PubChem CID: 5997
"""
from typing import Dict

# ============================================================================
# Glycerolipid scaffolds (MAG, DAG, TAG)
# ============================================================================
# Glycerol backbone: sn-1 (primary CH2), sn-2 (secondary CH, chiral), sn-3 (primary CH2)
# Natural sn-2 configuration: (R) → [C@@H] in these scaffolds
#
# TAG: all three positions occupied (triacylglycerol)
# DAG: two positions occupied (diacylglycerol)
# MAG: one position occupied (monoacylglycerol)

GLYCEROLIPID_SCAFFOLDS: Dict[str, str] = {
    # TAG: sn-1/sn-2/sn-3 all esterified
    "TAG": "{sn1}OC[C@@H](O{sn2})CO{sn3}",
    # DAG: sn-1 and sn-2 esterified, sn-3 = free OH
    "DAG_12": "{sn1}OC[C@@H](O{sn2})CO",
    # DAG: sn-1 and sn-3 esterified, sn-2 = free OH
    "DAG_13": "{sn1}OC[C@@H](O)CO{sn3}",
    # DAG: sn-2 and sn-3 esterified, sn-1 = free OH
    "DAG_23": "OC[C@@H](O{sn2})CO{sn3}",
    # MAG: sn-1 only
    "MAG_1": "{sn1}OCC(O)CO",
    # MAG: sn-2 only
    "MAG_2": "OC[C@@H](O{sn2})CO",
    # MAG: sn-3 only
    "MAG_3": "OCC(O)CO{sn3}",
}

# ============================================================================
# Glycerophospholipid scaffolds
# ============================================================================
# The phosphate group is at sn-3 (primary carbon).
# Each headgroup differs only in what is attached to the phosphate oxygen.
#
# Template format: {sn1}OC[C@@H](O{sn2})COP(=O)([O-])O[headgroup]
#
# Ionization state: phosphate is shown as P(=O)([O-]) (deprotonated, pH ~7).
# The choline nitrogen in PC is shown as [N+](C)(C)C (quaternary ammonium).

GLYCEROPHOSPHOLIPID_SCAFFOLDS: Dict[str, str] = {
    # PC: phosphatidylcholine
    "PC": "{sn1}OC[C@@H](O{sn2})COP(=O)([O-])OCC[N+](C)(C)C",
    # PE: phosphatidylethanolamine
    "PE": "{sn1}OC[C@@H](O{sn2})COP(=O)([O-])OCCN",
    # PA: phosphatidic acid
    "PA": "{sn1}OC[C@@H](O{sn2})COP(=O)([O-])[O-]",
    # PI: phosphatidylinositol (inositol ring with natural D-myo configuration)
    "PI": (
        "{sn1}OC[C@@H](O{sn2})COP(=O)([O-])"
        "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    ),
    # PS: phosphatidylserine
    "PS": "{sn1}OC[C@@H](O{sn2})COP(=O)([O-])OC[C@@H](N)C(=O)O",
    # PG: phosphatidylglycerol
    "PG": "{sn1}OC[C@@H](O{sn2})COP(=O)([O-])OC[C@@H](O)CO",
    # Lyso variants: one chain position = free OH
    # LPC sn-1 occupied, sn-2 = OH
    "LPC_SN1": "{sn1}OC[C@@H](O)COP(=O)([O-])OCC[N+](C)(C)C",
    # LPC sn-2 occupied, sn-1 = OH
    "LPC_SN2": "OC[C@@H](O{sn2})COP(=O)([O-])OCC[N+](C)(C)C",
    # LPE sn-1 occupied, sn-2 = OH
    "LPE_SN1": "{sn1}OC[C@@H](O)COP(=O)([O-])OCCN",
    # LPE sn-2 occupied, sn-1 = OH
    "LPE_SN2": "OC[C@@H](O{sn2})COP(=O)([O-])OCCN",
    # LPA sn-1
    "LPA_SN1": "{sn1}OC[C@@H](O)COP(=O)([O-])[O-]",
    "LPA_SN2": "OC[C@@H](O{sn2})COP(=O)([O-])[O-]",
    # LPI sn-1
    "LPI_SN1": (
        "{sn1}OC[C@@H](O)COP(=O)([O-])"
        "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    ),
    "LPI_SN2": (
        "OC[C@@H](O{sn2})COP(=O)([O-])"
        "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    ),
    # LPS sn-1
    "LPS_SN1": "{sn1}OC[C@@H](O)COP(=O)([O-])OC[C@@H](N)C(=O)O",
    "LPS_SN2": "OC[C@@H](O{sn2})COP(=O)([O-])OC[C@@H](N)C(=O)O",
    # LPG sn-1
    "LPG_SN1": "{sn1}OC[C@@H](O)COP(=O)([O-])OC[C@@H](O)CO",
    "LPG_SN2": "OC[C@@H](O{sn2})COP(=O)([O-])OC[C@@H](O)CO",
}

# ============================================================================
# Sphingolipid headgroup modifications at C1-OH of sphingoid base
# ============================================================================
# These are attached to the primary alcohol (C1-OH) of the sphingoid base.
# The template uses {ceramide_c1} for the C1-OH oxygen.

SPHINGOLIPID_HEADGROUPS: Dict[str, str] = {
    # Ceramide: no modification on C1-OH (just the free alcohol)
    "Cer": "O",  # C1-OH stays as OH
    # Sphingomyelin: phosphocholine on C1-OH
    "SM": "OP(=O)([O-])OCC[N+](C)(C)C",
    # Hexosylceramide (generic hexose)
    "HexCer": "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O",
    # Glucosylceramide (Glc = beta-D-glucose)
    "GlcCer": "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O",
    # Galactosylceramide (Gal = beta-D-galactose)
    "GalCer": "O[C@@H]1O[C@H](CO)[C@H](O)[C@@H](O)[C@@H]1O",
    # Lactosylceramide (Hex2Cer): Glc-Gal
    "Hex2Cer": (
        "O[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@H](O)[C@@H](O)[C@@H]2O)"
        "[C@@H](O)[C@H](O)[C@@H]1O"
    ),
    # Ceramide-1-phosphate
    "Cer1P": "OP(=O)([O-])[O-]",
}


def get_glycerolipid_scaffold(
    lipid_class: str, sn1_present: bool, sn2_present: bool, sn3_present: bool
) -> str:
    """Select the appropriate glycerolipid scaffold."""
    n_chains = sum([sn1_present, sn2_present, sn3_present])
    if n_chains == 3:
        return GLYCEROLIPID_SCAFFOLDS["TAG"]
    elif n_chains == 2:
        if sn1_present and sn2_present:
            return GLYCEROLIPID_SCAFFOLDS["DAG_12"]
        elif sn1_present and sn3_present:
            return GLYCEROLIPID_SCAFFOLDS["DAG_13"]
        else:
            return GLYCEROLIPID_SCAFFOLDS["DAG_23"]
    elif n_chains == 1:
        if sn1_present:
            return GLYCEROLIPID_SCAFFOLDS["MAG_1"]
        elif sn2_present:
            return GLYCEROLIPID_SCAFFOLDS["MAG_2"]
        else:
            return GLYCEROLIPID_SCAFFOLDS["MAG_3"]
    raise ValueError("No chains present")


def get_gp_scaffold(headgroup: str, sn1_present: bool, sn2_present: bool) -> str:
    """Select the appropriate glycerophospholipid scaffold."""
    if sn1_present and sn2_present:
        return GLYCEROPHOSPHOLIPID_SCAFFOLDS[headgroup]

    # Lyso variants
    lyso_key = f"L{headgroup}_SN{'1' if sn1_present else '2'}"
    if lyso_key in GLYCEROPHOSPHOLIPID_SCAFFOLDS:
        return GLYCEROPHOSPHOLIPID_SCAFFOLDS[lyso_key]

    raise KeyError(
        f"No scaffold found for headgroup {headgroup!r} with "
        f"sn1={sn1_present}, sn2={sn2_present}"
    )
