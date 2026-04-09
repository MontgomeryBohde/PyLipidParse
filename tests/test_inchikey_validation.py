"""InChIKey validation tests for all supported lipid classes.

Reference values sourced from two independent authoritative databases:
  - PubChem (NIH/NCBI): https://pubchem.ncbi.nlm.nih.gov/
  - LIPID MAPS Structure Database (LMSD): https://lipidmaps.org/databases/lmsd/

Where both sources were available, InChIKeys were confirmed to match exactly.
Compounds marked "(cross-verified)" were confirmed in both databases.

InChIKey anatomy (XXXXXXXXXXXXXX-XXXXXXXXXX-N):
  - Part 1 (14 chars): molecular connectivity (atom types + bonds)
  - Part 2 (10 chars): stereochemistry
  - Part 3 (1 char):   protonation/charge layer

Test failures are classified by the assert_inchikey_match helper:
  - CONNECTIVITY MISMATCH → structural bug in the builder (wrong atoms/bonds)
  - STEREO MISMATCH → correct connectivity but wrong stereocenters or double-bond geometry
  - CHARGE MISMATCH → connectivity and stereo correct but wrong protonation state
"""

import pytest

from pylipidparse import LipidConverter
from tests.conftest import assert_inchikey_match


@pytest.fixture(scope="module")
def conv():
    return LipidConverter()


# ---------------------------------------------------------------------------
# Fatty Acids
# ---------------------------------------------------------------------------


class TestFattyAcidInChIKeys:
    """InChIKey validation for free fatty acids.

    All references from PubChem; palmitic and oleic also cross-verified with
    LIPID MAPS (LMFA01010001 and LMFA01030002 respectively).
    """

    @pytest.mark.parametrize(
        "lipid_name,expected_inchikey,pubchem_cid",
        [
            # Acetic acid
            ("FA 2:0", "QTBSBXVTEAMEQO-UHFFFAOYSA-N", 176),
            # Butyric acid
            ("FA 4:0", "FERIUCNNQQJTOY-UHFFFAOYSA-N", 264),
            # Lauric acid
            ("FA 12:0", "POULHZVOKOAJMA-UHFFFAOYSA-N", 3893),
            # Palmitic acid (cross-verified: PubChem CID 985 + LMFA01010001)
            ("FA 16:0", "IPCSVZSSVZVIGE-UHFFFAOYSA-N", 985),
            # Stearic acid
            ("FA 18:0", "QIQXTHQIDYTFRH-UHFFFAOYSA-N", 5281),
        ],
    )
    def test_saturated_fa_inchikeys(self, conv, lipid_name, expected_inchikey, pubchem_cid):
        """Saturated fatty acid InChIKeys from PubChem."""
        ik = conv.to_inchikey(lipid_name)
        assert_inchikey_match(ik, expected_inchikey, f"{lipid_name} (PubChem CID {pubchem_cid})")

    @pytest.mark.parametrize(
        "lipid_name,expected_inchikey,pubchem_cid",
        [
            # Oleic acid (cross-verified: PubChem CID 445639 + LMFA01030002)
            ("FA 18:1(9Z)", "ZQPPMHVWECSIRJ-KTKRTIGZSA-N", 445639),
            # Linoleic acid
            ("FA 18:2(9Z,12Z)", "OYHQOLUKZRVURQ-HZJYTTRNSA-N", 5280450),
            # Alpha-linolenic acid
            ("FA 18:3(9Z,12Z,15Z)", "DTOSIQBPPRVQHS-PDBXOOCHSA-N", 5280934),
            # Arachidonic acid
            ("FA 20:4(5Z,8Z,11Z,14Z)", "YZXBAPSDXZZRGB-DOFZRALJSA-N", 444899),
            # DHA (docosahexaenoic acid)
            ("FA 22:6(4Z,7Z,10Z,13Z,16Z,19Z)", "MBMBGCFOFBJSGT-KUBAVDMBSA-N", 445580),
        ],
    )
    def test_unsaturated_fa_inchikeys(self, conv, lipid_name, expected_inchikey, pubchem_cid):
        """Unsaturated fatty acid InChIKeys from PubChem."""
        ik = conv.to_inchikey(lipid_name)
        assert_inchikey_match(ik, expected_inchikey, f"{lipid_name} (PubChem CID {pubchem_cid})")

    def test_z_vs_e_isomers_have_same_connectivity(self, conv):
        """Oleic (9Z) and elaidic (9E) share the same connectivity block.

        Both are octadecenoic acid — same atoms, same bonds, different geometry.
        Their InChIKeys should share the first 14 characters but differ in part 2.
        """
        ik_z = conv.to_inchikey("FA 18:1(9Z)")  # oleic
        ik_e = conv.to_inchikey("FA 18:1(9E)")  # elaidic
        assert ik_z is not None
        assert ik_e is not None
        # Same connectivity block
        assert ik_z.split("-")[0] == ik_e.split("-")[0] == "ZQPPMHVWECSIRJ"
        # Different stereo block
        assert (
            ik_z.split("-")[1] != ik_e.split("-")[1]
        ), "9Z and 9E isomers should have different stereo blocks"

    def test_elaidic_acid_inchikey(self, conv):
        """Elaidic acid (9E) from PubChem CID 637517."""
        ik = conv.to_inchikey("FA 18:1(9E)")
        assert_inchikey_match(ik, "ZQPPMHVWECSIRJ-MDZDMXLPSA-N", "FA 18:1(9E) (PubChem CID 637517)")


# ---------------------------------------------------------------------------
# Glycerolipids
# ---------------------------------------------------------------------------


class TestGlycerolipidInChIKeys:
    """InChIKey validation for glycerolipids (MG, DG, TG).

    All references from PubChem.
    """

    def test_mg_inchikey(self, conv):
        """1-Palmitoyl-sn-glycerol (MG 16:0) from PubChem CID 3084463."""
        ik = conv.to_inchikey("MG 16:0/0:0/0:0")
        assert_inchikey_match(
            ik, "QHZLMUACJMDIAE-SFHVURJKSA-N", "MG 16:0/0:0/0:0 (PubChem CID 3084463)"
        )

    def test_dg_inchikey(self, conv):
        """1-Palmitoyl-2-oleoyl-sn-glycerol (DG 16:0/18:1(9Z)) from PubChem CID 5282283."""
        ik = conv.to_inchikey("DG 16:0/18:1(9Z)/0:0")
        assert_inchikey_match(
            ik,
            "YEJYLHKQOBOSCP-OZKTZCCCSA-N",
            "DG 16:0/18:1(9Z)/0:0 (PubChem CID 5282283)",
        )

    def test_tg_saturated_inchikey(self, conv):
        """Tripalmitin (TG 16:0/16:0/16:0) from PubChem CID 11147."""
        ik = conv.to_inchikey("TG 16:0/16:0/16:0")
        assert_inchikey_match(
            ik, "PVNIQBQSYATKKL-UHFFFAOYSA-N", "TG 16:0/16:0/16:0 (PubChem CID 11147)"
        )

    def test_tg_mixed_inchikey(self, conv):
        """1-Palmitoyl-2-oleoyl-3-linoleoyl-glycerol from PubChem CID 9544086."""
        ik = conv.to_inchikey("TG 16:0/18:1(9Z)/18:2(9Z,12Z)")
        assert_inchikey_match(
            ik,
            "KGLAHZTWGPHKFF-FBSASISJSA-N",
            "TG 16:0/18:1(9Z)/18:2(9Z,12Z) (PubChem CID 9544086)",
        )


# ---------------------------------------------------------------------------
# Glycerophospholipids
# ---------------------------------------------------------------------------


class TestGlycerophospholipidInChIKeys:
    """InChIKey validation for glycerophospholipids (PC, PE, PA, PI, PS, PG + lyso).

    PC and PE cross-verified with LIPID MAPS (LMGP01010005 and LMGP02010009).
    All others from PubChem only.
    """

    def test_pc_inchikey(self, conv):
        """POPC (PC 16:0/18:1(9Z)) — cross-verified: PubChem CID 5497103 + LMGP01010005."""
        ik = conv.to_inchikey("PC 16:0/18:1(9Z)")
        assert_inchikey_match(
            ik, "WTJKGGKOPKCXLL-VYOBOKEXSA-N", "PC 16:0/18:1(9Z) (PubChem CID 5497103)"
        )

    def test_pe_inchikey(self, conv):
        """POPE (PE 16:0/18:1(9Z)) — cross-verified: PubChem CID 5283496 + LMGP02010009."""
        ik = conv.to_inchikey("PE 16:0/18:1(9Z)")
        assert_inchikey_match(
            ik, "FHQVHHIBKUMWTI-OTMQOFQLSA-N", "PE 16:0/18:1(9Z) (PubChem CID 5283496)"
        )

    def test_pa_inchikey(self, conv):
        """POPA (PA 16:0/18:1(9Z)) from PubChem CID 5283523."""
        ik = conv.to_inchikey("PA 16:0/18:1(9Z)")
        assert_inchikey_match(
            ik, "OPVZUEPSMJNLOM-QEJMHMKOSA-N", "PA 16:0/18:1(9Z) (PubChem CID 5283523)"
        )

    def test_pi_inchikey(self, conv):
        """POPI (PI 16:0/18:1(9Z)) from PubChem CID 71296232."""
        ik = conv.to_inchikey("PI 16:0/18:1(9Z)")
        assert_inchikey_match(
            ik, "PDLAMJKMOKWLAJ-ZNHRTHKOSA-N", "PI 16:0/18:1(9Z) (PubChem CID 71296232)"
        )

    def test_ps_inchikey(self, conv):
        """POPS (PS 16:0/18:1(9Z)) from PubChem CID 5283499."""
        ik = conv.to_inchikey("PS 16:0/18:1(9Z)")
        assert_inchikey_match(
            ik, "OIWCYIUQAVBPGV-DAQGAKHBSA-N", "PS 16:0/18:1(9Z) (PubChem CID 5283499)"
        )

    def test_pg_inchikey(self, conv):
        """POPG (PG 16:0/18:1(9Z)) from PubChem CID 5283509."""
        ik = conv.to_inchikey("PG 16:0/18:1(9Z)")
        assert_inchikey_match(
            ik, "PAZGBAOHGQRCBP-DDDNOICHSA-N", "PG 16:0/18:1(9Z) (PubChem CID 5283509)"
        )

    def test_lpc_inchikey(self, conv):
        """1-Palmitoyl-sn-glycero-3-phosphocholine (LPC 16:0) from PubChem CID 460602."""
        ik = conv.to_inchikey("LPC 16:0")
        assert_inchikey_match(ik, "ASWBNKHCZGQVJV-HSZRJFAPSA-N", "LPC 16:0 (PubChem CID 460602)")

    def test_lpe_inchikey(self, conv):
        """1-Palmitoyl-sn-glycero-3-phosphoethanolamine (LPE 16:0) from PubChem CID 9547069."""
        ik = conv.to_inchikey("LPE 16:0")
        assert_inchikey_match(ik, "YVYMBNSKXOXSKW-HXUWFJFHSA-N", "LPE 16:0 (PubChem CID 9547069)")

    def test_lpa_inchikey(self, conv):
        """1-Palmitoyl-sn-glycerol 3-phosphate (LPA 16:0) from PubChem CID 6419701."""
        ik = conv.to_inchikey("LPA 16:0")
        assert_inchikey_match(ik, "YNDYKPRNFWPPFU-GOSISDBHSA-N", "LPA 16:0 (PubChem CID 6419701)")


# ---------------------------------------------------------------------------
# Sphingolipids
# ---------------------------------------------------------------------------


class TestSphingolipidInChIKeys:
    """InChIKey validation for sphingolipids (Cer, SM, HexCer).

    Ceramide and sphingomyelin cross-verified with LIPID MAPS
    (LMSP02010004 and LMSP03010003 respectively).
    """

    def test_ceramide_inchikey(self, conv):
        """Cer d18:1/16:0 — cross-verified: PubChem CID 5283564 + LMSP02010004."""
        ik = conv.to_inchikey("Cer 18:1;O2/16:0")
        assert_inchikey_match(
            ik, "YDNKGFDKKRUKPY-TURZORIXSA-N", "Cer 18:1;O2/16:0 (PubChem CID 5283564)"
        )

    def test_sphingomyelin_inchikey(self, conv):
        """SM d18:1/16:0 — cross-verified: PubChem CID 9939941 + LMSP03010003."""
        ik = conv.to_inchikey("SM 18:1;O2/16:0")
        assert_inchikey_match(
            ik, "RWKUXQNLWDTSLO-GWQJGLRPSA-N", "SM 18:1;O2/16:0 (PubChem CID 9939941)"
        )

    def test_hexcer_inchikey(self, conv):
        """GlcCer d18:1/16:0 from PubChem CID 14035030 (beta-D-glucosyl-N-hexadecanoylsphingosine).

        Note: pygoslin normalises HexCer to 'Cer' with an [X] functional group; the builder
        resolves this back to HexCer. The reference InChIKey corresponds to the
        beta-glucosylceramide stereochemistry of the sugar attachment.
        """
        ik = conv.to_inchikey("HexCer 18:1;O2/16:0")
        assert_inchikey_match(
            ik, "VJLLLMIZEJJZTE-NNTBDIJYSA-N", "HexCer 18:1;O2/16:0 (PubChem CID 14035030)"
        )


# ---------------------------------------------------------------------------
# Sterols
# ---------------------------------------------------------------------------


class TestSterolInChIKeys:
    """InChIKey validation for sterols (cholesterol and cholesterol esters).

    Cholesterol cross-verified: PubChem CID 5997 + LIPID MAPS LMST01010001.
    All CE variants from PubChem only.
    """

    def test_cholesterol_st_inchikey(self, conv):
        """Cholesterol (ST 27:1;O) — cross-verified: PubChem CID 5997 + LMST01010001."""
        ik = conv.to_inchikey("ST 27:1;O")
        assert_inchikey_match(ik, "HVYWMOMLDIMFJA-DPAQBDIFSA-N", "ST 27:1;O (PubChem CID 5997)")

    def test_cholesterol_fc_inchikey(self, conv):
        """Free cholesterol (FC) should produce the same InChIKey as ST 27:1;O."""
        ik = conv.to_inchikey("FC")
        assert_inchikey_match(ik, "HVYWMOMLDIMFJA-DPAQBDIFSA-N", "FC (PubChem CID 5997)")

    def test_fc_and_st_produce_same_inchikey(self, conv):
        """FC and ST 27:1;O are both free cholesterol — must produce identical InChIKeys."""
        ik_fc = conv.to_inchikey("FC")
        ik_st = conv.to_inchikey("ST 27:1;O")
        assert ik_fc == ik_st, (
            f"FC and ST 27:1;O should map to the same structure.\n"
            f"  FC:       {ik_fc}\n"
            f"  ST 27:1;O:{ik_st}"
        )

    @pytest.mark.parametrize(
        "lipid_name,expected_inchikey,pubchem_cid",
        [
            # Cholesteryl palmitate
            ("CE 16:0", "BBJQPKLGPMQWBU-JADYGXMDSA-N", 246520),
            # Cholesteryl oleate
            ("CE 18:1(9Z)", "RJECHNNFRHZQMU-RMUVNZEASA-N", 5283632),
            # Cholesteryl stearate
            ("CE 18:0", "XHRPOTDGOASDJS-XNTGVSEISA-N", 118246),
            # Cholesteryl arachidonate
            ("CE 20:4(5Z,8Z,11Z,14Z)", "IMXSFYNMSOULQS-BEDFLICRSA-N", 6479222),
        ],
    )
    def test_cholesterol_ester_inchikeys(self, conv, lipid_name, expected_inchikey, pubchem_cid):
        """Cholesterol ester InChIKeys from PubChem."""
        ik = conv.to_inchikey(lipid_name)
        assert_inchikey_match(ik, expected_inchikey, f"{lipid_name} (PubChem CID {pubchem_cid})")


# ---------------------------------------------------------------------------
# Cross-class consistency
# ---------------------------------------------------------------------------


class TestCrossClassConsistency:
    """Sanity checks that hold across all lipid classes."""

    @pytest.mark.parametrize(
        "lipid_name",
        [
            "FA 16:0",
            "FA 18:1(9Z)",
            "MG 16:0/0:0/0:0",
            "DG 16:0/18:1(9Z)/0:0",
            "TG 16:0/16:0/16:0",
            "PC 16:0/18:1(9Z)",
            "PE 16:0/18:1(9Z)",
            "PA 16:0/18:1(9Z)",
            "PS 16:0/18:1(9Z)",
            "Cer 18:1;O2/16:0",
            "SM 18:1;O2/16:0",
            "ST 27:1;O",
            "CE 16:0",
        ],
    )
    def test_all_inchikeys_are_valid_format(self, conv, lipid_name):
        """Every supported lipid should produce a properly formatted InChIKey."""
        ik = conv.to_inchikey(lipid_name)
        assert ik is not None, f"InChIKey is None for {lipid_name}"
        assert len(ik) == 27, f"InChIKey wrong length for {lipid_name}: {ik!r}"
        assert ik[14] == "-", f"Missing hyphen at pos 14 for {lipid_name}: {ik!r}"
        assert ik[25] == "-", f"Missing hyphen at pos 25 for {lipid_name}: {ik!r}"

    def test_all_reference_lipids_have_unique_inchikeys(self, conv):
        """Each distinct lipid should produce a distinct InChIKey.

        Note: FA 18:1(9Z) and FA 18:1(9E) share connectivity but must differ.
        """
        lipids = [
            "FA 2:0",
            "FA 4:0",
            "FA 12:0",
            "FA 16:0",
            "FA 18:0",
            "FA 18:1(9Z)",
            "FA 18:1(9E)",
            "FA 18:2(9Z,12Z)",
            "FA 20:4(5Z,8Z,11Z,14Z)",
            "FA 22:6(4Z,7Z,10Z,13Z,16Z,19Z)",
            "MG 16:0/0:0/0:0",
            "DG 16:0/18:1(9Z)/0:0",
            "TG 16:0/16:0/16:0",
            "TG 16:0/18:1(9Z)/18:2(9Z,12Z)",
            "PC 16:0/18:1(9Z)",
            "PE 16:0/18:1(9Z)",
            "PA 16:0/18:1(9Z)",
            "PS 16:0/18:1(9Z)",
            "Cer 18:1;O2/16:0",
            "SM 18:1;O2/16:0",
            "ST 27:1;O",
            "CE 16:0",
            "CE 18:0",
        ]
        seen: dict = {}
        for name in lipids:
            ik = conv.to_inchikey(name)
            assert ik is not None, f"InChIKey is None for {name}"
            assert (
                ik not in seen
            ), f"InChIKey collision: {name!r} and {seen[ik]!r} both produce {ik}"
            seen[ik] = name
