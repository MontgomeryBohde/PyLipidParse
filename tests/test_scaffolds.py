"""Unit tests for scaffold lookup functions in scaffolds/headgroups.py."""

import pytest

from pylipidparse.scaffolds.headgroups import (
    GLYCEROLIPID_SCAFFOLDS,
    get_glycerolipid_scaffold,
    get_gp_scaffold,
)


class TestGlycerolipidScaffolds:
    """Tests for get_glycerolipid_scaffold() edge cases."""

    def test_dag_12_scaffold(self):
        """DAG with sn-1 and sn-2 returns DAG_12 scaffold."""
        result = get_glycerolipid_scaffold(
            "DAG", sn1_present=True, sn2_present=True, sn3_present=False
        )
        assert result == GLYCEROLIPID_SCAFFOLDS["DAG_12"]

    def test_dag_13_scaffold(self):
        """DAG with sn-1 and sn-3 returns DAG_13 scaffold."""
        result = get_glycerolipid_scaffold(
            "DAG", sn1_present=True, sn2_present=False, sn3_present=True
        )
        assert result == GLYCEROLIPID_SCAFFOLDS["DAG_13"]

    def test_dag_23_scaffold(self):
        """DAG with sn-2 and sn-3 (sn-1 absent) returns DAG_23 scaffold."""
        result = get_glycerolipid_scaffold(
            "DAG", sn1_present=False, sn2_present=True, sn3_present=True
        )
        assert result == GLYCEROLIPID_SCAFFOLDS["DAG_23"]

    def test_mag_1_scaffold(self):
        """MAG at sn-1 returns MAG_1 scaffold."""
        result = get_glycerolipid_scaffold(
            "MAG", sn1_present=True, sn2_present=False, sn3_present=False
        )
        assert result == GLYCEROLIPID_SCAFFOLDS["MAG_1"]

    def test_mag_2_scaffold(self):
        """MAG at sn-2 returns MAG_2 scaffold."""
        result = get_glycerolipid_scaffold(
            "MAG", sn1_present=False, sn2_present=True, sn3_present=False
        )
        assert result == GLYCEROLIPID_SCAFFOLDS["MAG_2"]

    def test_mag_3_scaffold(self):
        """MAG at sn-3 returns MAG_3 scaffold."""
        result = get_glycerolipid_scaffold(
            "MAG", sn1_present=False, sn2_present=False, sn3_present=True
        )
        assert result == GLYCEROLIPID_SCAFFOLDS["MAG_3"]

    def test_tag_scaffold(self):
        """TAG with all sn positions returns TAG scaffold."""
        result = get_glycerolipid_scaffold(
            "TAG", sn1_present=True, sn2_present=True, sn3_present=True
        )
        assert result == GLYCEROLIPID_SCAFFOLDS["TAG"]

    def test_no_chains_raises(self):
        """No chains present raises ValueError."""
        with pytest.raises(ValueError, match="No chains present"):
            get_glycerolipid_scaffold(
                "MAG", sn1_present=False, sn2_present=False, sn3_present=False
            )


class TestGPScaffolds:
    """Tests for get_gp_scaffold() edge cases."""

    def test_full_pc_scaffold(self):
        """PC with both sn-1 and sn-2 returns PC scaffold."""
        result = get_gp_scaffold("PC", sn1_present=True, sn2_present=True)
        assert "{sn1}" in result
        assert "{sn2}" in result

    def test_lyso_pc_sn1_scaffold(self):
        """Lyso-PC at sn-1 returns LPC_SN1 scaffold."""
        result = get_gp_scaffold("PC", sn1_present=True, sn2_present=False)
        assert "{sn1}" in result

    def test_lyso_pc_sn2_scaffold(self):
        """Lyso-PC at sn-2 returns LPC_SN2 scaffold."""
        result = get_gp_scaffold("PC", sn1_present=False, sn2_present=True)
        assert "{sn2}" in result

    def test_nonexistent_headgroup_raises_key_error(self):
        """Non-existent headgroup and lyso combo raises KeyError."""
        # "NONEXISTENT" has no scaffold at all — both full and lyso lookups fail
        with pytest.raises(KeyError):
            get_gp_scaffold("NONEXISTENT", sn1_present=True, sn2_present=False)
