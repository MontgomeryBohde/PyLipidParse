"""Performance benchmarks for PyLipidParse.

Run with: pytest tests/test_performance.py --benchmark-only

Performance targets:
- Cold cache: < 100ms per lipid
- Hot cache: < 1ms per lipid (after first parse)
- Batch 100 unique: < 5 seconds
"""
import pytest

from pylipidparse import LipidConverter


@pytest.fixture
def conv():
    return LipidConverter()


@pytest.fixture
def warm_conv():
    """Converter with pre-populated cache."""
    c = LipidConverter()
    c.to_smiles("FA 18:1(9Z)")
    c.to_smiles("PC 16:0/18:1(9Z)")
    return c


@pytest.mark.benchmark(group="cold_cache")
def test_bench_fa_cold(benchmark, conv):
    benchmark(conv.to_smiles, "FA 18:1(9Z)")


@pytest.mark.benchmark(group="cold_cache")
def test_bench_pc_cold(benchmark, conv):
    benchmark(conv.to_smiles, "PC 16:0/18:1(9Z)")


@pytest.mark.benchmark(group="hot_cache")
def test_bench_fa_hot(benchmark, warm_conv):
    # Already in cache from fixture
    benchmark(warm_conv.to_smiles, "FA 18:1(9Z)")


@pytest.mark.benchmark(group="hot_cache")
def test_bench_pc_hot(benchmark, warm_conv):
    benchmark(warm_conv.to_smiles, "PC 16:0/18:1(9Z)")


@pytest.mark.benchmark(group="batch")
def test_bench_batch_fa_50(benchmark, conv):
    """50 unique fatty acids."""
    lipids = [f"FA {n}:0" for n in range(4, 54)]

    def run():
        for name in lipids:
            conv.to_smiles(name)

    benchmark(run)


@pytest.mark.benchmark(group="batch")
def test_bench_inchikey(benchmark, conv):
    benchmark(conv.to_inchikey, "FA 16:0")


@pytest.mark.benchmark(group="batch")
def test_bench_to_smiles_cached_100(benchmark):
    """100 calls to same lipid (all cached after first)."""
    c = LipidConverter()
    c.to_smiles("PC 16:0/18:1(9Z)")  # warm cache

    def run():
        for _ in range(100):
            c.to_smiles("PC 16:0/18:1(9Z)")

    benchmark(run)
