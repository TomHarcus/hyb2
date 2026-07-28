"""Golden parity: plot_hybrids_3 (+ the swap prep) vs the legacy plot_hybrids_3.awk.

This is the deterministic contact-binning core of hyb2_coverage (before it shells
out to the R plotting scripts, which aren't golden-tested). Goldens captured from
bin/plot_hybrids_3.awk on the real test.ua.hyb -- see
scripts/generate_coverage_baseline.sh.

The awk emits its bins in hash order and the Python port in insertion order, so the
output *set* is the parity contract, not byte order -- compared as sorted lines.
"""

from hyb2.stages.plot_hybrids_3 import plot_hybrids_3, swap_gene1_to_arm1


def test_single_gene_matches_legacy(coverage_fixtures_dir, tier2_fixtures_dir):
    with open(tier2_fixtures_dir / "test.ua.hyb") as fh:
        out = plot_hybrids_3(fh, "Zika_virusRNA", "Zika_virusRNA", bin_size=10)
    golden = (coverage_fixtures_dir / "plot_hybrids_3.single.golden").read_text()
    assert sorted(out.splitlines()) == sorted(golden.splitlines())


def test_two_gene_with_swap_matches_legacy(coverage_fixtures_dir, tier2_fixtures_dir):
    with open(tier2_fixtures_dir / "test.ua.hyb") as fh:
        swapped = swap_gene1_to_arm1(fh, "Zika_virusRNA")
    out = plot_hybrids_3(swapped, "Zika_virusRNA", "18S_rRNA", bin_size=10)
    golden = (coverage_fixtures_dir / "plot_hybrids_3.twogene.golden").read_text()
    assert sorted(out.splitlines()) == sorted(golden.splitlines())
