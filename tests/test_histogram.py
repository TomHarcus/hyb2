"""Golden-diff parity: histogram vs fixtures/.../test.ua.hyb_stats_by_gene.txt.

The histogram input isn't captured on its own; it's pipeline step 7's
`cut -f4,10 test.ua.hyb`, derived here from the golden .ua.hyb. The golden
output is ordered by descending count, so this also pins histogram's ordering.
Remove xfail once implemented. (The single-read tophit histogram is covered by
the end-to-end pipeline test.)
"""

import pytest

from hyb2.stages.histogram import histogram


@pytest.mark.xfail(raises=NotImplementedError, strict=False, reason="Tier 1: stage not yet ported")
def test_histogram_gene_pair_matches_golden(fixtures_dir):
    ua = (fixtures_dir / "test.ua.hyb").read_text().splitlines()
    # pipeline step 7: cut -f4,10  (1-indexed cols 4 and 10 -> py indices 3, 9)
    cut = [f"{c[3]}\t{c[9]}\n" for c in (line.split("\t") for line in ua)]
    produced = "".join(histogram(cut))
    assert produced == (fixtures_dir / "test.ua.hyb_stats_by_gene.txt").read_text()
