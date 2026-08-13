"""Golden-diff parity: bp2hyb vs fixtures/folding_run/test.ranked_interactions.txt.

Golden captured from bin/bp2hyb.sh on the real fragment scores
(test.fragment_scores.txt: the ct2bps | histogram base-pair scores after the
comradesMakeConstraints_2 line-64 fragment window + printed<=1000 cap, 1001
lines -> 258 ranked-interaction lines; see scripts/generate_folding_baseline.sh).

The ~1000-line cap is load-bearing, not incidental: bp2hyb reformats every base
pair into an all-"RNA" chimera, so every record lands in one combine bucket and
combine runs O(n^2). The pipeline caps the input for exactly this reason.

bp2hyb() is a stdin-style filter: it takes an iterable of "i<TAB>j<TAB>count"
lines and yields the ranked, merged .hyb rows. Output is deterministic (both
GNU sorts fall back to a full-line tiebreak, reproduced in the port).
"""

from hyb2.folding.bp2hyb import bp2hyb


def test_bp2hyb_matches_golden(folding_fixtures_dir):
    with (folding_fixtures_dir / "test.fragment_scores.txt").open() as f:
        produced = "".join(bp2hyb(f))
    assert produced == (folding_fixtures_dir / "test.ranked_interactions.txt").read_text()
