"""Golden-diff parity: bp2hyb vs fixtures/folding_run/test.ranked_interactions.txt.

Golden captured from bin/bp2hyb.sh on the real fragment scores
(test.fragment_scores.txt: the ct2bps | histogram base-pair scores after the
comradesMakeConstraints_2 line-64 fragment window + printed<=1000 cap, 1001
lines -> 258 ranked-interaction lines; see scripts/generate_folding_baseline.sh).

The ~1000-line cap is load-bearing, not incidental: bp2hyb reformats every base
pair into an all-"RNA" chimera, so every record lands in one combine bucket and
combine runs O(n^2). The pipeline caps the input for exactly this reason.
"""

from hyb2.folding.bp2hyb import bp2hyb


def test_bp2hyb_matches_golden(folding_fixtures_dir):
    with (folding_fixtures_dir / "test.fragment_scores.txt").open() as f:
        produced = "".join(bp2hyb(f))
    assert produced == (folding_fixtures_dir / "test.ranked_interactions.txt").read_text()


def test_bp2hyb_trims_weak_helix_ends():
    from hyb2.folding.bp2hyb import bp2hyb
    # helix 10-50 .. 14-46 well supported, plus a 1-count pair 15-45 extending it, and a weak 5' end 9-51
    pairs = [(9, 51, 2), (10, 50, 100), (11, 49, 100), (12, 48, 90), (13, 47, 80), (14, 46, 100), (15, 45, 1)]
    lines = [f"{i}\t{j}\t{c}\n" for i, j, c in pairs]
    legacy = [l.split("\t") for l in bp2hyb(lines)]
    assert [(f[6], f[7], f[12], f[13]) for f in legacy] == [("9", "15", "45", "51")]
    trimmed = [l.split("\t") for l in bp2hyb(lines, trim_fraction=0.05)]
    assert [(f[6], f[7], f[12], f[13]) for f in trimmed] == [("10", "14", "46", "50")]
    assert float(trimmed[0][2]) == 470   # total support of the kept pairs
