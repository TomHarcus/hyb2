"""Golden-diff parity: make_varna_scores vs fixtures/folding_run/make_varna.*.

Golden captured from bin/make_VARNA_scores_2.sh on real inputs: the experimental
basepair_scores and a .bps from a constraint-honouring fold of the same fragment
(see scripts/generate_folding_baseline.sh step 7). 301 positions, 36 nonzero.

Params match how the golden was built: length=10300, tail=301, max_score=1000000.
Only the core transform is tested; the -s/-p filename derivation is CLI-only glue.
"""

from hyb2.folding.make_VARNA_scores_2 import make_varna_scores


def test_make_varna_scores_matches_golden(folding_fixtures_dir):
    with (folding_fixtures_dir / "make_varna.basepair_scores.txt").open() as bp, \
         (folding_fixtures_dir / "make_varna.bps").open() as bps:
        produced = make_varna_scores(bp, bps, length=10300, tail=301, max_score=1000000)
    expected = (folding_fixtures_dir / "make_varna.VARNA_scores.golden").read_text()
    assert "".join(produced) == expected
