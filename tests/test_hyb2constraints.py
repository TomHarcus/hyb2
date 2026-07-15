"""Golden-diff parity: hyb2constraints vs fixtures/tier2_run/test.constraints.

Golden captured from bin/hyb2constraints.pl on test.hyb (see
scripts/generate_tier2_baseline.sh). test.hyb deliberately carries a
'#'-comment header -- that header was the edge case that crashed an early
version of the port, so exercising it here is intentional.
"""

from hyb2.stages.hyb2constraints import hyb2constraints


def test_hyb2constraints_matches_golden(tier2_fixtures_dir):
    with (tier2_fixtures_dir / "test.hyb").open() as f:
        produced = "".join(hyb2constraints(f))
    assert produced == (tier2_fixtures_dir / "test.constraints").read_text()
