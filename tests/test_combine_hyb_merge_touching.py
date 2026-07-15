"""Golden-diff parity: combine_hyb_merge_touching vs fixtures/tier2_run/.

Golden captured from bin/combine_hyb_merge_touching.pl on the full real
test.ua.hyb (10,428 records -> 10,425 output lines; see
scripts/generate_tier2_baseline.sh). Run with default switches, matching the
legacy script -- whose FOO=bar switch parsing is commented out, so TARGET /
GUIDE / TWO_WAY_MERGE / PRINT_SEQ_IDS are always at their defaults there. This
also transitively exercises the Hybrid2 class (Hybrid_long_2.pm port) end to
end, including the deliberately-reproduced join(",")/split(";") quirk.
"""

from hyb2.stages.combine_hyb_merge_touching import combine_hyb_merge_touching


def test_combine_hyb_merge_touching_matches_golden(tier2_fixtures_dir):
    with (tier2_fixtures_dir / "test.ua.hyb").open() as f:
        produced = "".join(combine_hyb_merge_touching(f))
    assert produced == (tier2_fixtures_dir / "test.combine.hyb").read_text()
