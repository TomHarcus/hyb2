"""Golden parity: stages/ct2b_gk3.py vs the legacy Ct2B_GK_3.pl (HYBRID_SS_MIN=1).

Input is a REAL hybrid-ss-min .ct (a folded 200-nt Zika fragment); the golden is
Ct2B_GK_3.pl's vienna output on it. Byte-exact -- ct2b_gk3 is fully deterministic
(no binary needed at test time, only the captured .ct + golden). Fixtures from
scripts/generate_unafold_baseline.sh.
"""

from hyb2.folding.ct2b_gk3 import ct2b_gk3


def test_matches_legacy(unafold_fixtures_dir):
    ct_text = (unafold_fixtures_dir / "ct2b.ct").read_text()
    produced = ct2b_gk3(ct_text, hybrid_ss_min=True)
    golden = (unafold_fixtures_dir / "ct2b_gk3.vienna.golden").read_text()
    assert produced == golden
