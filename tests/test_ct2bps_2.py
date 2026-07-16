"""Golden-diff parity: ct2bps_2 vs fixtures/folding_run/test.ua.cofold.bps.

Golden captured from bin/ct2bps_2.awk on the real RNAcofold .ct
(test.ua.cofold.ct, 10,428 folded chimeras -> 244,323 base-pair lines; see
scripts/generate_folding_baseline.sh). This is the .ct-parsing step of the
comradesMakeConstraints_2 / CoupleFold support-matrix path.

ct2bps_2() takes the .ct contents as a single string (it splitlines internally)
and yields "start<TAB>end\\n" base-pair rows.
"""

from hyb2.stages.ct2bps_2 import ct2bps_2


def test_ct2bps_2_matches_golden(folding_fixtures_dir):
    ct_text = (folding_fixtures_dir / "test.ua.cofold.ct").read_text()
    produced = "".join(ct2bps_2(ct_text))
    assert produced == (folding_fixtures_dir / "test.ua.cofold.bps").read_text()
