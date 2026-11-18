"""Golden-diff parity: ct2bps_2 vs fixtures/folding_run/test.ua.cofold.bps.

Golden captured from bin/ct2bps_2.awk on the real RNAcofold .ct
(test.ua.cofold.ct, 10,428 folded chimeras -> 244,323 base-pair lines; see
scripts/generate_folding_baseline.sh). This is the .ct-parsing step of the
comradesMakeConstraints_2 / CoupleFold support-matrix path.

ct2bps_2() takes the .ct contents as a single string (it splitlines internally)
and yields "start<TAB>end\\n" base-pair rows.
"""

from hyb2.folding.ct2bps_2 import ct2bps_2


def test_ct2bps_2_matches_golden(folding_fixtures_dir):
    ct_text = (folding_fixtures_dir / "test.ua.cofold.ct").read_text()
    produced = "".join(ct2bps_2(ct_text, legacy=True))
    assert produced == (folding_fixtures_dir / "test.ua.cofold.bps").read_text()


def test_ct2bps_2_skips_intra_arm_and_orders_pairs():
    # arm 1 = ref 101-104 (GGGG), arm 2 = ref 11-14 (CCCC); local 1-8.
    # local 1-8 and 2-7 are inter-arm pairs; local 3-4 would be an intra-arm pair (fake, for the test)
    ct = (
        "8\tdG = -5.0\tx_RNA_101_104-x_RNA_11_14\n"
        "1 G 0 2 8 1\n2 G 1 3 7 2\n3 G 2 4 4 3\n4 G 3 5 3 4\n"
        "5 C 4 6 0 5\n6 C 5 7 0 6\n7 C 6 8 2 7\n8 C 7 0 1 8\n"
    )
    assert list(ct2bps_2(ct)) == ["14\t101\n", "13\t102\n"]
    legacy = list(ct2bps_2(ct, legacy=True))
    assert legacy[:2] == ["101\t14\n", "102\t13\n"] and len(legacy) == 4
