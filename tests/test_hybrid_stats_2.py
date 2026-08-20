import subprocess

from hyb2.legacy_leaves.hybrid_stats_2 import hybrid_stats_2


def _awk_reference(text: str) -> str:
    result = subprocess.run(
        ["awk", "-f", "legacy_bin/hybrid_stats_2"],
        input=text,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def test_matches_legacy_awk_output_on_real_fixture():
    with open("fixtures/legacy_run/main_short_range/test.hyb") as f:
        hyb_text = f.read()

    assert hybrid_stats_2(hyb_text) == _awk_reference(hyb_text)


def test_type_pairs_are_not_alphabetically_sorted():

    with open("fixtures/legacy_run/main_short_range/test.hyb") as f:
        result = hybrid_stats_2(f.read())

    assert "virusRNA:rRNA" in result
    assert "rRNA:virusRNA" in result


def test_sorted_descending_by_collapsed_count():
    with open("fixtures/legacy_run/main_short_range/test.hyb") as f:
        result = hybrid_stats_2(f.read())

    rows = result.splitlines()[1:]
    counts = [int(row.split("\t")[1]) for row in rows]
    assert counts == sorted(counts, reverse=True)


def test_header_printed_on_empty_input():
    assert hybrid_stats_2("") == "type1:type2\tcollapsed\tall\n"
