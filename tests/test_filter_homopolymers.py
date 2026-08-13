import subprocess

from hyb2.legacy_leaves.filter_homopolymers import filter_homopolymers


def _awk_reference(text: str) -> str:
    result = subprocess.run(
        ["awk", "-f", "bin/filter_homopolymers.awk"],
        input=text,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def test_matches_legacy_awk_output():
    with open("data/Zika_18S_formatted.fasta") as f:
        fasta_text = f.read()

    assert filter_homopolymers(fasta_text) == _awk_reference(fasta_text)


def test_line_under_threshold_is_kept():
    text = "ACGTACGTACGTACGT\n"
    assert filter_homopolymers(text) == _awk_reference(text)


def test_15_run_is_filtered_but_14_run_is_kept():
    text = "AAAAAAAAAAAAAA\nAAAAAAAAAAAAAAA\n"
    assert filter_homopolymers(text) == _awk_reference(text) == "AAAAAAAAAAAAAA\n"


def test_u_runs_are_not_filtered():
    text = "UUUUUUUUUUUUUUUU\n"
    assert filter_homopolymers(text) == _awk_reference(text) == "UUUUUUUUUUUUUUUU\n"


def test_mixed_case_run_is_not_filtered():
    text = "AaAaAaAaAaAaAaA\nAAAAAAAAAAAAAAa\n"
    assert filter_homopolymers(text) == _awk_reference(text)


def test_lowercase_run_is_filtered():
    text = "aaaaaaaaaaaaaaa\n"
    assert filter_homopolymers(text) == _awk_reference(text) == ""
