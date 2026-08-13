import subprocess

from hyb2.legacy_leaves.txt2hyb import add_sequences_to_hyb

TAB_FILE = "READ1\tACGTACGT\nREAD2\tTTTTGGGG\n"
HYB_FILE = (
    "READ1\tXXXXXXXX\tcol3\tcol4\n"
    "READ3\tYYYYYYYY\tcol3\tcol4\n"
    "READ2\tZZZZZZZZ\tcol3\tcol4\n"
)


def _awk_reference(tab_file: str, hyb_file: str, tmp_path) -> str:
    tab_path = tmp_path / "tab.txt"
    hyb_path = tmp_path / "hyb.txt"
    tab_path.write_text(tab_file)
    hyb_path.write_text(hyb_file)

    result = subprocess.run(
        ["awk", "-f", "bin/txt2hyb.awk", str(tab_path), str(hyb_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def test_matches_legacy_awk_output(tmp_path):
    assert add_sequences_to_hyb(TAB_FILE, HYB_FILE) == _awk_reference(
        TAB_FILE, HYB_FILE, tmp_path
    )


def test_substitutes_second_column_and_keeps_trailing_tab():
    result = add_sequences_to_hyb(TAB_FILE, HYB_FILE)
    assert result == (
        "READ1\tACGTACGT\tcol3\tcol4\t\n"
        "READ2\tTTTTGGGG\tcol3\tcol4\t\n"
    )


def test_unmatched_read_id_is_dropped():
    result = add_sequences_to_hyb(TAB_FILE, HYB_FILE)
    assert "READ3" not in result
