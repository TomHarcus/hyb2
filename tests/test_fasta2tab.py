import subprocess

from hyb2.common.fasta2tab import fasta_to_tab


def test_matches_legacy_awk_output():
    with open("data/Zika_18S_formatted.fasta") as f:
        fasta_text = f.read()

    awk_result = subprocess.run(
        ["awk", "-f", "bin/fasta2tab.awk", "data/Zika_18S_formatted.fasta"],
        capture_output=True,
        text=True,
        check=True,
    )

    assert fasta_to_tab(fasta_text) == awk_result.stdout


def test_multiline_sequence_is_concatenated():
    fasta = ">seq1\nACGT\nACGT\n>seq2\nTTTT\n"
    assert fasta_to_tab(fasta) == "seq1\tACGTACGT\nseq2\tTTTT\n"


def test_comment_and_blank_lines_are_dropped():
    fasta = "# a comment\n>seq1\n\nACGT\n# another comment\nACGT\n"
    assert fasta_to_tab(fasta) == "seq1\tACGTACGT\n"


def test_header_with_spaces_is_kept_verbatim():
    fasta = ">seq1 some description\nACGT\n"
    assert fasta_to_tab(fasta) == "seq1 some description\tACGT\n"
