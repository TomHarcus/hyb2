import subprocess

from hyb2.stages.solexa2fasta import solexa_to_fasta

def test_matches_legacy_awk_output():
    with open("data/synthetic/test.fastq") as f:
        fastq_text = f.read()

    awk_result = subprocess.run(
        ["awk", "-f", "bin/solexa2fasta.awk", "data/synthetic/test.fastq"],
        capture_output=True,
        text=True,
        check=True,
    )

    assert solexa_to_fasta(fastq_text) == awk_result.stdout