"""Golden-diff parity: hyb2fasta_bits_allRNAs vs fixtures/tier2_run/.

Goldens captured from bin/hyb2fasta_bits_allRNAs.awk on ref.tab (built from the
underscore-formatted Zika/18S fasta) + test.ua.hyb, producing the two paired
FASTA files (see scripts/generate_tier2_baseline.sh). The generator yields
(bit_1, bit_2) tuples, record N in one aligning to record N in the other, so
both halves are diffed against their respective goldens.
"""

from hyb2.stages.hyb2fasta_bits_allRNAs import hyb2fasta_bits_allRNAs


def test_hyb2fasta_bits_allRNAs_matches_golden(tier2_fixtures_dir):
    bit1, bit2 = [], []
    with (tier2_fixtures_dir / "ref.tab").open() as tab, \
         (tier2_fixtures_dir / "test.ua.hyb").open() as hyb:
        for bit1_rec, bit2_rec in hyb2fasta_bits_allRNAs(tab, hyb):
            bit1.append(bit1_rec)
            bit2.append(bit2_rec)
    assert "".join(bit1) == (tier2_fixtures_dir / "test.ua.bit_1.fasta").read_text()
    assert "".join(bit2) == (tier2_fixtures_dir / "test.ua.bit_2.fasta").read_text()
