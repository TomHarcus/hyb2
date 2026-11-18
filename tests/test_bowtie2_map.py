"""FASTA input handling in bowtie2_map: raw reads are collapsed to <nr>_<count>
IDs before mapping (as fastq input is); already-collapsed input is passed through.
bowtie2 itself is stubbed out.
"""

import gzip

import pytest

from hyb2.mapping import bowtie2_map as bm
from hyb2.chimera.sam_composition import _check_read_ids


@pytest.fixture
def mapped(monkeypatch):
    """Stub out database building + bowtie2; record which reads file was mapped."""
    calls = {}
    monkeypatch.setattr(bm, "make_hyb_db_2", lambda db: None)
    monkeypatch.setattr(bm, "_bowtie2", lambda db, out, reads, reproducible, **kw: calls.update(reads=reads, **kw))
    return calls


RAW = ">@HWUSI-EAS721:8:1:4:358#0/1\nACGT\n>@HWUSI-EAS721:8:1:4:398#0/1\nGGCC\n>@HWUSI-EAS721:8:1:4:400#0/1\nACGT\n"


def test_raw_fasta_is_collapsed(tmp_path, mapped):
    reads = tmp_path / "reads.fasta"
    reads.write_text(RAW)

    bm.bowtie2_map(str(reads), "ref.fasta", str(tmp_path / "out"), adapter=None)

    assert mapped["reads"] == str(tmp_path / "out_comp.fasta")
    assert (tmp_path / "out_comp.fasta").read_text() == ">1_2\nACGT\n>2_1\nGGCC\n"


def test_raw_fasta_gz_is_collapsed(tmp_path, mapped):
    reads = tmp_path / "reads.fasta.gz"
    with gzip.open(reads, "wt") as f:
        f.write(RAW)

    bm.bowtie2_map(str(reads), "ref.fasta", str(tmp_path / "out"), adapter=None)

    assert (tmp_path / "out_comp.fasta").read_text() == ">1_2\nACGT\n>2_1\nGGCC\n"


def test_collapsed_fasta_is_mapped_as_is(tmp_path, mapped):
    reads = tmp_path / "reads.fasta"
    reads.write_text(">1_523\nACGT\n>2-3_7\nGGCC\n")

    bm.bowtie2_map(str(reads), "ref.fasta", str(tmp_path / "out"), adapter=None)

    assert mapped["reads"] == str(reads)
    assert not (tmp_path / "out_comp.fasta").exists()


@pytest.mark.parametrize("read_id", ["@HWUSI-EAS721:8:1:4:358#0/1", "HWUSI-EAS721:8:1:4:358#0/1"])
def test_sam_with_uncollapsed_ids_fails_clearly(tmp_path, read_id):
    sam = tmp_path / "x.sam"
    sam.write_text(f"@HD\tVN:1.5\n@SQ\tSN:g\tLN:10\n{read_id}\t0\tg\t1\t40\t4M\t*\t0\t0\tACGT\tIIII\n")

    with pytest.raises(ValueError, match="collapsed"):
        _check_read_ids(str(sam))


def test_sam_with_collapsed_ids_passes(tmp_path):
    sam = tmp_path / "x.sam"
    sam.write_text("@HD\tVN:1.5\n@SQ\tSN:g\tLN:10\n1_523\t0\tg\t1\t40\t4M\t*\t0\t0\tACGT\tIIII\n")

    _check_read_ids(str(sam))


ADAPTER = "TGGAATTCTCGGGTGCCAAGG"


@pytest.mark.skipif(__import__("shutil").which("cutadapt") is None, reason="cutadapt not installed")
def test_adapter_trimmed_before_collapsing(tmp_path, mapped):
    insert = "ACGTACGTTTGCAAGCTTGCA"            # 21 nt, kept (>= 18 after trimming)
    short = "ACGTACGTACG"                        # 11 nt, dropped by the length filter
    reads = tmp_path / "reads.fasta"
    reads.write_text(f">r1\n{insert}{ADAPTER[:12]}\n>r2\n{insert}{ADAPTER}AAAA\n>r3\n{short}{ADAPTER}\n")

    bm.bowtie2_map(str(reads), "ref.fasta", str(tmp_path / "out"))

    assert (tmp_path / "out_comp.fasta").read_text() == f">1_2\n{insert}\n"
    assert mapped["seed_length"] == 10 and mapped["score_min"] == "L,15,0"


def test_collapsed_input_is_not_trimmed(tmp_path, mapped, monkeypatch):
    monkeypatch.setattr(bm, "_trim_adapters", lambda *a: pytest.fail("collapsed reads must not be trimmed"))
    reads = tmp_path / "reads.fasta"
    reads.write_text(">1_2\nACGT\n")
    bm.bowtie2_map(str(reads), "ref.fasta", str(tmp_path / "out"))
    assert mapped["reads"] == str(reads)


def test_bowtie2_command_uses_settings(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(bm.subprocess, "run", lambda cmd, **kw: seen.setdefault("cmd", cmd))
    bm._bowtie2("ref.fasta", str(tmp_path / "o"), "r.fasta", False, seed_length=12,
                score_min="L,16,0", bowtie2_args="--norc -D 30")
    cmd = seen["cmd"]
    assert cmd[cmd.index("-L") + 1] == "12" and cmd[cmd.index("--score-min") + 1] == "L,16,0"
    assert cmd[-3:] == ["--norc", "-D", "30"]
