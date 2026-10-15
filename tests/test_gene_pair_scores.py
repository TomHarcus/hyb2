from pathlib import Path

import pytest

from hyb2.folding.gene_pair_scores import paired_coordinates, score_gene_pair
from hyb2.folding.comrades_make_constraints import main, run


def ct_block(token, intra=False):
    partners = (2, 1, 0, 0) if intra else (4, 3, 2, 1)
    return f"4 dG = -1 {token}a-{token}b\n" + "".join(
        f"{i} G {i-1} {i+1} {partner} {i}\n"
        for i, partner in enumerate(partners, 1)
    )


def test_coordinate_orientation_and_intra_arm_exclusion():
    records = {
        "fa-fb": (10, 11, 20, 21, False),
        "ra-rb": (20, 21, 10, 11, True),
        "ia-ib": (10, 11, 20, 21, False),
    }
    rows = list(paired_coordinates(ct_block("f") + ct_block("r") + ct_block("i", True), records))
    assert rows == ["10\t21\n", "11\t20\n", "11\t20\n", "10\t21\n"]


def hyb_row(read, first, second, a=10, b=20, width=2):
    return "\t".join(map(str, [read, ".", ".", first, 1, width, a, a+width-1,
                                ".", second, width+1, width*2, b, b+width-1, "."])) + "\n"


def test_mixed_genes_and_delimiters(tmp_path, monkeypatch):
    from hyb2.folding import comrades_make_constraints as pipeline
    ref = tmp_path / "ref.fa"
    ref.write_text("".join(f">{g}\n{'G'*40}\n" for g in ("gene_A-one", "gene_B-two", "other")))
    hyb = tmp_path / "reads.hyb"
    hyb.write_text(
        hyb_row("read", "gene_A-one", "gene_B-two")
        + hyb_row("many_parts_read_5", "gene_B-two", "gene_A-one", 20, 10)
        + hyb_row("irrelevant", "other", "gene_B-two")
    )

    def fake_fold(first, second, ct, vienna_bin):
        headers = Path(first).read_text().splitlines()[::2]
        assert headers == [">hyb2pair1a", ">hyb2pair2a"]
        Path(ct).write_text(ct_block("hyb2pair1") + ct_block("hyb2pair2"))

    monkeypatch.setattr(pipeline, "_fold_vienna", fake_fold)
    output = run(str(hyb), str(ref), gene1="gene_A-one", gene2="gene_B-two")
    assert Path(output).read_text() == "gene_A-one\tgene_B-two\tscore\n10\t21\t2\n11\t20\t2\n"
    assert not list(tmp_path.glob("*folding_constraints*"))


@pytest.mark.parametrize("reuse", [False, True])
@pytest.mark.parametrize("rows", [
    hyb_row("read", "A", "B"),
    hyb_row("read1", "A", "A") + hyb_row("read2", "B", "B"),
])
def test_existing_routine_rejects_mixed_genes_before_writing(tmp_path, rows, reuse):
    hyb = tmp_path / "mixed.hyb"
    hyb.write_text(rows)
    with pytest.raises(ValueError, match="Mixed-gene HYB"):
        run(str(hyb), str(tmp_path / "unused.fa"), 1, 40,
            support_matrix="unused.txt" if reuse else None)
    assert list(tmp_path.iterdir()) == [hyb]


@pytest.mark.parametrize("reverse", [False, True])
def test_two_rna_fold_preprocessing_remains_compatible(tmp_path, vienna_bin, reverse):
    from hyb2.folding.hyb2_fold import _transform, _fasta_extraction_two_region
    ref = tmp_path / "ref.fa"
    ref.write_text(">A\nGGGGGGGGGG\n>B\nCCCCCCCCCC\n")
    hyb = tmp_path / "raw.hyb"
    hyb.write_text(hyb_row("read_1", "B" if reverse else "A",
                           "A" if reverse else "B", 1, 1, 10))
    transformed = tmp_path / "selected.hyb"
    fragment = tmp_path / "fragment.fa"
    _transform(str(hyb), str(transformed), "A", 1, 1, 10, "B", 1, 1, 10, 10)
    _fasta_extraction_two_region(str(ref), "A", 1, 1, 10, str(fragment), "B")
    constraints = run(str(transformed), str(fragment), 1, 120, vienna_bin=vienna_bin)
    assert Path(constraints).exists()
    scores = (tmp_path / "selected.basepair_scores.txt").read_text().splitlines()
    assert scores
    for row in scores:
        x, y, count = map(int, row.split())
        if reverse:
            x, y = y, x
        assert 1 <= x <= 10 and 111 <= y <= 120 and count == 1


def test_empty_pair_and_self_pair(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">A\nGGGG\n>B\nCCCC\n")
    hyb = tmp_path / "empty.hyb"
    hyb.write_text("")
    output = score_gene_pair(str(hyb), str(ref), "A", "A")
    assert Path(output).read_text() == "A\tA\tscore\n"
    with pytest.raises(ValueError, match="Gene not found"):
        score_gene_pair(str(hyb), str(ref), "missing", "B")


@pytest.mark.parametrize("options", [["-1", "A"], ["-1", "A", "-2", "B", "-b", "1"], []])
def test_cli_rejects_ambiguous_modes(options):
    with pytest.raises(SystemExit):
        main(["-i", "unused.hyb", "-f", "unused.fa", *options])


def test_top_level_scores_command(tmp_path, monkeypatch):
    from hyb2.cli import main as cli_main
    monkeypatch.chdir(tmp_path)
    Path("ref.fa").write_text(">A\nGGGG\n>B\nCCCC\n")
    Path("reads.hyb").write_text("")
    assert cli_main(["scores", "-i", "reads.hyb", "-f", "ref.fa",
                     "-1", "A", "-2", "B", "-o", "pair.tsv"]) == 0
    assert Path("pair.tsv").read_text() == "A\tB\tscore\n"


@pytest.mark.parametrize("backend", ["vienna", "unafold"])
def test_real_folding_preserves_opaque_ids(tmp_path, monkeypatch, backend):
    import shutil
    if not shutil.which("RNAcofold" if backend == "vienna" else "hybrid-min"):
        pytest.skip("folding tool unavailable")
    ref = tmp_path / "ref.fa"
    ref.write_text(">A_one-two\nGGGGGGGGGG\n>B_three-four\nCCCCCCCCCC\n")
    hyb = tmp_path / "reads.hyb"
    hyb.write_text(hyb_row("arbitrary_read_name", "A_one-two", "B_three-four", 1, 1, 10))
    result = score_gene_pair(str(hyb), str(ref), "A_one-two", "B_three-four", fold=backend)
    rows = Path(result).read_text().splitlines()
    assert rows[0] == "A_one-two\tB_three-four\tscore"
    assert len(rows) > 1
    assert all(1 <= int(x) <= 10 and 1 <= int(y) <= 10 and int(n) == 1
               for x, y, n in (row.split("\t") for row in rows[1:]))
