"""Golden/parity tests for the hyb2_fold orchestrator (mode 1, short-range).

Two kinds of coverage, per the folding tier's mixed provenance:

* transform + fragment extraction are legacy-faithful awk (bin/hyb2_fold lines
  73/76), so they are diffed byte-for-byte against PARITY goldens captured from
  the legacy awk.
* the fold is comrades_fold Option B (constraint-honouring, deliberately != the
  buggy legacy comradesFold2 no-op), so .vienna/.ct/log2scores are REGRESSION
  snapshots captured from the port itself.

Goldens come from scripts/generate_folding_baseline.sh step 10 (window
Zika_virusRNA 3900, len 150). The unit tests exercise the two awk-faithful
helpers directly (no ViennaRNA needed); the end-to-end test runs the whole
run() dispatch and needs ViennaRNA + the VARNA jar.
"""

from pathlib import Path

from hyb2.folding.hyb2_fold import run, _transform, _fasta_extraction

REPO = Path(__file__).resolve().parents[1]
REF = REPO / "data" / "Zika_18S_formatted.fasta"
GENE = "Zika_virusRNA"
X, L = 3900, 150
X2 = X + L - 1


def test_transform_matches_legacy(folding_fixtures_dir, tier2_fixtures_dir, tmp_path):
    out = tmp_path / "t.hyb"
    # mode-1 call shape from run(): GENE_2/y_coord/Y1/Y2/length all None
    _transform(str(tier2_fixtures_dir / "test.ua.hyb"), str(out),
               GENE, X, X, X2, None, None, None, None, None)
    assert out.read_text() == (folding_fixtures_dir / "hyb2_fold.mode1.transform.golden").read_text()


def test_fragment_matches_legacy(folding_fixtures_dir, tmp_path):
    out = tmp_path / "f.fasta"
    _fasta_extraction(str(REF), GENE, X, L, str(out))
    assert out.read_text() == (folding_fixtures_dir / "hyb2_fold.mode1.frag.golden").read_text()


def test_mode1_end_to_end(folding_fixtures_dir, tier2_fixtures_dir, vienna_bin,
                          varna_jar, tmp_path, monkeypatch):
    # run() copies the fasta into cwd and writes every intermediate there
    monkeypatch.chdir(tmp_path)
    (tmp_path / "hf.hyb").write_text((tier2_fixtures_dir / "test.ua.hyb").read_text())

    run("hf.hyb", GENE, None, str(REF), X, None, L, varna_jar, False, 1,
        vienna_bin=vienna_bin)

    g = folding_fixtures_dir
    # transform + fragment: parity vs legacy awk (proves run() wires the right helpers)
    assert (tmp_path / f"hf_{GENE}_{X}-{X2}.hyb").read_text() == (g / "hyb2_fold.mode1.transform.golden").read_text()
    assert (tmp_path / f"{GENE}_{X}-{X2}.fasta").read_text() == (g / "hyb2_fold.mode1.frag.golden").read_text()
    # fold outputs: Option B regression snapshots
    assert (tmp_path / f"{GENE}_{X}-{X2}.fasta.vienna").read_text() == (g / "hyb2_fold.mode1.vienna").read_text()
    assert (tmp_path / f"{GENE}_{X}-{X2}.fasta.ct").read_text() == (g / "hyb2_fold.mode1.ct").read_text()
    # log2scores: proves the %.6g formatting + the whole postfold chain
    assert (tmp_path / f"hf__{GENE}_{X}-{X2}.fasta.VARNA_log2scores.txt").read_text() == (g / "hyb2_fold.mode1.log2scores").read_text()
    # VARNA render produced the final plot
    assert (tmp_path / f"hf.{GENE}_{X}-{X2}.fasta_plot.svg").exists()
