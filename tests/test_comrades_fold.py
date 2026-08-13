"""Regression golden for comrades_fold (port of comradesFold2), vienna backend.

comrades_fold implements the CORRECTED Option B greedy loop (adds ranked
base-pair constraints one at a time, keeps only those the fold honours). This
deliberately diverges from the legacy comradesFold2, whose constraint loop is a
silent no-op (see CLAUDE.md), so there is no legacy oracle to diff against --
the golden is a snapshot captured FROM THE PORT (Grzegorz-approved re-baseline).
The fold is deterministic (shuffling=False), so the .vienna/.ct are stable.

Inputs (constraints + fragment fasta) and outputs come from
scripts/generate_folding_baseline.sh step 7 (mini window 10000-10300). Needs
ViennaRNA (RNAfold/b2ct); skips otherwise.
"""

from hyb2.folding.comrades_fold import run


def test_comrades_fold_vienna_matches_snapshot(folding_fixtures_dir, vienna_bin,
                                               tmp_path, monkeypatch):
    g = folding_fixtures_dir
    (tmp_path / "c.txt").write_text((g / "comrades_fold.constraints").read_text())
    (tmp_path / "frag.fasta").write_text((g / "comrades_fold.frag.fasta").read_text())

    # comrades_fold writes .vienna/.ct/.aux next to in_fasta
    monkeypatch.chdir(tmp_path)
    vienna, ct = run("c.txt", "frag.fasta", fold="vienna", vienna_bin=vienna_bin)

    assert (tmp_path / vienna).read_text() == (g / "comrades_fold.vienna").read_text()
    assert (tmp_path / ct).read_text() == (g / "comrades_fold.ct").read_text()
