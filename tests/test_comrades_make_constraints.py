"""Golden-diff parity: comrades_make_constraints pipeline vs the legacy
bin/comradesMakeConstraints_2, end-to-end (fasta2tab -> bit fastas -> RNAcofold
-> ct2bps -> histogram -> fragment cap -> bp2hyb -> hyb2constraints).

Golden captured from the legacy shell on a small 20-chimera input over a wide
window (comrades_mini.hyb, begin=1 end=10298 -> 75 folding constraints; see
scripts/generate_folding_baseline.sh step 6). The input is deliberately small
so the printed<=1000 fragment cap never truncates -- with no truncation the
histogram tie-order is irrelevant and the whole pipeline is deterministic, so
the byte-exact golden is stable. Requires ViennaRNA to fold.

The pipeline writes intermediates next to its inputs, so both inputs are copied
into tmp_path to keep the repo's data/ and fixtures/ clean.
"""

from pathlib import Path

from hyb2.folding.comrades_make_constraints import run


def test_comrades_make_constraints_matches_golden(
    folding_fixtures_dir, vienna_bin, tmp_path, repo_root
):
    hyb = tmp_path / "mini.hyb"
    hyb.write_text((folding_fixtures_dir / "comrades_mini.hyb").read_text())
    ref = tmp_path / "ref.fasta"
    ref.write_text((repo_root / "data" / "Zika_18S_formatted.fasta").read_text())

    constr = run(str(hyb), str(ref), 1, 10298, vienna_bin=vienna_bin, legacy=True)

    produced = Path(constr).read_text()
    golden = (folding_fixtures_dir / "comrades_mini.folding_constraints").read_text()
    assert produced == golden


def test_cofold_chunk_retried_after_crash(tmp_path, monkeypatch):
    import subprocess
    from hyb2.folding import comrades_make_constraints as m
    real = subprocess.run
    calls = {"n": 0}

    def flaky(cmd, **kw):
        if cmd[0].endswith("RNAcofold"):
            calls["n"] += 1
            if calls["n"] == 1:
                return subprocess.CompletedProcess(cmd, -11, "", "")
        return real(cmd, **kw)

    monkeypatch.setattr(m.subprocess, "run", flaky)
    (tmp_path / "b1.fa").write_text(">a_RNA_1_4\nGGGG\n")
    (tmp_path / "b2.fa").write_text(">b_RNA_11_14\nCCCC\n")
    m._fold_vienna(str(tmp_path / "b1.fa"), str(tmp_path / "b2.fa"), str(tmp_path / "x.ct"), None, workers=1)
    assert calls["n"] == 2 and (tmp_path / "x.ct").read_text().strip()
