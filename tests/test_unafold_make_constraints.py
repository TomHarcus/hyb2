"""Legacy parity: comrades_make_constraints fold="unafold" vs the legacy
bin/comradesMakeConstraints_2 -r 0 (hybrid-min), on the 20-chimera comrades_mini
window. Compares basepair_scores.

CONTENT-identical, not byte-identical: the .ct (from the same hybrid-min) and the
ct2bps_2 extraction match, but histogram tie-order differs between the Perl
histogram.pl (hash order) and the Python port (dict/insertion order) for
equal-count rows -- the documented "ordering fragility" family. So we diff the
sorted lines. Runs hybrid-min at test time, hence gated on unafold_env.

Golden + input from scripts/generate_unafold_baseline.sh.
"""

import shutil

from hyb2.folding.comrades_make_constraints import run


def test_basepair_scores_match_legacy(unafold_fixtures_dir, unafold_env, repo_root,
                                      tmp_path, monkeypatch):
    shutil.copy(unafold_fixtures_dir / "mini.hyb", tmp_path / "mini.hyb")
    shutil.copy(repo_root / "data" / "Zika_18S_formatted.fasta", tmp_path / "ref.fasta")
    monkeypatch.chdir(tmp_path)

    run("mini.hyb", "ref.fasta", 1, 10298, fold="unafold")

    produced = (tmp_path / "mini.basepair_scores.txt").read_text().splitlines()
    golden = (unafold_fixtures_dir / "mini.basepair_scores.golden").read_text().splitlines()
    assert sorted(produced) == sorted(golden)
