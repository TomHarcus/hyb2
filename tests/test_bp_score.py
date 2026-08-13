"""Golden-diff parity: bp_score vs the legacy bp_score.sh.

bp_score reads a manifest (one row per folded structure -> 'ct offset len tail'),
extracts each structure's base pairs from its .ct into a .bps (offset-shifted
local->genome), then drives make_varna_scores. Inputs are the frozen constrained
fold (make_varna.ct) + experimental support (make_varna.basepair_scores.txt);
the outputs must reproduce the make_varna .bps and VARNA_scores goldens (bp_score
is exactly the driver that produced those).

bp_score writes files next to its inputs, so it runs in tmp_path with bare names
(the legacy assumes cwd + bare filenames).
"""

from hyb2.folding.bp_score import bp_score


def test_bp_score_matches_golden(folding_fixtures_dir, tmp_path, monkeypatch):
    (tmp_path / "frag.ct").write_text((folding_fixtures_dir / "make_varna.ct").read_text())
    (tmp_path / "s.basepair_scores.txt").write_text(
        (folding_fixtures_dir / "make_varna.basepair_scores.txt").read_text())
    # manifest: ct_name  offset(=begin-1)  len  tail
    (tmp_path / "coords.txt").write_text("frag.ct\t9999\t10300\t301\n")

    monkeypatch.chdir(tmp_path)
    with open("coords.txt") as c:
        bp_score(c, "s.basepair_scores.txt", None)

    assert (tmp_path / "frag.bps").read_text() == \
        (folding_fixtures_dir / "make_varna.bps").read_text()
    assert (tmp_path / "s__frag.VARNA_scores.txt").read_text() == \
        (folding_fixtures_dir / "make_varna.VARNA_scores.golden").read_text()
