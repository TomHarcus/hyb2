"""Golden-diff parity: plot_VARNA (headless) vs the legacy plot_VARNA.

Golden captured from bin/plot_VARNA on a real .ct + VARNA_scores (see
scripts/generate_folding_baseline.sh step 9). This test actually renders with
VARNA, so it needs java + the jar (skips otherwise). Runs in tmp_path because
plot_VARNA writes the SVG + <prefix>.<name>_plot.svg next to its inputs.

The interactive/GUI mode is not tested (nothing to diff); only the headless
.ct -> SVG -> svg_mod_coord path.
"""

from hyb2.pipelines.plot_VARNA import plot_VARNA


def test_plot_VARNA_headless_matches_golden(folding_fixtures_dir, varna_jar, tmp_path, monkeypatch):
    (tmp_path / "frag.ct").write_text((folding_fixtures_dir / "plot_varna.ct").read_text())
    # scores filename with '__' so the -p prefix (SCORE%%__*txt) resolves to "s"
    (tmp_path / "s__frag.VARNA_scores.txt").write_text(
        (folding_fixtures_dir / "plot_varna.scores.txt").read_text())

    # bare filenames (the legacy assumes cwd) -> output is s.frag_plot.svg
    monkeypatch.chdir(tmp_path)
    out = plot_VARNA(
        "frag.ct", "s__frag.VARNA_scores.txt", varna_jar,
        x_coord=100, y_coord=None, length=150, interactive=False,
    )

    produced = (tmp_path / out).read_text()
    assert produced == (folding_fixtures_dir / "plot_varna.golden").read_text()
