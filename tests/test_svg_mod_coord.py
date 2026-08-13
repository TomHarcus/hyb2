"""Golden-diff parity: svg_mod_coord vs the legacy svg_mod_coord.sh, both branches.

Golden captured from bin/svg_mod_coord.sh on a real VARNA-rendered SVG
(svg_mod.input.svg, from make_varna.ct; see scripts/generate_folding_baseline.sh
step 8). Branch 1 = single strand (-x 100); branch 3 = two strand
(-x 100 -y 5000 -l 150). The port writes <in>_plot.svg next to the input, so it
runs in tmp_path.
"""

from hyb2.folding.svg_mod_coord import svg_mod_coord


def _run(tmp_path, y_coord, length):
    svg = tmp_path / "frag.svg"
    svg.write_text((tmp_path / "_src.svg").read_text())
    svg_mod_coord(str(svg), 100, y_coord, length, None)
    return (tmp_path / "frag_plot.svg").read_text()


def test_svg_mod_coord_branch1(folding_fixtures_dir, tmp_path):
    (tmp_path / "_src.svg").write_text((folding_fixtures_dir / "svg_mod.input.svg").read_text())
    produced = _run(tmp_path, None, None)          # branch 1: Y unset
    assert produced == (folding_fixtures_dir / "svg_mod.branch1.golden").read_text()


def test_svg_mod_coord_branch3(folding_fixtures_dir, tmp_path):
    (tmp_path / "_src.svg").write_text((folding_fixtures_dir / "svg_mod.input.svg").read_text())
    produced = _run(tmp_path, 5000, 150)           # branch 3: both set
    assert produced == (folding_fixtures_dir / "svg_mod.branch3.golden").read_text()
