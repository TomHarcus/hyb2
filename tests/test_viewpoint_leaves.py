"""Golden parity for the two deterministic leaves of plot_viewpoint:
hyb2blast (vs hyb2blast.awk) and blast2gplot (vs blast2gplot.pl).

Both are fully deterministic, so these are byte-exact diffs. viewpoint_graph.R
(the actual plotting) isn't golden-tested, same as the VARNA/CDM rendering.
Goldens from scripts/generate_viewpoint_baseline.sh.
"""

from hyb2.viewpoint.hyb2blast import hyb2blast
from hyb2.viewpoint.blast2gplot import blast2gplot


def test_hyb2blast_matches_legacy(viewpoint_fixtures_dir, tier2_fixtures_dir):
    with open(tier2_fixtures_dir / "test.ua.hyb") as f:
        produced = hyb2blast(f)
    assert produced == (viewpoint_fixtures_dir / "hyb2blast.golden").read_text()


def test_blast2gplot_matches_legacy(viewpoint_fixtures_dir, tmp_path, monkeypatch):
    g = viewpoint_fixtures_dir
    (tmp_path / "ref.blast").write_text((g / "blast2gplot.ref.blast").read_text())
    (tmp_path / "in.blast").write_text((g / "blast2gplot.blast").read_text())
    (tmp_path / "lengths.txt").write_text((g / "blast2gplot.lengths.txt").read_text())

    # blast2gplot writes <exp>_<gene>.gplot into cwd
    monkeypatch.chdir(tmp_path)
    blast2gplot(exp="exp", n_genes=1, ref_blast_file="ref.blast",
                blast_file="in.blast", gene_lengths_file="lengths.txt")

    produced = (tmp_path / "exp_Zika_virusRNA.gplot").read_text()
    assert produced == (g / "blast2gplot.gplot.golden").read_text()
