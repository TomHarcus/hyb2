"""Smoke test for the CPLfold backend wiring (comrades_fold -r cplfold).

Not a golden diff: cplfold output depends on the recompiled HotKnots binaries
(architecture-specific), so exact bytes aren't reproducible across machines. This
just drives the full parameter chain end-to-end with NON-DEFAULT bonus + search
flags and asserts it produces a valid, non-empty structure -- catching wiring
regressions like the max_phase2 / FOLD-vs-fold ones (params accepted by run() but
dropped before the fold, or a valid-looking flag that crashes).

Inputs: the mini fragment + its basepair scores from
scripts/generate_folding_baseline.sh (window 10000-10300, len 301). Skips unless
ViennaRNA + CPLfold + a runnable HotKnots binary are all present.
"""

from hyb2.folding.comrades_fold import run


def test_cplfold_backend_runs_with_nondefault_params(
    folding_fixtures_dir, vienna_bin, cplfold_env, tmp_path, monkeypatch
):
    g = folding_fixtures_dir
    (tmp_path / "frag.fasta").write_text((g / "comrades_fold.frag.fasta").read_text())
    (tmp_path / "scores.txt").write_text((g / "make_varna.basepair_scores.txt").read_text())

    monkeypatch.chdir(tmp_path)
    # non-default alpha/normalize + every search-breadth knob -> exercises the
    # whole run() -> _fold_cplfold -> two_phase_pseudoknot_fold param chain
    vienna, ct = run(
        "unused", "frag.fasta", fold="cplfold",
        basepair_scores="scores.txt", begin=10000, end=10300,
        vienna_bin=vienna_bin,
        alpha=0.3, beta=0.1, normalize="raw", beam_size=80,
        energy_delta=8.0, max_phase1=12, max_phase2=6, energy_model="DP09",
    )

    vlines = (tmp_path / vienna).read_text().splitlines()
    assert len(vlines) >= 3                       # header, sequence, structure(+energy)
    structure = vlines[2].split()[0]
    assert structure and set(structure) <= set("().[]{}<>")   # valid (pseudoknot-aware) dot-bracket
    # energy guard turns a HotKnots failure into a raise, so a returned .vienna
    # must carry a real numeric energy, never "None"
    assert "None" not in vlines[2]

    # .ct is non-trivial (header line + one row per base)
    assert len((tmp_path / ct).read_text().splitlines()) > len(structure)
