"""plot_all_cdm is a thin wrapper: glob *dg.hyb -> plot_hybrids_3 -> <stem><gene>.contact.txt
-> contact_density_map.R. The binning core (plot_hybrids_3) is golden-tested separately in
test_plot_hybrids_3.py, and the R plotting isn't golden-tested (same as the rest of the CDM
path). So this test pins the wrapper's own glue on the real test.ua.hyb:

  * the *dg.hyb glob picks the file up,
  * the ${f/hybrids_ua_dg.hyb/<gene>.contact.txt} substring rename is reproduced,
  * the written contact.txt equals plot_hybrids_3's (already-golden) single-gene output,
  * contact_density_map.R is invoked with (gene, limit) and no file arg.

The R call is stubbed so the test needs neither R nor the plotting scripts.
"""

import shutil

from hyb2.coverage import plot_all_cdm as mod
from hyb2.coverage.plot_hybrids_3 import plot_hybrids_3


def test_wrapper_glue(tier2_fixtures_dir, tmp_path, monkeypatch):
    # legacy plot_all_cdm operates on files matching *dg.hyb whose names embed the
    # literal "hybrids_ua_dg.hyb" substring the rename targets.
    hyb = tmp_path / "sample_hybrids_ua_dg.hyb"
    shutil.copy(tier2_fixtures_dir / "test.ua.hyb", hyb)

    calls = []
    monkeypatch.setattr(mod.subprocess, "run", lambda *a, **k: calls.append((a, k)))
    monkeypatch.chdir(tmp_path)

    mod.plot_all_cdm("Zika_virusRNA", 0.95)

    # ${f/hybrids_ua_dg.hyb/Zika_virusRNA.contact.txt}
    contact = tmp_path / "sample_Zika_virusRNA.contact.txt"
    assert contact.exists()

    with open(tier2_fixtures_dir / "test.ua.hyb") as fh:
        expected = plot_hybrids_3(fh, "Zika_virusRNA", "Zika_virusRNA", bin_size=10)
    assert contact.read_text() == expected

    # one R call: contact_density_map.R <gene> <limit>, no file argument
    assert len(calls) == 1
    argv = calls[0][0][0]
    assert argv[0] == "Rscript"
    assert argv[-2:] == ["Zika_virusRNA", "0.95"]


def test_cli_parses():
    args = mod.build_parser().parse_args(["-g", "Zika_virusRNA", "-q", "0.9"])
    assert args.gene == "Zika_virusRNA"
    assert args.limit == 0.9
