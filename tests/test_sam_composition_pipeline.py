"""End-to-end golden-diff parity for the whole spine.

Runs the Python pipeline on the real test.sam and diffs every intermediate
against fixtures/sam_composition_run/. Catches wiring bugs the per-stage tests
can't. Remove xfail once pipelines.sam_composition.run is implemented.

Note: step 9 (composition pies) is not asserted -- hyb2_composition_pies.py
isn't available locally (get from Grzegorz).
"""

import shutil

import pytest

from hyb2.pipelines.sam_composition import run

# .hyb carries a comment header (compare data rows only); everything else is byte-exact.
GOLDEN_OUTPUTS = [
    "test.blast",
    "test.collapse.blast",
    "test_mtophits.blast",
    "test_mtophits.ref",
    "test.hyb",            # header-stripped
    "test.ua.hyb",
    "test.ua.hyb_stats_by_gene.txt",
    "test_tophit_by_gene.txt",
]


def test_sam_composition_pipeline_matches_golden(fixtures_dir, tmp_path, strip_hyb_header):
    sam = tmp_path / "test.sam"
    shutil.copy(fixtures_dir / "test.sam", sam)

    run(str(sam))

    for name in GOLDEN_OUTPUTS:
        produced = (tmp_path / name).read_text()
        golden = (fixtures_dir / name).read_text()
        if name == "test.hyb":
            produced, golden = strip_hyb_header(produced), strip_hyb_header(golden)
        assert produced == golden, f"mismatch in {name}"
