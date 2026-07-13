"""Golden-diff parity: get_mtop_hybrids vs fixtures/.../test.hyb.

Compared on data rows only -- the .hyb comment header embeds the invocation
path/args, which differ between the Perl and the port (see conftest.data_lines).
Params match what hyb2_sam_composition.sh passes. Remove xfail once implemented.
"""

import pytest

from hyb2.stages.get_mtop_hybrids import get_mtop_hybrids


@pytest.mark.xfail(raises=NotImplementedError, strict=False, reason="Tier 1: stage not yet ported")
def test_get_mtop_hybrids_matches_golden(fixtures_dir, strip_hyb_header):
    with (fixtures_dir / "test.blast").open() as f:
        produced = "".join(
            get_mtop_hybrids(f, blast_threshold=0.1, mode=2, max_overlap=4, max_hits=10)
        )
    expected = (fixtures_dir / "test.hyb").read_text()
    assert strip_hyb_header(produced) == strip_hyb_header(expected)
