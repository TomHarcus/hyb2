"""Golden-diff parity: remove_duplicate_hybrids vs fixtures/.../test.ua.hyb.

Takes the ranking ref + the .hyb; the dedup output has no comment header, so a
plain byte-diff is the contract. Remove xfail once implemented (also implement
models.Hybrid, which this stage depends on).
"""

import pytest

from hyb2.stages.remove_duplicate_hybrids import remove_duplicate_hybrids


@pytest.mark.xfail(raises=NotImplementedError, strict=False, reason="Tier 1: stage not yet ported")
def test_remove_duplicate_hybrids_matches_golden(fixtures_dir):
    with (fixtures_dir / "test_mtophits.ref").open() as ref, (fixtures_dir / "test.hyb").open() as hyb:
        produced = "".join(remove_duplicate_hybrids(ref, hyb, prefer_mim=True))
    assert produced == (fixtures_dir / "test.ua.hyb").read_text()
