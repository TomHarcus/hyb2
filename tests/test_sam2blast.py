"""Golden-diff parity: sam2blast vs fixtures/sam_composition_run/test.blast.

Remove the xfail marker once bin/sam2blast_3 is adopted into the stage.
"""

import pytest

from hyb2.stages.sam2blast import sam2blast


def test_sam2blast_matches_golden(fixtures_dir):
    with (fixtures_dir / "test.sam").open() as f:
        produced = "".join(sam2blast(f))
    assert produced == (fixtures_dir / "test.blast").read_text()
