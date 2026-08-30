"""Golden-diff parity: sam2blast vs fixtures/sam_composition_run/test.blast.

"""

import pytest

from hyb2.chimera.sam2blast import sam2blast


def test_sam2blast_matches_golden(fixtures_dir):
    with (fixtures_dir / "test.sam").open() as f:
        produced = "".join(sam2blast(f))
    assert produced == (fixtures_dir / "test.blast").read_text()
