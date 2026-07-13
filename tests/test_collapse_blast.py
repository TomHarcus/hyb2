"""Golden-diff parity: collapse_blast vs fixtures/.../test.collapse.blast.

Remove the xfail marker once the stage is implemented.
"""

import pytest

from hyb2.stages.collapse_blast import collapse_blast


def test_collapse_blast_matches_golden(fixtures_dir):
    with (fixtures_dir / "test.blast").open() as f:
        produced = "".join(collapse_blast(f))
    assert produced == (fixtures_dir / "test.collapse.blast").read_text()
