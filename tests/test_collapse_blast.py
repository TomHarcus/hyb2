"""Golden-diff parity: collapse_blast vs fixtures/.../test.collapse.blast.

Remove the xfail marker once the stage is implemented.
"""

import pytest

from hyb2.chimera.collapse_blast import collapse_blast


def test_collapse_blast_matches_golden(fixtures_dir, tmp_path):
    out = tmp_path / "collapse.blast"
    collapse_blast(str(fixtures_dir / "test.blast"), str(out))
    assert out.read_text() == (fixtures_dir / "test.collapse.blast").read_text()
