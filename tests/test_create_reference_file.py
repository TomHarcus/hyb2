"""Golden-diff parity: create_reference_file vs fixtures/.../test_mtophits.ref.

"""

import pytest

from hyb2.chimera.create_reference_file import create_reference


def test_create_reference_matches_golden(fixtures_dir):
    with (fixtures_dir / "test_mtophits.blast").open() as f:
        produced = "".join(create_reference(f))
    assert produced == (fixtures_dir / "test_mtophits.ref").read_text()
