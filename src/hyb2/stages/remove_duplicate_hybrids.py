"""Port of bin/remove_duplicate_hybrids_hOH5_2.pl (Tier 1 #5).

Takes the ranking reference (from create_reference_file) plus a .hyb, keeps
exactly one best hybrid per read ID. Ranking criteria (in order): sum of the
two bits' e-values; then miRNA-mRNA hybrids (when PREFER_MIM); then rank of the
higher- then lower-ranked bit in the reference; else keep the first seen.

Uses the shared Hybrid data class (models.Hybrid), the one cross-file
dependency in the Tier 1 set.

Parity: fixtures/sam_composition_run/test.ua.hyb
"""

from __future__ import annotations

import sys
from typing import Iterable, Iterator

from hyb2.models import Hybrid  # noqa: F401  (used once implemented)


def remove_duplicate_hybrids(
    ref_lines: Iterable[str],
    hyb_lines: Iterable[str],
    *,
    prefer_mim: bool = True,
) -> Iterator[str]:
    """Yield the single best .hyb line per read ID."""
    raise NotImplementedError


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    # TODO(tier1): parse leading PREFER_MIM=.. switch, then <ref> <hyb> paths.
    raise NotImplementedError


if __name__ == "__main__":
    sys.exit(main())
