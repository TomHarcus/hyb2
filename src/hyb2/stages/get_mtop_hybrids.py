"""Port of bin/get_mtop_hybrids.pl (Tier 1 #6) -- the core chimera caller.

Reads a .blast (raw, NOT collapsed), groups by read ID, applies MODE / overlap
/ e-value cutoffs, and emits .hyb records for reads with >1 valid bit. Self-
contained in the original (defines its own min/max/overlap subs).

Legacy FOO=bar switches become keyword args with the same defaults. Note the
sam_composition pipeline overrides these to blast_threshold=0.1, mode=2,
max_overlap=4, max_hits=10.

Parity: fixtures/sam_composition_run/test.hyb
"""

from __future__ import annotations

import sys
from typing import Iterable, Iterator


def get_mtop_hybrids(
    lines: Iterable[str],
    *,
    blast_threshold: float = 0.001,
    mode: int = 2,
    max_overlap: int = 4,
    max_hits: int = 10,
    output_format: str = "HYB",
) -> Iterator[str]:
    """Call chimeras from blast lines. Streams input; yields output lines
    (including the leading '#'-comment header the Perl prints)."""
    raise NotImplementedError


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    # TODO(tier1): accept BLAST_THRESHOLD=.. MODE=.. etc. as leading FOO=bar
    # args for drop-in compatibility, then the trailing .blast path.
    raise NotImplementedError


if __name__ == "__main__":
    sys.exit(main())
