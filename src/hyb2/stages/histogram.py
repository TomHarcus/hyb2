"""Port of bin/histogram.pl (Tier 1 #2).

Tallies counts of identical input lines/keys and emits a frequency table. Used
twice by the sam_composition pipeline (gene-pair counts, single-read top-hit
counts) and again in Tier 2. The Perl uses only core/CPAN modules
(List::Util, etc.) -> Python builtins; no project-file dependencies.

Parity: fixtures/sam_composition_run/test.ua.hyb_stats_by_gene.txt
        fixtures/sam_composition_run/test_tophit_by_gene.txt
"""

from __future__ import annotations

import sys
from typing import Iterable, Iterator


def histogram(lines: Iterable[str]) -> Iterator[str]:
    """Yield '<key>\\t<count>' rows. Confirm ordering/tie behaviour against the
    Perl before finalising (histogram.pl has several output options)."""
    raise NotImplementedError


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    src = open(argv[0]) if argv else sys.stdin
    try:
        for out in histogram(src):
            sys.stdout.write(out)
    finally:
        if src is not sys.stdin:
            src.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
