"""Port of bin/histogram.pl 

Tallies counts of identical input lines/keys and emits a frequency table. Used
twice by the sam_composition pipeline (gene-pair counts, single-read top-hit
counts) and again in Tier 2. The Perl uses only core/CPAN modules
(List::Util, etc.) -> Python builtins; no project-file dependencies.

"""

import sys
from typing import Iterable, Iterator


def histogram(lines: Iterable[str]) -> Iterator[str]:
    """Count identical input lines and yield '<key>\\t<count>' rows, ordered by
    count descending.

    Only histogram.pl's default (no-flag) mode is ported - every HYB2 call site
    invokes it bare. The --numeric/--step/--noise/--tail/--percentual/
    --include-zero options are intentionally omitted (unused dead code).

    Tie note: the Perl orders equal counts by hash order (non-deterministic), so
    the legacy output isn't stable on ties. This keeps first-appearance order for
    ties (deterministic); the golden fixtures have no ties, so parity holds.
    """
    hits: dict[str, int] = {}
    for line in lines:
        key = line.rstrip("\n")  # Perl chomp
        hits[key] = hits.get(key, 0) + 1

    for key, count in sorted(hits.items(), key=lambda kv: kv[1], reverse=True):
        yield f"{key}\t{count}\n"


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
