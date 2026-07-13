"""Port of collapse_blast_2.sh (Tier 1 #3).

Collapses reads by (gene name, mapped sequence), tallying duplicate counts into
the read-id field, keeping a stable first occurrence per (gene, sequence).

Grzegorz's newer version (in hyb2_scripts/) is a memory-safe, disk-backed
rewrite with byte-identical output to the original; the parity oracle was built
with it. Decide during porting: reimplement in Python, or keep shelling out to
that script. On the real 17 MB blast the dedup is the memory-sensitive step, so
match its streaming behaviour rather than loading everything into a dict.

Parity: fixtures/sam_composition_run/test_mtophits.blast (after mtophits)
"""

from __future__ import annotations

import sys
from typing import Iterable, Iterator


def collapse_blast(lines: Iterable[str]) -> Iterator[str]:
    raise NotImplementedError


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: collapse_blast <in.blast>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        for out in collapse_blast(f):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
