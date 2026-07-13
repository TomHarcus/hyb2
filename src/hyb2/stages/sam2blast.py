"""Adopt bin/sam2blast_3 (Tier 1, already Python).

bin/sam2blast_3 is already Python 3 -- this is an ADOPT, not a translation:
move its logic here (SAM -> blast, with the mapped-sequence column and the
bit-score / e-value computation), add the iterable/main shape used by the other
stages, and keep byte-parity with the current script.

Parity: fixtures/sam_composition_run/test.blast
"""

from __future__ import annotations

import sys
from typing import Iterable, Iterator


def sam2blast(lines: Iterable[str]) -> Iterator[str]:
    # TODO(tier1): port bin/sam2blast_3.print_line / main here.
    raise NotImplementedError


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: sam2blast <in.sam>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        for out in sam2blast(f):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
