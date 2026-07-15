""" Port of hyb2constraints.pl

Takes hyb file as input, predicte RNA stem coordinates
Output file in form "F 1 100 5", which can be used as a list of constraints for hybrid-ss-min

***UNTESTED RIGHT NOW***

"""

import sys
from typing import Iterable, Iterator


def hyb2constraints(lines: Iterable[str]) -> Iterator[str]:
    for line in lines:
        line = line.rstrip("\n")
        # skip comment/blank/malformed lines
        columns = line.split("\t")
        if not line or len(columns) < 14:
            continue
        arm1_start, arm1_end = int(columns[6]), int(columns[7])
        arm2_start, arm2_end = int(columns[12]), int(columns[13])

        if overlap(arm1_start, arm1_end, arm2_start, arm2_end) > 0:
            continue
        
        coords = sorted([arm1_start, arm1_end, arm2_start, arm2_end])

        yield f"F\t{coords[0]}\t{coords[3]}\t{coords[1] - coords[0] + 1}\n"


def overlap(a, b, c, d):
    return 1 + min(b, d) - max(a, c)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: hyb2constraints <in.hyb>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        for out in hyb2constraints(f):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())