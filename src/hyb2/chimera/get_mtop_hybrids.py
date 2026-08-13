"""Port of bin/get_mtop_hybrids.pl - the core chimera caller.

Reads a .blast (raw, NOT collapsed), groups by read ID, applies MODE / overlap
/ e-value cutoffs, and emits .hyb records for reads with >1 valid bit. Self-
contained in the original (defines its own min/max/overlap subs).

Legacy FOO=bar switches become keyword args with the same defaults. Note the
sam_composition pipeline overrides these to blast_threshold=0.1, mode=2,
max_overlap=4, max_hits=10.

"""

import re
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

    name = None
    arr = []
    nbits = 0
    min_start = max_end = 0
    curr_overlap = 0
    curr_blast_threshold = 0

    for line in lines:
        fld = line.split("\t")

        if fld[0] != name:
            yield from _flush(arr, nbits, max_hits, max_overlap)

            arr = []
            name = fld[0]
            min_start = int(fld[6])
            max_end = int(fld[7])
            curr_overlap = 0
            curr_blast_threshold = 0
            nbits = 0
        
        else:
            curr_overlap = overlap(min_start, max_end, int(fld[6]), int(fld[7]))
            if float(fld[10]) <= blast_threshold:
                min_start = min(min_start, int(fld[6]))
                max_end = max(max_end, int(fld[7]))
        
        if int(fld[8]) > int(fld[9]) and mode in (2, 3):
            continue
        # Faithful quirk: this push happens before the overlap
        # e-value filters below, so a line can land in arr even
        # though a later filter continues it. Combined with the second push at the 
        # != curr_blast_threshold branch below, a strictly improving e-value gets
        # appended twice. Both reproduced deliberately - they change which chimera
        # pairs get emitted
        if float(fld[10]) <= curr_blast_threshold:
            arr.append(line)

        if float(fld[10]) > blast_threshold:
            continue

        if mode == 0 and curr_overlap > max_overlap:
            continue

        if mode == 1 and (curr_overlap > max_overlap or curr_overlap < -max_overlap):
            continue

        if mode == 2 and (curr_overlap > max_overlap or curr_overlap < -max_overlap):
            continue

        if mode == 3 and curr_overlap > max_overlap:
            continue

        nbits += 1

        if float(fld[10]) != curr_blast_threshold:
            curr_blast_threshold = float(fld[10])
            arr.append(line)

    yield from _flush(arr, nbits, max_hits, max_overlap)


def print_hyb_format(ordered, max_overlap):
    a = [L.rstrip("\n").split("\t") for L in ordered]

    n = len(a)
    for i in range(n):
        for j in range(i+1, n):
            abs_overlap = abs(overlap(int(a[i][6]), int(a[i][7]), int(a[j][6]), int(a[j][7])))

            if abs_overlap <= max_overlap:
                yield "\t".join([
                    a[i][0], ".", ".",
                    a[i][1], a[i][6], a[i][7], a[i][8], a[i][9], a[i][10],
                    a[j][1], a[j][6], a[j][7], a[j][8], a[j][9], a[j][10],
                ]) + "\n"

def _flush(arr, nbits, max_hits, max_overlap):
    if nbits > 1 and len(arr) <= max_hits:
        ordered = sorted(arr, key=lambda L: int(L.split("\t")[6]))
        yield from print_hyb_format(ordered, max_overlap)

def overlap(a, b, c, d):
    return 1 + min(b, d) - max(a, c)


_SWITCHES = {
    "BLAST_THRESHOLD": ("blast_threshold", float),
    "MODE": ("mode", int),
    "MAX_OVERLAP": ("max_overlap", int),
    "MAX_HITS_PER_SEQUENCE": ("max_hits", int),
    "OUTPUT_FORMAT": ("output_format", str),
}


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    kwargs = {}
    while argv and re.match(r"^[A-Za-z_0-9]+=", argv[0]):   # leading FOO=bar switches
        key, _, val = argv[0].partition("=")
        if key in _SWITCHES:
            name, conv = _SWITCHES[key]
            kwargs[name] = conv(val)
        argv = argv[1:]
    if len(argv) != 1:
        print("usage: get_mtop_hybrids [BLAST_THRESHOLD=..] [MODE=..] [MAX_OVERLAP=..] "
              "[MAX_HITS_PER_SEQUENCE=..] [OUTPUT_FORMAT=..] <in.blast>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        for out in get_mtop_hybrids(f, **kwargs):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
