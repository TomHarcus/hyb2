""" Port of bp2hyb.sh

Turns a list of scored base pairs (i<TAB>j<TAB>count, from ct2bps | histogram)
into a ranked .hyb of merged duplexes. 

trim_fraction (new, not in the legacy script; 0 = legacy behaviour): after merging,
base pairs at either end of a duplex are dropped while their count is below
trim_fraction x the duplex's best-supported pair. Merging joins any diagonally
adjacent pairs regardless of support, so one weakly supported pair could extend a
well-supported helix and, as a forced constraint, block a neighbouring helix (e.g.
U2 G1097-U1145, 1 chimera, extending the 1088-1154 helix of ~7,000).
Only ends are trimmed; weak pairs inside a helix are kept.

"""

import sys
from hyb2.folding.combine_hyb_merge_touching import combine_hyb_merge_touching

def bp2hyb(lines, trim_fraction=0.0):
    
    rows = [ln.rstrip("\n") for ln in lines if ln.strip() != ""]
    support = {}
    for ln in rows:
        i, j, count = ln.split("\t")[:3]
        support[(int(i), int(j))] = float(count)
    rows.sort(key=lambda ln: (int(ln.split("\t")[0]), ln))

    hyb = []
    for nr, ln in enumerate(rows, start=1):
        i, j, count = ln.split("\t")[:3]
        hyb.append("\t".join(
            [str(nr), ".", count, "RNA", ".", ".", i, i, ".",
             "RNA", ".", ".", j, j, "."]
        ))

    merged = combine_hyb_merge_touching(hyb)

    rescaled = []
    for ml in merged:
        f = ml.rstrip("\n").split("\t")
        if trim_fraction > 0:
            f = _trim_ends(f, support, trim_fraction)
        score = float(f[2]) * (1 + int(f[7]) - int(f[6]))
        f[2] = f"{score:.6g}"
        rescaled.append("\t".join(f))

    rescaled.sort(key=lambda ln: (-float(ln.split("\t")[2]), ln))

    for ln in rescaled:
        yield ln + "\n"


def _trim_ends(f, support, trim_fraction):
    """Drop weak pairs from both ends of one merged duplex (hyb fields f)."""
    s1, e1, e2 = int(f[6]), int(f[7]), int(f[13])
    counts = [support.get((s1 + k, e2 - k), 0.0) for k in range(e1 - s1 + 1)]
    cutoff = trim_fraction * max(counts)
    lo, hi = 0, len(counts) - 1
    while lo < hi and counts[lo] < cutoff:
        lo += 1
    while hi > lo and counts[hi] < cutoff:
        hi -= 1
    if (lo, hi) == (0, len(counts) - 1):
        return f
    kept = counts[lo:hi + 1]
    f = f[:]
    f[6], f[7] = str(s1 + lo), str(s1 + hi)
    f[12], f[13] = str(e2 - hi), str(e2 - lo)
    f[2] = f"{sum(kept) / len(kept):.6g}"   # mean count; x length below = total support
    return f


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    src = open(argv[0]) if argv else sys.stdin
    try:
        sys.stdout.writelines(bp2hyb(src))
    finally:
        if argv:
            src.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())