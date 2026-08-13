""" Port of bp2hyb.sh

Turns a list of scored base pairs (i<TAB>j<TAB>count, from ct2bps | histogram)
into a ranked .hyb of merged duplexes. 

"""

import sys
from hyb2.folding.combine_hyb_merge_touching import combine_hyb_merge_touching

def bp2hyb(lines):
    
    rows = [ln.rstrip("\n") for ln in lines if ln.strip() != ""]
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
        score = float(f[2]) * (1 + int(f[7]) - int(f[6]))
        f[2] = f"{score:.6g}"
        rescaled.append("\t".join(f))

    rescaled.sort(key=lambda ln: (-float(ln.split("\t")[2]), ln))

    for ln in rescaled:
        yield ln + "\n"


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