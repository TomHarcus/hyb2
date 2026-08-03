""" Port of bp_score.sh

For each folded structure listed in the coordinates manifest, extract its base
pairs from the .ct (offset-shifted local -> genome) into a .bps, then drive
make_varna_scores to produce that structure's per-position colour scores.
"""

import argparse
import sys

from hyb2.stages.make_VARNA_scores_2 import make_varna_scores, _output_name

def bp_score(coordinates, score, prefix):

    for line in coordinates:
        parts = line.split()

        if len(parts) < 4:
            continue

        nm, c1, c2, c3 = parts[0], parts[1], parts[2], parts[3]

        print(f"Name: {nm}")
        print(f"C1: {c1}")
        print(f"C2: {c2}")
        print(f"C3: {c3}")

        offset = int(c1)
        bps_path = nm.replace(".ct", ".bps")

        with open(nm) as ct, open(bps_path, "w") as out:
            for k, ctline in enumerate(ct):
                cols = ctline.split()
                if k == 0:
                    out.write(f"{cols[0]}\t{cols[4]}\n")
                elif int(cols[0]) < int(cols[4]):
                    out.write(f"{int(cols[0]) + offset}\t{int(cols[4]) + offset}\n")
        
        prefix_arg = prefix.replace(".hyb", "", 1) if prefix else None

        with open(score) as sc, open(bps_path) as bps:
            lines = make_varna_scores(sc, bps, length=int(c2), tail=int(c3), max_score=1000000)

        out_name = _output_name(score, bps_path, "", prefix_arg)
        with open(out_name, "w") as f:
            f.writelines(lines)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bp-score",
        description="Extract .bps from folded .ct(s) and drive make_varna_scores "
        "(port of bp_score.sh).",
    )
    p.add_argument("-c", dest="coord", required=True, metavar="coordinates.txt",
                   help="manifest: one row per structure -> 'ct_name offset len tail'")
    p.add_argument("-s", dest="score", required=True, metavar="basepair_scores.txt",
                   help="basepair_scores file passed through to make_varna_scores")
    p.add_argument("-p", dest="prefix", default=None, metavar="input.hyb",
                   help="optional output-name prefix (.hyb stripped)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    with open(args.coord) as c:
        bp_score(c, args.score, args.prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())
    