""" Port of make_VARNA_scores_2.sh

Maps a folded structure's base pairs to per-nucleotide colour scores (the
experimental support of the pair each nucleotide is in), for VARNA shading.
Core transform is pure; the -s/-p suffix/prefix only shape the output filename.
"""

import argparse
import sys


def make_varna_scores(basepair_scores, bps_lines, *, length, tail, max_score):
    fnd = {}

    for line in basepair_scores:
        i, j, c = line.split()[:3]
        c = int(c)
        fnd[(i, j)] = fnd.get((i, j), 0) + c
        fnd[(j, i)] = fnd.get((j, i), 0) + c

    score = [0] * (length + 1)

    for line in bps_lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        i, j = parts[0], parts[1]

        s = fnd.get((i, j))
        if not s:
            continue
        s = min(s, max_score)

        for pos in (int(i), int(j)):
            if 1 <= pos <= length:
                score[pos] = s

    out = [f"{score[p]}\n" for p in range(1, length + 1)]
    return out[-tail:] if tail > 0 else out


def _output_name(input_path: str, bps_path: str, suffix: str, prefix: str | None) -> str:
    base = prefix if prefix else input_path.replace(".basepair_scores.txt", "", 1)
    bps_base = bps_path.replace(".bps", "", 1)
    return f"{base}__{bps_base}{suffix}.VARNA_scores.txt"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="make-varna-scores",
        description="Per-nucleotide VARNA colour scores from folded pairs + "
        "experimental support (port of make_VARNA_scores_2.sh).",
    )
    p.add_argument("-i", dest="input", required=True, metavar="basepair_scores.txt",
                   help="basepair_scores file (i j count)")
    p.add_argument("-b", dest="bps", required=True, metavar="in.bps",
                   help="folded structure's base pairs (i j)")
    p.add_argument("-l", dest="length", type=int, required=True, help="sequence length")
    p.add_argument("-m", dest="max_score", type=int, required=True, help="cap scores at this value")
    p.add_argument("-t", dest="tail", type=int, default=0, help="keep only the last TAIL positions (0=all)")
    p.add_argument("-s", dest="suffix", default="", help="suffix appended to the output filename")
    p.add_argument("-p", dest="prefix", default=None, help="prefix for the output filename (else derived from -i)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    with open(args.input) as fin, open(args.bps) as bin_:
        lines = make_varna_scores(fin, bin_, length=args.length, tail=args.tail,
                                  max_score=args.max_score)
    out_name = _output_name(args.input, args.bps, args.suffix, args.prefix)
    with open(out_name, "w") as fout:
        fout.writelines(lines)
    return 0


if __name__ == "__main__":
    sys.exit(main())


