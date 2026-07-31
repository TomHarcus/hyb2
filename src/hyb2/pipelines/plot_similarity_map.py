""" Port of plot_similarity_map

"""

import subprocess, argparse, sys
from pathlib import Path
from hyb2.stages.similarity import similarity_contact
from hyb2 import config

def plot_similarity_map(a1, a2, a3, a4, b1, b2, b3, b4, limit):

    files = [f for f in (a1, a2, a3, a4, b1, b2, b3, b4) if f]

    out = f"{a1.replace('.contact.txt', '')}-{b1.replace('.contact.txt', '')}"

    merged = "".join(open(f).read() for f in files)
    Path(f"{out}.merge.txt").write_text(merged)

    out_rows = similarity_contact(open(f"{out}.merge.txt").read().splitlines())
    Path(f"{out}.contact.txt").write_text("\n".join(out_rows) + "\n")

    subprocess.run(["Rscript", config.rscript("similarity_heatmap.R"),
                    f"{out}.contact.txt", str(limit)], check=True)

def build_parser() -> argparse.ArgumentParser:

    p = argparse.ArgumentParser(
        prog="plot-similarity-map",
        description="locate overlapping coordinates and print the minimum value",
        add_help=False,
    )
    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-a", dest="a1", required=True, metavar="A1.contact.txt", help="condition one, replicate 1 (required)")
    p.add_argument("-b", dest="a2", required=True, metavar="A2.contact.txt", help="condition one, replicate 2 (required)")
    p.add_argument("-c", dest="a3", default=None, metavar="A3.contact.txt", help="condition one, replicate 3")
    p.add_argument("-d", dest="a4", default=None, metavar="A4.contact.txt", help="condition one, replicate 4")
    p.add_argument("-i", dest="b1", required=True, metavar="B1.contact.txt", help="condition two, replicate 1 (required)")
    p.add_argument("-j", dest="b2", required=True, metavar="B2.contact.txt", help="condition two, replicate 2 (required)")
    p.add_argument("-k", dest="b3", default=None, metavar="B3.contact.txt", help="condition two, replicate 3")
    p.add_argument("-l", dest="b4", default=None, metavar="B4.contact.txt", help="condition two, replicate 4")
    p.add_argument("-p", dest="limit", type=float, default=0.95, help="upper quantile for heatmap contrast (default=0.95)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    plot_similarity_map(args.a1, args.a2, args.a3, args.a4,
                        args.b1, args.b2, args.b3, args.b4, args.limit)
    return 0


if __name__ == "__main__":
    sys.exit(main())