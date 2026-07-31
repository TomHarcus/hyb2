""" Port of plot_differential_map

Finds significant differences between two conditions 
"""

import argparse, sys
from pathlib import Path
from hyb2.pipelines.hyb2_compare import _differential_map, _num
from hyb2.stages.make_hybrid_annotation_table import make_hybrid_annotation_table


def plot_differential_map(a1, a2, a3, a4, b1, b2, b3, b4, minimum=2, range=100, limit=0.95):

    files = [f for f in (a1, a2, a3, a4, b1, b2, b3, b4) if f]

    out = f"{a1.replace('.contact.txt', '')}-{b1.replace('.contact.txt', '')}"

    cond_one = [f for f in (a1, a2, a3, a4) if f]
    cond_two = [f for f in (b1, b2, b3, b4) if f]
    conditions = ["condition_one"] * len(cond_one) + ["condition_two"] * len(cond_two)

    
    for f in files:
        name = f.replace(".contact.txt", "")
        rows = [l.rstrip("\n").split("\t") for l in open(f) if not l.startswith("#")]
        rows.sort(key=lambda c: _num(c[2]), reverse=True)
        Path(f"{name}.forTable.txt").write_text(
            "\n".join(f"{c[0]}_{c[1]}\t{name}={c[2]}" for c in rows) + "\n"
        )

    fortable_files = [f.replace(".contact.txt", ".forTable.txt") for f in files]
    names = [f.replace(".contact.txt", "") for f in files]

   
    header = ("#seq_ID\t" + "\t".join(names)).replace("-", "_")
    Path(f"{out}.table.txt").write_text(header + "\n")
    body = [r.replace("NA", "0")
            for r in make_hybrid_annotation_table(fortable_files).splitlines()[1:]]
    with open(f"{out}.table.txt", "a") as fh:
        fh.write("\n".join(body) + "\n")

  
    names_table = "\n".join(
        f"{name}\t{cond}".replace("-", "_") for name, cond in zip(names, conditions)
    ) + "\n"
    Path(f"{out}_names.table").write_text(names_table)

    # NOTE (deviation): _differential_map calls differential_coverage_map.R with
    # (heatmap, out); the legacy plot_differential_map passed (heatmap, A1, B1). Same
    # deviation carried by hyb2_compare -- parameterize the helper's trailing args if
    # byte-parity with the legacy R call matters.
    _differential_map(out, minimum, range, limit, names)


def build_parser() -> argparse.ArgumentParser:

    p = argparse.ArgumentParser(
        prog="plot-differential-map",
        description="finds significant differences between two sets of data",
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
    p.add_argument("-m", dest="minimum", type=int, default=2, help="minimum reads for DESeq2 (default=2)")
    p.add_argument("-r", dest="range", type=int, default=100, help="range to define an interaction (default=100)")
    p.add_argument("-p", dest="limit", type=float, default=0.95, help="upper quantile for heatmap contrast (default=0.95)")
    return p


def main(argv: list[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    plot_differential_map(a.a1, a.a2, a.a3, a.a4, a.b1, a.b2, a.b3, a.b4,
                          a.minimum, a.range, a.limit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
