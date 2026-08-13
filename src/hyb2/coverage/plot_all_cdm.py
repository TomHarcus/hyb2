""" Port of plot_all_cdm

"""

import glob, subprocess, argparse, sys
from pathlib import Path
from hyb2.coverage.plot_hybrids_3 import plot_hybrids_3
from hyb2.tools import config

def plot_all_cdm(gene, limit=0.95):

    for file in glob.glob("*dg.hyb"):
        with open(file) as f:
            output = plot_hybrids_3(f, gene_1=gene, gene_2=gene, bin_size=10, threshold_num=0)

        file_out = file.replace("hybrids_ua_dg.hyb", f"{gene}.contact.txt")

        Path(file_out).write_text(output)

    subprocess.run(["Rscript", config.rscript("contact_density_map.R"),
                            gene,
                            str(limit)],
                            check=True)


def build_parser() -> argparse.ArgumentParser:

    p = argparse.ArgumentParser(
        prog="plot-all-cdm",
        description="plot contact-density maps for every *dg.hyb in the current directory",
        add_help=False,
    )
    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-g", dest="gene", required=True, metavar="GENE", help="gene to plot (required)")
    p.add_argument("-q", dest="limit", type=float, default=0.95, help="upper quantile for heatmap contrast (default=0.95)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    plot_all_cdm(args.gene, args.limit)
    return 0


if __name__ == "__main__":
    sys.exit(main())