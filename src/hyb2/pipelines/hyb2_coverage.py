""" Port of hyb2_coverage


"""
import argparse
import subprocess
import sys
from pathlib import Path
from hyb2.stages.plot_hybrids_3 import plot_hybrids_3, swap_gene1_to_arm1
from hyb2 import config

def hyb2_coverage(in_hyb, gene_1, gene_2, limit, x1, x2, y1, y2):

    if not gene_2 and not x1 and not x2 and not y1 and not y2:

        print(f"Plotting contact density map of {gene_1}...")

        contact = in_hyb.replace(".hyb", f".{gene_1}.contact.txt")

        with open(in_hyb) as fin:
            awk_out = plot_hybrids_3(fin, gene_1, gene_1, bin_size=10)

        rows = awk_out.splitlines()

        Path(contact).write_text(rows[0] + "\nx\ty\tcount\n" + "\n".join(rows[1:]) + "\n")

        subprocess.run(["Rscript", config.rscript("contact_density_map_indiv.R"),
                        contact,
                        str(limit)],
                        check=True)
        

        print(f"Contact density map of {gene_1} saved")

    elif gene_2 and not x1 and not x2 and not y1 and not y2:

        print(f"Plotting contact density map of {gene_1} and {gene_2}...")

        contact = in_hyb.replace(".hyb", f".{gene_1}-{gene_2}.contact.txt")

        with open(in_hyb) as fin:
            swapped = swap_gene1_to_arm1(fin, gene_1)
            
        awk_out = plot_hybrids_3(swapped, gene_1, gene_2, bin_size=10)

        rows = awk_out.splitlines()
        
        Path(contact).write_text(rows[0] + "\nx\ty\tcount\n" + "\n".join(rows[1:]) + "\n")

        subprocess.run(["Rscript", config.rscript("cdm_2genes.R"),
                        contact,
                        str(limit)],
                        check=True)

        print(f"Contact density map of {gene_1} and {gene_2} saved")

    if not gene_2 and x1 and x2 and y1 and y2:

        print(f"Plotting zoomed in contact density map of {gene_1}...")

        contact = in_hyb.replace(".hyb", f".{gene_1}.contact.txt")
    
        subprocess.run(["Rscript", config.rscript("cdm_indiv_zoom.R"),
                        contact, str(x1), str(x2), str(y1), str(y2),
                        str(limit)],
                        check=True)

        print(f"Zoomed in contact density map of {gene_1} saved")

    elif gene_2 and x1 and x2 and y1 and y2:

        print(f"Plotting zoomed in contact density map of {gene_1} and {gene_2}...")

        contact = in_hyb.replace(".hyb", f".{gene_1}-{gene_2}.contact.txt")
            
        subprocess.run(["Rscript", config.rscript("cdm_indiv_zoom.R"),
                        contact,str(x1), str(x2), str(y1), str(y2),
                        str(limit), gene_1, gene_2],
                        check=True)

        print(f"Zoomed in contact density map of {gene_1} and {gene_2} saved")


def build_parser() -> argparse.ArgumentParser:
    
    p = argparse.ArgumentParser(
        prog="hyb2-coverage",
        description="plot contact-density maps (wraps the R plotting scripts)",
        add_help=False,
    )
    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-i", dest="in_hyb", required=True, metavar="INPUT.HYB", help="Input HYB (required)")
    p.add_argument("-a", dest="gene_1", required=True, metavar="GENE_1", help="gene to plot (required)")
    p.add_argument("-b", dest="gene_2", default=None, metavar="GENE_2", help="second gene (two-gene map)")
    p.add_argument("-q", dest="limit", type=float, default=0.95, help="upper quantile for heatmap contrast (default=0.95)")
    p.add_argument("-w", dest="x1", type=int, default=None, help="zoom window X start")
    p.add_argument("-x", dest="x2", type=int, default=None, help="zoom window X end")
    p.add_argument("-y", dest="y1", type=int, default=None, help="zoom window Y start")
    p.add_argument("-z", dest="y2", type=int, default=None, help="zoom window Y end")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    hyb2_coverage(
        args.in_hyb, args.gene_1, args.gene_2, args.limit,
        args.x1, args.x2, args.y1, args.y2,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())









