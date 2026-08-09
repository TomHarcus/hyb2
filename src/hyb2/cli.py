"""CLI skeleton reserved for the top-level `hyb2` mega-orchestrator port.

Flags here mirror bin/hyb2's getopts string ("i:d:o:v:m:h:a:b:q:x:y:l:j:e:r:")

"""


import argparse
import sys


def build_parser() -> argparse.ArgumentParser:
    # add_help=False because the legacy bin/hyb2 getopts string uses -h for
    # max-hits-per-sequence, not help; --help still works via the explicit
    # long-only alias added below.
    parser = argparse.ArgumentParser(
        prog="hyb2-py",
        description="Python port of the HYB2 RNA proximity-ligation pipeline (in progress).",
        add_help=False,
    )
    parser.add_argument("--help", action="help", help="Show this help message and exit")
    parser.add_argument("-i", dest="in_file", metavar="INPUT", help="Input fastq/SAM/hyb file")
    parser.add_argument("-d", dest="db_1", metavar="FASTA", help="Reference fasta used for mapping")
    parser.add_argument("-o", dest="out", metavar="OUTPUT_ID", help="Output ID")
    parser.add_argument("-v", dest="hval", type=float, default=0.1, help="BLAST threshold (default=0.1)")
    parser.add_argument("-m", dest="gmax", type=int, default=4, help="Maximum overlap (default=4)")
    parser.add_argument("-h", dest="hmax", type=int, default=10, help="Maximum hits per sequence (default=10)")
    parser.add_argument("-a", dest="gene_1", metavar="GENE_ID", help="Gene ID of interest")
    parser.add_argument("-b", dest="gene_2", metavar="GENE_ID", help="Second gene ID of interest")
    parser.add_argument("-q", dest="limit", type=float, default=0.95, help="Heatmap upper quantile (default=0.95)")
    parser.add_argument("-x", dest="x_coord", type=int, help="Start coordinate of 1st strand/gene")
    parser.add_argument("-y", dest="y_coord", type=int, help="Start coordinate of 2nd strand/gene")
    parser.add_argument("-l", dest="length", type=int, help="Length of fragments")
    parser.add_argument("-j", dest="varna", metavar="VARNAcmd.jar", help="Path to VARNAcmd.jar")
    parser.add_argument(
        "-e",
        dest="energy",
        type=int,
        default=0,
        help="Calculate folding energy: 1 on, 0 off (default=0)",
    )
    parser.add_argument(
        "-r",
        dest="fold",
        default="cplfold",
        choices=["cplfold", "vienna", "unafold", "0", "1"]
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    from hyb2.pipelines import hyb2

    if args.in_file is None:
        return hyb2.print_help()

    hyb2.run(args.in_file, args.db_1, args.out,
             hmax=args.hmax, blast_threshold=args.hval, max_overlap=args.gmax,
             gene_1=args.gene_1, gene_2=args.gene_2, limit=args.limit,
             x_coord=args.x_coord, y_coord=args.y_coord, length=args.length,
             varna=args.varna, fold=args.fold)

    return 0


if __name__ == "__main__":
    sys.exit(main())
