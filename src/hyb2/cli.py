"""CLI skeleton reserved for the top-level `hyb2` mega-orchestrator port.

This is the placeholder for porting `bin/hyb2` specifically -- the big
orchestrator that fans out to mapping, sam2hyb, folding, coverage, compare,
etc. (Tier 3 / dead-last in MIGRATION_ORDER.md). It is NOT the chimera-calling
spine: that is `hyb2_sam_composition.sh`, already ported and runnable via
`hyb2.pipelines.sam_composition` (entry point `hyb2-sam-composition`).

Flags here mirror bin/hyb2's getopts string ("i:d:o:v:m:h:a:b:q:x:y:l:j:e:r:")
so the eventual replacement is a drop-in for users. `main()` is a stub
(raises NotImplementedError) until Tier 3 begins; `test_cli.py` pins the
intended flag surface in the meantime.
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
        default="vienna",
        choices=["vienna", "unafold", "couplefold"],
        help="Folding backend (default=couplefold once integrated; vienna for now)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.in_file is None:
        parser.print_help()
        return 0

    raise NotImplementedError(
        "hyb2-py does not run the pipeline yet -- stages are being ported "
        "incrementally into hyb2.stages. Use bin/hyb2 for now."
    )


if __name__ == "__main__":
    sys.exit(main())
