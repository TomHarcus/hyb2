"""CLI skeleton reserved for the top-level `hyb2` mega-orchestrator port.

Flags here mirror bin/hyb2's getopts string ("i:d:o:v:m:h:a:b:q:x:y:l:j:e:r:")

+ the hy2_fold args specifically cplfold customisation

"""


import argparse
import sys
import yaml
import logging

from hyb2.config import CPL_DEFAULTS


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
    parser.add_argument("-V", "--verbose", action="store_true", help="show detailed stage output")
    parser.add_argument("--config", default=None, metavar="RUN.YML", help="YAML config of args, CLI flags override it")
    parser.add_argument("-i", "--input", dest="input", metavar="INPUT", help="Input fastq/SAM/hyb file")
    parser.add_argument("-d", "--reference", dest="reference", metavar="FASTA", help="Reference fasta used for mapping")
    parser.add_argument("-o", "--output-id", dest="output_id", metavar="OUTPUT_ID", help="Output ID")
    parser.add_argument("-v", "--blast-threshold", dest="blast_threshold", type=float, default=0.1, help="BLAST threshold (default=0.1)")
    parser.add_argument("-m", "--max-overlap", dest="max_overlap", type=int, default=4, help="Maximum overlap (default=4)")
    parser.add_argument("-h", "--max-hits", dest="max_hits", type=int, default=10, help="Maximum hits per sequence (default=10)")
    parser.add_argument("-a", "--gene-1", dest="gene_1", metavar="GENE_ID", help="Gene ID of interest")
    parser.add_argument("-b", "--gene-2", dest="gene_2", metavar="GENE_ID", help="Second gene ID of interest")
    parser.add_argument("-q", "--heatmap-quantile", dest="heatmap_quantile", type=float, default=0.95, help="Heatmap upper quantile (default=0.95)")
    parser.add_argument("-x", "--x-start", dest="x_start", type=int, help="Start coordinate of 1st strand/gene")
    parser.add_argument("-y", "--y-start", dest="y_start", type=int, help="Start coordinate of 2nd strand/gene")
    parser.add_argument("-l", "--length", dest="length", type=int, help="Length of fragments")
    parser.add_argument("-j", "--varna-jar", dest="varna_jar", metavar="VARNAcmd.jar", help="Path to VARNAcmd.jar")
    parser.add_argument(
        "-e", "--calc-energy", # for the add_dG_hyb2_2 port
        dest="calc_energy",
        type=int,
        default=0,
        help="Calculate folding energy: 1 on, 0 off (default=0)",
    ) 
    
    parser.add_argument(
        "-r", "--fold-backend",
        dest="fold_backend",
        default="cplfold",
        choices=["cplfold", "vienna", "unafold", "0", "1"]
    )
    parser.add_argument("-0", "--interactive", dest="interactive", default=None, help="1 to launch the interactive VARNA GUI")
    parser.add_argument("--alpha", dest="alpha", type=float, default=CPL_DEFAULTS["alpha"], help="cplfold bonus weight")
    parser.add_argument("--beta", dest="beta", type=float, default=CPL_DEFAULTS["beta"], help="cplfold bonus weight")
    parser.add_argument("--normalize", dest="normalize", choices=["raw", "log"], default=CPL_DEFAULTS["normalize"], help="cplfold bonus normalization")
    parser.add_argument("--beam-size", dest="beam_size", type=int, default=CPL_DEFAULTS["beam_size"], help="cplfold beam size")
    parser.add_argument("--energy-delta", dest="energy_delta", type=float, default=CPL_DEFAULTS["energy_delta"], help="cplfold energy delta")
    parser.add_argument("--max-phase1", dest="max_phase1", type=int, default=CPL_DEFAULTS["max_phase1"], help="cplfold max phase 1")
    parser.add_argument("--max-phase2", dest="max_phase2", type=int, default=CPL_DEFAULTS["max_phase2"], help="cplfold max phase 2")
    parser.add_argument("--energy-model", dest="energy_model", choices=["DP09", "DP03", "CC06", "CC09", "RE"], default=CPL_DEFAULTS["energy_model"], help="cplfold energy model")

    return parser

from hyb2.logsetup import configure_logging

def main(argv: list[str] | None = None) -> int:

    argv = sys.argv[1:] if argv is None else argv

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=None)
    pre_args, _ = pre.parse_known_args(argv)

    parser = build_parser()

    # parse the yaml file if present
    if pre_args.config:
        with open(pre_args.config) as fh:
            cfg = yaml.safe_load(fh) or {}

        valid = {a.dest for a in parser._actions if a.dest not in ("help", "config")}
        unknown = set(cfg) - valid

        if unknown:
            parser.error(f"unknown config keys: {sorted(unknown)}")
        parser.set_defaults(**cfg)

    args = parser.parse_args(argv)

    configure_logging(args.verbose)

    from hyb2.pipelines import hyb2

    if args.input is None:
        return hyb2.print_help()

    hyb2.run(args.input, args.reference, args.output_id,
             hmax=args.max_hits, blast_threshold=args.blast_threshold, max_overlap=args.max_overlap,
             gene_1=args.gene_1, gene_2=args.gene_2, limit=args.heatmap_quantile,
             x_coord=args.x_start, y_coord=args.y_start, length=args.length,
             varna=args.varna_jar, fold=args.fold_backend,
             interactive=(str(args.interactive) == "1"),
             alpha=args.alpha, beta=args.beta, normalize=args.normalize,
             beam_size=args.beam_size, energy_delta=args.energy_delta,
             max_phase1=args.max_phase1, max_phase2=args.max_phase2,
             energy_model=args.energy_model)

    return 0


if __name__ == "__main__":
    sys.exit(main())
