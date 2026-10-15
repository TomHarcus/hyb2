""" Port of comradesMakeConstraints_2

Folds each chimeric read pairs arm sequences with RNAcofold, tallies basepair frequencies across 
all chimeras in a given reference window, and outputs the top-N most frequent basepairs as folding 
constraints for that region

"""

import argparse, subprocess, sys, os, shutil
from hyb2.common.fasta2tab import fasta_to_tab
from hyb2.folding.hyb2fasta_bits_allRNAs import hyb2fasta_bits_allRNAs
from hyb2.folding.ct2bps_2 import ct2bps_2
from hyb2.common.histogram import histogram
from hyb2.folding.bp2hyb import bp2hyb
from hyb2.folding.hyb2constraints import hyb2constraints
from hyb2.tools import config
from hyb2.tools.ui import spinner
import logging
import concurrent.futures

log = logging.getLogger(__name__)

def _validate_single_rna(in_hyb):
    """Coordinate-only constraints require one RNA coordinate system.

    hyb2_fold has already selected and remapped two-RNA interactions to a
    synthetic single RNA. Raw intergenic input must instead select a gene pair.
    """
    gene = None
    with open(in_hyb) as source:
        for number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                raise ValueError(f"Invalid HYB row {number}: expected at least 14 columns")
            if gene is None:
                gene = fields[3]
            if fields[3] != gene or fields[9] != gene:
                raise ValueError(
                    "Mixed-gene HYB input cannot be scored as one RNA: "
                    "supply gene1/gene2 (--gene-1/--gene-2) for gene-pair scores, "
                    "or use hyb2 fold -a RNA_A -b RNA_B with region coordinates "
                    "for structure prediction."
                )


def run(in_hyb, ref_fasta, begin=None, end=None, *, num_constraints=75, fold="vienna",
        vienna_bin=None, support_matrix=None, gene1=None, gene2=None, output=None):
    """Build single-RNA constraints or explicitly selected gene-pair scores.

    Returns the constraint path in single-RNA mode, or the labelled score-table
    path when gene1/gene2 are supplied. The latter retains original coordinates
    and must not be interpreted as single-sequence folding constraints.
    """
    if gene1 is not None or gene2 is not None:
        if not gene1 or not gene2:
            raise ValueError("gene1 and gene2 must be supplied together")
        if begin is not None or end is not None or support_matrix is not None:
            raise ValueError("gene-pair scoring requires original HYB input without a window or support_matrix")
        from hyb2.folding.gene_pair_scores import score_gene_pair
        return score_gene_pair(in_hyb, ref_fasta, gene1, gene2, fold=fold,
                               vienna_bin=vienna_bin, output=output)
    if begin is None or end is None:
        raise ValueError("single-RNA constraints require begin and end")
    if output is not None:
        raise ValueError("output is only available with gene-pair scoring")
    # Validate before writing intermediates, including when reusing a matrix.
    _validate_single_rna(in_hyb)

    ref_bare = ref_fasta.rsplit(".", 1)[0]
    ref_tab = ref_bare + ".tab"
    ref_frag = f"{ref_bare}_{begin}-{end}.fasta"
    bit1 = in_hyb.replace(".hyb", ".bit_1.fasta")
    bit2 = in_hyb.replace(".hyb", ".bit_2.fasta")
    ct = f"{bit1}-{os.path.basename(bit2)}.ct"
    bp_scores = in_hyb.replace(".hyb", ".basepair_scores.txt")
    frag_scr = in_hyb.replace(".hyb", ".fragment_scores.txt")
    ranked = in_hyb.replace(".hyb", "") + f".{begin}-{end}_ranked_interactions.txt"
    constr = ranked.replace("ranked_interactions.txt", "folding_constraints.txt")

    log.debug(f"Input hybrids file: {in_hyb}")
    log.debug(f"Reference fasta file: {ref_fasta}")
    log.debug(f"Fragment start coordinate: {begin}")
    log.debug(f"Fragment end coordinate: {end}")


    with open(ref_fasta) as fin, open(ref_tab, "w") as fout:
        fout.write(fasta_to_tab(fin.read()))

    with open(ref_tab) as fin, open(ref_frag, "w") as fout:
        for line in fin:
            elements = line.rstrip("\n").split("\t")
            seq = elements[1][begin - 1:end]

            fout.write(f">{elements[0]}\n{seq}\n")


    if support_matrix is None:

        with open(ref_tab) as tab, open(in_hyb) as hyb, open(bit1, "w") as o1, open(bit2, "w") as o2:
            for bit1_rec, bit2_rec in hyb2fasta_bits_allRNAs(tab, hyb):
                o1.write(bit1_rec)
                o2.write(bit2_rec)

        with spinner("building support matrix "):
            if fold in ("vienna", "cplfold"):
                _fold_vienna(bit1, bit2, ct, vienna_bin)

            elif fold == "unafold":
                _fold_unafold(bit1, bit2, ct)
        
            else:
                raise ValueError(f"unknown folder: {fold}")
    
        with open(ct) as fin:
            scores = histogram(ct2bps_2(fin.read()))

        with open(bp_scores, "w") as fout:
            fout.writelines(scores)

    else:
        if not (os.path.isfile(support_matrix) and os.path.getsize(support_matrix) > 0):
            raise FileNotFoundError(f"support matrix not found or empty : {support_matrix}")
        log.info(f"reusing support matrix: {support_matrix}")
        if os.path.abspath(support_matrix) != os.path.abspath(bp_scores):
            shutil.copyfile(support_matrix, bp_scores)

    with open(bp_scores) as fin, open(frag_scr, "w") as fout:
        printed = 0
        for line in fin:
            if printed > 1000:
                break

            elements = line.rstrip("\n").split("\t")

            if begin <= int(elements[0]) <= end and begin <= int(elements[1]) <= end:
                fout.write(line)
                printed += 1

    with open(frag_scr) as fin, open(ranked, "w") as fout:
        fout.writelines(bp2hyb(fin))

    with open(ranked) as fin, open(constr, "w") as fout:
        for count, output in enumerate(hyb2constraints(fin)):
            if count >= num_constraints:
                break

            elements = output.rstrip("\n").split("\t")
            elements[1] = str(int(elements[1]) - (begin - 1))
            elements[2] = str(int(elements[2]) - (begin - 1))

            fout.write("\t".join(elements) + "\n")
        

    return constr


def _fold_vienna(bit1, bit2, ct, vienna_bin, workers=None):
    b = (vienna_bin.rstrip("/") + "/") if vienna_bin else ""
    rnacofold = b + "RNAcofold"

    with open(bit1) as f1, open(bit2) as f2:
        pasted = "".join(a.rstrip("\n") + "&" + b2 for a, b2 in zip(f1, f2))

    lines = pasted.splitlines(keepends=True)
    records = ["".join(lines[i:i + 2]) for i in range(0, len(lines), 2)]

    if workers is None:
        try:
            workers = len(os.sched_getaffinity(0))
        except AttributeError:
            workers = os.cpu_count() or 1

    def _fold_chunk(chunk):
        return subprocess.run([rnacofold, "--noconv", "--noPS"],
                              input=chunk, capture_output=True, text=True,
                              check=True).stdout

    if workers <= 1 or len(records) < 2 * workers:
        cofold = _fold_chunk("".join(records))
    else:
        size = -(-len(records) // workers)

        chunks = ["".join(records[i:i + size]) for i in range(0, len(records), size)]

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            cofold = "".join(ex.map(_fold_chunk, chunks))

    vienna = cofold.replace("&>", "-").replace("&", "")

    ctdata = subprocess.run([b + "b2ct"], input=vienna,
                            capture_output=True, text=True, check=True).stdout

    ctdata = ctdata.replace("ENERGY =", "dG =")

    with open(ct, "w") as fout:
        fout.write(ctdata)

def _fold_unafold(bit1, bit2, ct):

    # hybrid-min derives output names from its input arguments. Run alongside
    # the paired FASTAs so absolute paths do not become embedded in filenames.
    workdir = os.path.dirname(os.path.abspath(bit1))
    first = os.path.basename(bit1)
    second = os.path.relpath(os.path.abspath(bit2), workdir)
    # run unafold
    subprocess.run(["hybrid-min", first, second], cwd=workdir,
                   capture_output=True, text=True, check=True,
                   env={**os.environ, "UNAFOLDDAT": str(config.UNAFOLD_DIR)})

    result = os.path.join(workdir, f"{first}-{second}.ct")

    if os.path.abspath(result) != os.path.abspath(ct):
        shutil.move(result, ct)

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hyb2 scores",
        description="Calculate gene-pair COMRADES scores, or single-sequence constraints with -b/-e",
        add_help=False
    )

    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-i", dest="in_hyb", required=True, metavar="INPUT.HYB", help="Input HYB (required)")
    p.add_argument("-f", dest="ref_fasta", required=True, metavar="INPUT.FASTA", help="Input FASTA (required)")
    p.add_argument("-b", dest="begin", type=int, help="start coordinate (constraint mode)")
    p.add_argument("-e", dest="end", type=int, help="end coordinate (constraint mode)")
    p.add_argument("-1", "--gene-1", dest="gene1", help="exact first reference ID (gene-pair scoring)")
    p.add_argument("-2", "--gene-2", dest="gene2", help="exact second reference ID (gene-pair scoring)")
    p.add_argument("-o", "--output", help="gene-pair score table output path")
    p.add_argument("-n", dest="num_constraints", type=int, default=75, help="number of folding constraints (default=75)")
    p.add_argument("-r", dest="fold", type=str, default="vienna", help="folder to be used")

    return p

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.gene1 is not None or args.gene2 is not None:
        if not args.gene1 or not args.gene2:
            parser.error("--gene-1 and --gene-2 must be supplied together")
        if args.begin is not None or args.end is not None:
            parser.error("gene-pair scoring uses full RNA coordinates; omit -b/-e")
        run(args.in_hyb, args.ref_fasta, gene1=args.gene1, gene2=args.gene2,
            fold=args.fold, output=args.output)
        return 0
    if args.begin is None or args.end is None:
        parser.error("constraint mode requires -b and -e")
    if args.output is not None:
        parser.error("--output is only available with gene-pair scoring")
    run(
        args.in_hyb,
        args.ref_fasta,
        begin=args.begin,
        end=args.end,
        num_constraints=args.num_constraints,
        fold=args.fold
    )

    return 0

if __name__ == "__main__":
    sys.exit(main())
