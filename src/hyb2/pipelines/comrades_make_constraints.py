""" Port of comradesMakeConstraints_2

Folds each chimeric read pairs arm sequences with RNAcofold, tallies basepair frequencies across 
all chimeras in a given reference window, and outputs the top-N most frequent basepairs as folding 
constraints for that region

"""

import argparse, subprocess, sys, os, shutil
from hyb2.stages.fasta2tab import fasta_to_tab
from hyb2.stages.hyb2fasta_bits_allRNAs import hyb2fasta_bits_allRNAs
from hyb2.stages.ct2bps_2 import ct2bps_2
from hyb2.stages.histogram import histogram
from hyb2.stages.bp2hyb import bp2hyb
from hyb2.stages.hyb2constraints import hyb2constraints
from hyb2 import config
import logging

log = logging.getLogger(__name__)

def run(in_hyb, ref_fasta, begin, end, *, num_constraints=75, fold="vienna",
        vienna_bin=None):

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


    with open(ref_tab) as tab, open(in_hyb) as hyb, open(bit1, "w") as o1, open(bit2, "w") as o2:
        for bit1_rec, bit2_rec in hyb2fasta_bits_allRNAs(tab, hyb):
            o1.write(bit1_rec)
            o2.write(bit2_rec)

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

def _fold_vienna(bit1, bit2, ct, vienna_bin):
    b = (vienna_bin.rstrip("/") + "/") if vienna_bin else ""

    with open(bit1) as f1, open(bit2) as f2:
        pasted = "".join(a.rstrip("\n") + "&" + b2 for a, b2 in zip(f1, f2))

        # run vienna cofold
        cofold = subprocess.run([b + "RNAcofold", "--noconv", "--noPS"],
                                input=pasted, capture_output=True, text=True,
                                check=True).stdout
        
        vienna = cofold.replace("&>", "-").replace("&", "")

        ctdata = subprocess.run([b + "b2ct"], input=vienna,
                                capture_output=True, text=True, check=True).stdout
        
        ctdata = ctdata.replace("ENERGY =", "dG =")

        with open(ct, "w") as fout:
            fout.write(ctdata)

def _fold_unafold(bit1, bit2, ct):

    # run unafold
    subprocess.run(["hybrid-min", bit1, bit2],
                   capture_output=True, text=True, check=True,
                   env={**os.environ, "UNAFOLDDAT": str(config.UNAFOLD_DIR)})

    result = f"{bit1}-{bit2}.ct"

    if result != ct:
        shutil.move(result, ct)

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="comrades-make-constraints",
        description="create constraints",
        add_help=False
    )

    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-i", dest="in_hyb", required=True, metavar="INPUT.HYB", help="Input HYB (required)")
    p.add_argument("-f", dest="ref_fasta", required=True, metavar="INPUT.FASTA", help="Input FASTA (required)")
    p.add_argument("-b", dest="begin", required=True, type=int, help="start coordinate")
    p.add_argument("-e", dest="end", required=True, type=int, help="end coordinate")
    p.add_argument("-n", dest="num_constraints", type=int, default=75, help="number of folding constraints (default=75)")
    p.add_argument("-r", dest="fold", type=str, default="vienna", help="folder to be used")

    return p

def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
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