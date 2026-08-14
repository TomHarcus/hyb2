""" Port of plot_viewpoint

"""

import re, subprocess, os, argparse, sys, contextlib
from pathlib import Path
from hyb2.viewpoint.hyb2blast import hyb2blast
from hyb2.viewpoint.blast2gplot import blast2gplot
from hyb2.tools import config
from hyb2.tools.logsetup import is_quiet
from hyb2.tools.ui import spinner

import logging

log = logging.getLogger(__name__)

def plot_viewpoint(in_hyb, db_1, gene_1, gene_2):

    quiet = is_quiet()

    # legacy DB_2 is unreachable: always == DB_1

    lines = open(db_1).read().splitlines()
    len_1 = None

    stem = in_hyb.replace(".hyb", "")

    for i, line in enumerate(lines):
        if re.search(gene_1, line):
            len_1 = len(lines[i+1])
            break

    with open(f"{gene_1}.length.txt", "w") as f:
        f.write(f"{gene_1}\t{len_1}\n")

    with open(in_hyb) as f:
        blast_rows = hyb2blast(f).splitlines()

    ref = next((l for l in blast_rows if re.search(gene_1, l)), "")
    Path(f"{gene_1}_ref.blast").write_text(ref + "\n")


    with open(in_hyb) as f:
        lines = f.read().splitlines()
        filtered = []
        for line in lines:
            columns = line.split("\t")
            if re.search(gene_1, columns[3]) and re.search(gene_1, columns[9]):
                filtered.append(line)

    out_blast = f"{stem}_{gene_1}.blast"
    Path(out_blast).write_text(hyb2blast(filtered))

    

    blast2gplot(exp=stem, n_genes=1, 
                ref_blast_file=f"{gene_1}_ref.blast",
                blast_file=out_blast,
                gene_lengths_file=f"{gene_1}.length.txt")

    with spinner("rendering viewpoint graph ") if quiet else contextlib.nullcontext():
        if not quiet:
            print("rendering viewpoint graph")
        subprocess.run(["Rscript", config.rscript("viewpoint_graph.R"),
                        f"{stem}_{gene_1}.gplot"],
                        stderr=subprocess.DEVNULL if quiet else None,
                        check=True)

    homodimers = in_hyb.replace(".hyb", "_homodimers.hyb")

    if os.path.isfile(homodimers) and os.path.getsize(homodimers) > 0:
        with open(homodimers) as f:
            lines = f.read().splitlines()
            filtered = []
            for line in lines:
                columns = line.split("\t")
                if re.search(gene_1, columns[3]) and re.search(gene_1, columns[9]):
                    filtered.append(line)
        
        out_homodimers_blast = f"{stem}_{gene_1}_homodimers.blast"
    
        Path(out_homodimers_blast).write_text(hyb2blast(filtered))

        blast2gplot(exp=f"{stem}_homodimers", n_genes=1, 
                    ref_blast_file=f"{gene_1}_ref.blast",
                    blast_file=out_homodimers_blast,
                    gene_lengths_file=f"{gene_1}.length.txt")

        with spinner("rendering viewpoint graph (homodimers)") if quiet else contextlib.nullcontext():
            if not quiet:
                print("rendering viewpoint graph (homodimers)")
            subprocess.run(["Rscript", config.rscript("viewpoint_graph.R"),
                            f"{stem}_homodimers_{gene_1}.gplot"],
                            stderr=subprocess.DEVNULL if quiet else None,
                            check=True)
        log.info(f"Viewpoint graph of {gene_1} saved")

    if gene_2:

        with open(in_hyb) as f:
            filtered = [l for l in f if re.search(gene_1, l) and re.search(gene_2, l)]
        
        blast_g1_g2 = f"{stem}_{gene_1}_{gene_2}.blast"
        Path(blast_g1_g2).write_text(hyb2blast(filtered))

        blast2gplot(exp=f"{stem}_{gene_2}", n_genes=1, 
                    ref_blast_file=f"{gene_1}_ref.blast",
                    blast_file=blast_g1_g2,
                    gene_lengths_file=f"{gene_1}.length.txt")

        lines = open(db_1).read().splitlines()
        len_2 = None

        for i, line in enumerate(lines):
            if re.search(gene_2, line):
                len_2 = len(lines[i+1])
                break

        with open(f"{gene_2}.length.txt", "w") as f:
            f.write(f"{gene_2}\t{len_2}\n")


        with open(in_hyb) as f:
            blast_rows = hyb2blast(f).splitlines()
        
        ref = next((l for l in blast_rows if re.search(gene_2, l)), "")
        Path(f"{gene_2}_ref.blast").write_text(ref + "\n")

        blast2gplot(exp=f"{stem}_{gene_1}", n_genes=1, 
                    ref_blast_file=f"{gene_2}_ref.blast",
                    blast_file=blast_g1_g2,
                    gene_lengths_file=f"{gene_2}.length.txt")

        with spinner(f"rendering viewpoint graph ({gene_2}) ") if quiet else contextlib.nullcontext():
            if not quiet:
                print(f"rendering viewpoint graph ({gene_2})")
            subprocess.run(["Rscript", config.rscript("viewpoint_graph.R"),
                            f"{stem}_{gene_2}_{gene_1}.gplot"],
                            stderr=subprocess.DEVNULL if quiet else None,
                            check=True)

        with spinner(f"rendering viewpoint graph ({gene_1}) ") if quiet else contextlib.nullcontext():
            if not quiet:
                print(f"rendering viewpoint graph ({gene_1})")
            subprocess.run(["Rscript", config.rscript("viewpoint_graph.R"),
                            f"{stem}_{gene_1}_{gene_2}.gplot"],
                            stderr=subprocess.DEVNULL if quiet else None,
                            check=True)
        log.info(f"Viewpoint graph of {gene_2} saved")
        
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="plot-viewpoint", add_help=False)
    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-V", "--verbose", action="store_true", help="show detailed output")
    p.add_argument("-i", dest="in_hyb", required=True, metavar="INPUT.HYB")
    p.add_argument("-d", dest="db_1", required=True, metavar="REFERENCE.FASTA")
    p.add_argument("-a", dest="gene_1", required=True, metavar="GENE_1")
    p.add_argument("-b", dest="gene_2", default=None, metavar="GENE_2")

    return p

from hyb2.tools.logsetup import configure_logging
def main(argv=None):
    a = build_parser().parse_args(argv)
    configure_logging(a.verbose)
    plot_viewpoint(a.in_hyb, a.db_1, a.gene_1, a.gene_2)
    return 0

if __name__ == "__main__":
    sys.exit(main())








