""" Port of hyb2_fold


"""

import os, shutil, re, glob, math

from hyb2.config import varna_jar

def run(in_hyb, GENE_1, GENE_2, FASTA_1, x_coord, y_coord, length, VARNA, interactive, FOLD, vienna_bin=None):

    if FOLD is None:
        FOLD = 1

    if VARNA is None:
        VARNA = varna_jar()

    fasta_1 = FASTA_1
    fasta_2 = None

    shutil.copy(fasta_1, "./")
    if fasta_2:
        shutil.copy(fasta_2, "./")
    else:
        fasta_2 = fasta_1

    fasta_1 = os.path.basename(fasta_1)
    fasta_2 = os.path.basename(fasta_2)

    if GENE_2 is None and y_coord is not None and (x_coord > y_coord):
        raise ValueError("Intramolecular Folding x_coord must be lower than y_coord")
    
    if y_coord is None and GENE_2 is None:
        X1 = x_coord
        X2 = x_coord + length - 1
        out_file = in_hyb.replace(".hyb", f"_{GENE_1}_{X1}-{X2}.hyb")
        out_fasta = f"{GENE_1}_{X1}-{X2}.fasta"
        span = length

        _transform(in_hyb, out_file, GENE_1, x_coord, X1, X2)
        _fasta_extraction(fasta_1, GENE_1, X1, length, out_fasta)
        _fold(out_file, out_fasta, span, FOLD, vienna_bin, basepair_scores=None)
        _postfold(in_hyb, out_file, out_fasta, span, x_coord, y_coord, length,
              VARNA, interactive)

    


def _transform(in_hyb, out_file, GENE, x_coord, X1, X2):

    with open(in_hyb) as fin, open(out_file, "w") as fout:
        for line in fin:
            c = line.rstrip("\n").split("\t")

            if (re.search(GENE, c[3]) and re.search(GENE, c[9])
                and int(c[6]) >= X1 and int(c[7]) <= X2
                and int(c[12]) >= X1 and int(c[13]) <= X2):

                out = [
                    c[0], c[1], c[2], c[3], c[4], c[5],
                    str(int(c[6]) - x_coord + 1),
                    str(int(c[7]) - x_coord + 1),
                    c[8], c[9], c[10],
                    "12",
                    str(int(c[12]) - x_coord + 1),
                    str(int(c[13]) - x_coord + 1)
                ]

                fout.write("\t".join(out) + "\n")

def _fasta_extraction(fasta_1, GENE, X1, length, out_fasta):
    
    records, cur = {}, None

    for l in open(fasta_1).read().splitlines():
        if l.startswith(">"):
            cur = l
            records[cur] = ""

        elif cur is not None:
            records[cur] += l.strip()
    
    header, seq = next((h, s) for h, s in records.items() if GENE in h)
    name = header.split()[0]
    fragment = seq[X1 - 1: X1 - 1 + length]

    with open(out_fasta, "w") as f:
        f.write(f"{name}\n{fragment}\n")

def _fold(out_file, out_fasta, span, fold, vienna_bin, basepair_scores=None):
    
    from hyb2.pipelines import comrades_make_constraints, comrades_fold

    comrades_make_constraints.run(out_file, out_fasta, 1, span,
                                  fold=fold, vienna_bin=vienna_bin)
    
    if fold == "cplfold":
        bp_scores = out_file.replace(".hyb", ".basepair_scores.txt")

        comrades_fold.run("unused", out_fasta, fold="cplfold",
                          basepair_scores=bp_scores, begin=1, end=span,
                          vienna_bin=vienna_bin)
        
    else:
        constraints = out_file.replace(".hyb", f".1-{span}_folding_constraints.txt")

        comrades_fold.run(constraints, out_fasta, fold=fold, vienna_bin=vienna_bin)

  
def _postfold(in_hyb, out_file, out_fasta, span, x_coord, y_coord, length, VARNA, interactive):
    
    coords = out_file.replace(".hyb", "_structure_coordinates.txt")
    cts = sorted(glob.glob(out_fasta.replace(".fasta", "") + "*.ct"))

    with open(coords, "w") as f:
        for ct in cts:
            f.write(f"{ct}\t0\t{span}\t{span}\n")

    bp_scores = out_file.replace(".hyb", ".basepair_scores.txt")

    from hyb2.stages.bp_score import bp_score

    with open(coords) as c:
        bp_score(c, bp_scores, in_hyb)

    scores = sorted(glob.glob(in_hyb.replace(".hyb", "") + "__"
                              + out_fasta.replace(".fasta", "") + "*VARNA_scores.txt"))
    
    ranked = sorted(scores,
                    key=lambda fn: sum(float(x) for x in open(fn) if x.strip()),
                    reverse=True)
    
    name = ranked[0]
    vname = name.split("__", 1)[1]

    log2 = name.replace("_scores.txt", "_log2scores.txt")

    with open(name) as fin, open(log2, "w") as fout:
        for x in fin:
            v = math.log2(float(x) + 1) if x.strip() else 0
            fout.write(f"{0 if v < 0 else v}\n")

    ct_top = vname.replace(".VARNA_scores.txt", ".ct")
    from hyb2.pipelines.plot_VARNA import plot_VARNA

    plot_VARNA(ct_top, log2, VARNA, x_coord=x_coord, y_coord=y_coord,
               length=length, interactive=interactive)
    
    print("randomized parallel RNA folding (cluster): qsub comradesFold2 ... -s 1")
    print("assign scores: comradesScore -i <basepair_scores> -f <out_fasta>")


    

