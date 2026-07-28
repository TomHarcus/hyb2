""" Port of plot_viewpoint

"""

import re, subprocess
from pathlib import Path
from hyb2.stages.hyb2blast import hyb2blast
from hyb2.stages.blast2gplot import blast2gplot
from hyb2 import config

def plot_viewpoint(in_hyb, db_1, gene_1, gene_2):

    # legacy DB_2 is unreachable: always == DB_1

    lines = open(db_1).read().splitlines()
    len_1 = None

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

    out_blast = in_hyb.replace(".hyb", f"_{gene_1}.blast")
    Path(out_blast).write_text(hyb2blast(filtered))

    blast2gplot(exp=in_hyb.replace(".hyb", ""), n_genes=1, 
                ref_blast_file=f"{gene_1}_ref.blast",
                blast_file=out_blast,
                gene_lengths_file=f"{gene_1}.length.txt")


    subprocess.run(["Rscript", config.rscript("viewpoint_graph.R"),
                    f"{in_hyb.replace(".hyb", "")}_{gene_1}.gplot"],
                    check=True)


    """
    still got all the branches to add
    """









