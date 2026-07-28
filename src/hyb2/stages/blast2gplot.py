""" Port of blast2gplot.pl

"""

from collections import defaultdict

def blast2gplot(exp="exp1", n_genes=10, blast_file="input.blast",
                ref_blast_file="", gene_lengths_file="human_OH_3_gene_lengths.txt"):

    if not ref_blast_file:
        ref_blast_file = blast_file

    try:
        with open(ref_blast_file, "r") as f:
            gi = {}
            lines = f.read().splitlines()

            for line in lines:
                fields = line.split("\t")

                if int(fields[9]) > int(fields[8]):
                    gi[fields[1]] = gi.get(fields[1], 0) + 1

    except OSError as e:
        raise ValueError(f"Cannot open file: {ref_blast_file}") from e


    sorted_gi = sorted(gi, key=gi.get, reverse=True)
    top_gi = {}

    for i in range(n_genes):
        top_gi[sorted_gi[i]] = top_gi.get(sorted_gi[i], 0) + 1

    try:
        with open(gene_lengths_file, "r") as f:
            len_map = {}
            lines = f.read().splitlines()

            for line in lines:
                fields = line.split("\t")
                curr_gi = fields[0]
                len_map[curr_gi] = int(fields[1])
                

    except OSError as e:
        raise ValueError(f"Cannot open file: {gene_lengths_file}") from e


    try:
        with open(blast_file, "r") as f:
            data = defaultdict(lambda: defaultdict(int))
            lines = f.read().splitlines()

            for line in lines:
                fields = line.split("\t")

                curr_gi = fields[1]

                if (top_gi.get(curr_gi) and int(fields[9]) > int(fields[8])):
                    for i in range(int(fields[8]), int(fields[9])+1):
                        data[curr_gi][i] += 1

    except OSError as e:
        raise ValueError(f"Cannot open file: {blast_file}") from e


    for i in range(n_genes):

        curr_gene = sorted_gi[i]
        curr_filename = f"{exp}_{curr_gene}.gplot"
        curr_gplot_tag = f"#{exp}_{curr_gene}_{len_map[curr_gene]}\n"

        with open(curr_filename, "w") as f:
            f.write(curr_gplot_tag)
            for j in range(1, len_map[curr_gene]+1):
                val = data[curr_gene][j]

                f.write(f"{j}\t{val}\n")

    return None

