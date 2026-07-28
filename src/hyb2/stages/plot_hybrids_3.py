""" Port of plot_hybrids_3.awk

Returns the awk's output text: a '# hybrids...' 
comment line then 'i\tj\tcount' rows

"""

def plot_hybrids_3(lines, gene_1, gene_2, bin_size=10, threshold_num=0):

    c, d = {}, {}

    for line in lines:

        f = line.rstrip("\n").split("\t")

        if len(f) < 14:
            continue

        a = bin_size * (int(f[6]) // bin_size) # $7 arm1 start
        b = bin_size * (int(f[7]) // bin_size) # $8 arm1 end

        e = bin_size * (int(f[12]) // bin_size) # $13 arm2 start
        g = bin_size * (int(f[13]) // bin_size) # $14 arm2 end

        if f[3] == gene_1 and f[9] == gene_2:
            for i in range(a, b + 1, bin_size):
                for j in range(e, g + 1, bin_size):
                    c[(i,j)] = c.get((i,j), 0) + 1

        if f[3] == gene_2 and f[9] == gene_1 and gene_1 != gene_2:
            for i in range(a, b + 1, bin_size):
                for j in range(e, g + 1, bin_size):
                    d[(i, j)] = d.get((i,j), 0) + 1

    out = [f"# hybrids between {gene_1} and {gene_2}"]
    out += [f"{i}\t{j}\t{n}" for (i,j), n in c.items() if n > threshold_num]

    if gene_1 != gene_2:
        out.append(f"# hybrids between {gene_2} and {gene_1}")
        out += [f"{i}\t{j}\t{n}" for (i, j), n in d.items() if n > threshold_num]

    return "\n".join(out) + "\n"

def swap_gene1_to_arm1(lines, gene_1):
    """Put gene_1 in arm1"""
    out = []

    for line in lines:
        c = line.rstrip("\n").split("\t")
        if len(c) < 15:
            continue

        if c[3] == gene_1:
            out.append("\t".join(c[:15]))

        if c[9] == gene_1:
            out.append("\t".join([c[0], c[1], c[2],
                                  c[9], c[10], c[11], c[12], c[13], c[14],
                                  c[3], c[4], c[5], c[6], c[7], c[8]]))

    return out