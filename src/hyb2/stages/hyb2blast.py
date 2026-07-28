""" Port of hyb2blast.awk

Emits two 12 column BLAST rose per hyb line

"""

def hyb2blast(lines):

    out = []

    for line in lines:
        columns = line.split("\t")

        out.append("\t".join(
            columns[0], columns[3], ".", ".", ".", ".",
            columns[4], columns[5], columns[6], columns[7], 
            columns[8], "."
        ))

        out.append("\t".join(
                    columns[0], columns[9], ".", ".", ".", ".",
                    columns[10], columns[11], columns[12], columns[13], 
                    columns[14], "."
                ))


    return "\t".join(out) + "\n"


