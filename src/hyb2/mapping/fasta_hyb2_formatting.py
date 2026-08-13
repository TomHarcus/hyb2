""" Port of fasta_hyb2_formatting

"""

def fasta_hyb2_formatting(text):

    lines = text.splitlines()

    for i, line in enumerate(lines):

        line = line.replace("_", "-")
        line = line.replace("protein-coding", "mRNA")
        line = line.replace("|", "_")
        line = line.replace("__", "_unknown_")

        lines[i] = line

    buf = ""

    for i, line in enumerate(lines):
        if i == 0:
            buf += line + "\n"

        elif line.startswith(">"):
            buf += "\n" + line + "\n"

        else:
            buf += line

    out, skip = [], False

    for line in buf.splitlines():
        if skip:
            skip = False
            continue
        if "unknown" in line or "pseudogene" in line:
            skip = True
            continue

        out.append(line)

    return "\n".join(out)
