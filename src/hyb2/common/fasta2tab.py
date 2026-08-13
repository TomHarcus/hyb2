"""Port of bin/fasta2tab.awk.

Converts multi-line FASTA to one tab-separated "header\\tsequence" record
per line
"""


import sys


def fasta_to_tab_lines(lines):
    header, seq = None, []

    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#") or line == "":
            continue
        if line.startswith(">"):
            if header is not None:
                yield header + "\t" + "".join(seq)
            header, seq = line[1:], []
        else:
            seq.append(line)

    if header is not None:
        yield header + "\t" + "".join(seq)

def fasta_to_tab(text: str) -> str:
    return "\n".join(fasta_to_tab_lines(text.splitlines())) + "\n"


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: fasta2tab <fasta_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        sys.stdout.write(fasta_to_tab(f.read()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
