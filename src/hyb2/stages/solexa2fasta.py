""" Port of solexa2fasta.awk

Converts fastq file to fasta file

"""

from __future__ import annotations

import sys

def solexa_to_fasta(text: str) -> str:
    lines = text.splitlines()
    out: list[str] = []
    i = 0

    while i < len(lines):
        if lines[i].startswith("@"):
            # possible bug in solexa2fasta.awk where it keeps the @ symbol
            out.append(">@" + lines[i][1:])
            out.append(lines[i+1])
            i += 4
        else:
            i += 1
    return "\n".join(out) + "\n"

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: solexa2fasta <fastq_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        sys.stdout.write(solexa_to_fasta(f.read()))
    return 0

if __name__ == "__main__":
    sys.exit(main())