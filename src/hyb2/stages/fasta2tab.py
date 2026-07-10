"""Port of bin/fasta2tab.awk.

Converts multi-line FASTA to one tab-separated "header\\tsequence" record
per line
"""

from __future__ import annotations

import sys


def fasta_to_tab(text: str) -> str:
    out: list[str] = []
    first_record = True
    for line in text.splitlines():
        if line.startswith("#") or line == "":
            continue
        if line.startswith(">"):
            header = line[1:]
            out.append(f"{header}\t" if first_record else f"\n{header}\t")
            first_record = False
        else:
            out.append(line)
    out.append("\n")
    return "".join(out)


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
