""" Port of filter_homopolymers.awk

Filters out sequences that contain sequences of 15 or more consecutive
identical nucleotides (AAAAAAAAAAAAAAA, CCCCCCCCCCCCCCC, etc)

"""

import sys
import re

def filter_homopolymers(text: str) -> str:
    out: list[str] = []

    filter_length = 15

    for line in text.splitlines():
        # regex to filter consecutive sequence runs of A/C/G/T (upper or lower case not mixed) not U
        match = re.search(rf'([AGCTagct])\1{{{filter_length-1}}}', line)
        if match:
            continue
        
        out.append(line)

    return "".join(f"{line}\n" for line in out)

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: filter_homopolymers <fasta_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        sys.stdout.write(filter_homopolymers(f.read()))
    return 0

if __name__ == "__main__":
    sys.exit(main())