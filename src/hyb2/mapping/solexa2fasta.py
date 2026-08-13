""" Port of solexa2fasta.awk

Converts fastq file to fasta file

"""


import sys

def solexa_to_fasta_lines(lines):
    it = iter(lines)
    for line in it:
        line = line.rstrip("\n")
        if line.startswith("@"):
            seq = next(it).rstrip("\n")
            yield ">@" + line[1:]
            yield seq
            next(it, None)
            next(it, None)

def solexa_to_fasta(text: str) -> str:
    return "\n".join(solexa_to_fasta_lines(text.splitlines())) + "\n"

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