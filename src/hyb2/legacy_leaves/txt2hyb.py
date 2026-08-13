""" Port of txt2hyb.awk

Takes in tab file and hyb file

Looks up read ID from tab file and substitutes into hyb files second column

"""


import sys

def add_sequences_to_hyb(tab_file: str, hyb_file: str) -> str:
    seq_lookup = {}

    for line in tab_file.splitlines():
        columns = line.split("\t")

        seq_lookup[columns[0]] = columns[1]

    out: list[str] = []
    for line in hyb_file.splitlines():
        columns = line.split("\t")
        read_id = columns[0]

        if read_id in seq_lookup:
            columns[1] = seq_lookup[read_id]
            # trailing tab before the new line, same behaviour as the legacy awk
            out.append("\t".join(columns) + "\t")

    return "\n".join(out) + "\n"

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print("usage: txt2hyb <tab_file> <hyb_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f1, open(argv[1]) as f2:
        sys.stdout.write(add_sequences_to_hyb(f1.read(), f2.read()))
    return 0


if __name__ == "__main__":
    sys.exit(main())