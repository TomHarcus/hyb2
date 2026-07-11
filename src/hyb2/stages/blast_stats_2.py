""" Port of blast_stats_2.awk

Generates a summary table with headers:
biotype   collapsed   all

"""

import sys

def blast_stats_2(text: str) -> str:

    comp = {}
    dec = {}

    for line in text.splitlines():
        columns = line.split("\t")
        
        cnt = columns[0].split("_")
        a = columns[1].split("|")
        b = a[-1].split("_")
        name = b[-1]

        comp[name] = comp.get(name, 0) + 1
        dec[name] = dec.get(name, 0) + int(cnt[1])

    out = ["biotype\tcollapsed\tall"]

    for name in sorted(comp, key=lambda x: comp[x], reverse=True):
        out.append(f"{name}\t{comp[name]}\t{dec[name]}")

    return "\n".join(out) + "\n"

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: blast_stats_2 <hyb_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        sys.stdout.write(blast_stats_2(f.read()))
    return 0


if __name__ == "__main__":
    sys.exit(main())