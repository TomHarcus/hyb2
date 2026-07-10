""" Port of hybrid_stats_2.awk

Generates a summary table with headers:
type1:type2   collapsed   all

"""

from __future__ import annotations

import sys

"""
SORTED is hardcoded to 0 in the original script and never reassigned
"""
 

def hybrid_stats_2(text: str) -> str:

    comp = {}
    dec = {}

    sort_pairs = False

    for line in text.splitlines():
        columns = line.split("\t")

        cnt = columns[0].split("_")
        n = columns[3].split("|")
        b = n[-1].split("_")
        nm1 = b[-1]

        m = columns[9].split("|")
        y = m[-1].split("_")
        nm2 = y[-1]

        if sort_pairs:
            if nm1 < nm2:
                name = f"{nm1}:{nm2}"
            else:
                name = f"{nm2}:{nm1}"
        else:
            name = f"{nm1}:{nm2}"

        comp[name] = comp.get(name, 0) + 1
        dec[name] = dec.get(name, 0) + int(cnt[1])

    out = ["type1:type2\tcollapsed\tall"]

    for name in sorted(comp, key=lambda x: comp[x], reverse=True):
        out.append(f"{name}\t{comp[name]}\t{dec[name]}")

    return "\n".join(out) + "\n"

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: hybrid_stats_2 <hyb_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        sys.stdout.write(hybrid_stats_2(f.read()))
    return 0


if __name__ == "__main__":
    sys.exit(main())