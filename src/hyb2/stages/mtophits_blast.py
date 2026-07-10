""" Port of mtophits_blast.awk

Keeps all rows for a read ID as long as column 11 is consistent with the
first row seen for that ID. Rows with the same ID but a different column 11
value get dropped

"""

from __future__ import annotations

import sys

"""
Original AWK tracks the comparison e-value in one global (not per-ID)
variable, so an intervening different-ID row overwrites it. A later row
that would otherwise match its own ID's first e-value gets dropped.
"""

def deduplicate_by_second_fragement_start(text: str) -> str:
    counts: dict[str, int] = {}
    current_value = None
    out: list[str] = []

    for line in text.splitlines():
        columns = line.split("\t")
        read_id = columns[0]
        second_fragment = columns[10]

        counts[read_id] = counts.get(read_id, 0) + 1

        if counts[read_id] == 1:
            current_value = second_fragment
            out.append(line)
        elif second_fragment == current_value:
            out.append(line)

    return "\n".join(out) + "\n"

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: mtophits_blast <hyb_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        sys.stdout.write(deduplicate_by_second_fragement_start(f.read()))
    return 0


if __name__ == "__main__":
    sys.exit(main())