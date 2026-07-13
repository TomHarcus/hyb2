"""Port of collapse_blast_2.sh (Tier 1 #3).

WHAT IT DOES (plain English): removes duplicate rows, keeping exactly one row per
unique (gene = col 2, mapped sequence = col 13) pair. That is the whole job --
collapsing PCR/optical duplicates so each molecule is counted once.
Verified on real data: 102090 rows in -> 92299 out.

NOTE, despite the legacy "tallying reads to the ID" comment: NO count tallying
survives. The awk recomputes the read-id count, then `cut -f3-20` throws it away,
so every output read ID is unchanged. That path is dead code; this port omits it.

WHY THE CODE BELOW IS CONVOLUTED: collapse's exact ROW ORDER changes the
downstream numbers -- mtophits keeps the first-seen e-value per read id, so a
different order -> different reference sums (reversing the rows shifted Zika
244903 -> 243644). So this faithfully reproduces the legacy script's specific
order and representative row rather than deduping cleanly. If that fragility is
ever judged not worth preserving, all of stages 1-2 could collapse to a clean
"keep first per (gene, seq) in input order" -- but that changes results and needs
a re-baselined golden (a decision for Grzegorz).

The sort key is pinned to code points (locale-independent), which is actually
more reproducible than the legacy shell `sort`, whose order depended on locale.

Parity: fixtures/sam_composition_run/test.collapse.blast
"""


import sys
from typing import Iterable, Iterator


def collapse_blast(lines: Iterable[str]) -> Iterator[str]:
    # Load every row up front: this needs two passes over the data (the sort
    # below, then `rows` is reused whole in stage 3), so it can't stay lazy.
    rows = [l.rstrip("\n") for l in lines]

    # STAGE 1 -- sort by column 13 (the mapped sequence) so identical sequences
    # sit next to each other; duplicates can then be found by comparing each row
    # to the one before it. Ties break on the whole line (the `, L`).
    srt = sorted(rows, key=lambda L: (L.split("\t")[12], L))

    # STAGE 2 -- the legacy awk's roundabout way of collecting rows into `temp`.
    # Buffer each run of identical (gene, seq) in `a`; dump it at every boundary.
    # This block has quirks (drops some rows, odd counting) that DO NOT change the
    # final result -- stage 3 re-adds every original row. It exists only to
    # reproduce the legacy order/representative that downstream depends on.
    temp, a, n, l2, l13 = [], {}, 0, None, None

    for i, line in enumerate(srt):
        c = line.split("\t")
        g, s = c[1], c[12]          # g = gene (col 2), s = sequence (col 13)

        if i == 0:
            n = 0
            a = {0: line}

        # Same (gene, seq) as the previous row -> keep accumulating the run.
        if g == l2 and s == l13:
            n += 1
            a[n] = line

        else:
            # Boundary (different gene/seq): flush the buffered run into temp.
            # The current row is deliberately NOT added here -- a legacy quirk
            # that drops "singleton" rows; stage 3 recovers them from `rows`.
            for k in sorted(a):
                temp.append(a[k])
            n = 0
            a = {}

        l2, l13 = g, s

    # awk END clause: always emit the very last row.
    if srt:
        temp.append(srt[-1])

    # STAGE 3 -- where the real dedup happens. Concatenate: one prepended row
    # (the `sed -n '2p'` hack, which lets temp's 2nd row win its group) + temp +
    # ALL the original rows. Then keep the FIRST time each (gene, seq) appears.
    # Because every original row is appended, anything stage 2 dropped is still
    # caught here -> the output is exactly one row per (gene, seq).
    stream = temp[1:2] + temp + rows

    seen, out = set(), []

    for line in stream:
        c = line.split("\t")
        key = (c[1], c[12])
        if key not in seen:
            seen.add(key)
            out.append(line)

    for line in out:
        yield f"{line}\n"


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: collapse_blast <in.blast>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        for out in collapse_blast(f):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
