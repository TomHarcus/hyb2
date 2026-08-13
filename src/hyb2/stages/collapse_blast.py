"""Port of collapse_blast_2.sh 

removes duplicate rows, keeping exactly one row per
unique (gene = col 2, mapped sequence = col 13) pair. 

NOTE, despite the legacy "tallying reads to the ID" comment: NO count tallying
survives. The awk recomputes the read-id count, then `cut -f3-20` throws it away,
so every output read ID is unchanged. That path is dead code; this port omits it.

WHY THE CODE BELOW IS CONVOLUTED: collapse's exact ROW ORDER changes the
downstream numbers - mtophits keeps the first-seen e-value per read id, so a
different order -> different reference sums (reversing the rows shifted Zika
244903 -> 243644). So this faithfully reproduces the legacy script's specific
order and representative row rather than deduping cleanly. If that fragility is
ever judged not worth preserving, all of stages 1-2 could collapse to a clean
"keep first per (gene, seq) in input order" - but that changes results and needs
a re-baselined golden (a decision for Grzegorz).

The sort key is pinned to code points (locale-independent), which is actually
more reproducible than the legacy shell `sort`, whose order depended on locale.

"""


import sys, subprocess, tempfile, os
from typing import Iterable, Iterator
from hyb2 import config


def collapse_blast(in_path, out_path, *, sort_mem="4G", tmpdir=None):
    sort_bin = config.gnu_sort()
    tmpdir = tmpdir or tempfile.gettempdir()

    env = {**os.environ, "LC_ALL": "C"}     

    sort = [sort_bin, "-S", sort_mem, "-T", tmpdir, "-k13", in_path]

    sort_process = subprocess.Popen(
        sort,
        stdout=subprocess.PIPE,
        text=True,
        env=env,
        )

    current_run_buffer = []
    last_gene = None
    last_seq = None

    tmpf = tempfile.NamedTemporaryFile("w", dir=tmpdir, delete=False).name

    line = None
    with open(tmpf, "a") as temp_fin:
        for line in sort_process.stdout:
            elements = line.rstrip("\n").split("\t")

            if last_gene is None and last_seq is None:
                last_gene = elements[1]
                last_seq = elements[12]
                current_run_buffer.append(line)

            elif elements[1] == last_gene and elements[12] == last_seq:
                current_run_buffer.append(line)
                
            else:
                for l in current_run_buffer:
                    temp_fin.write(f"{l}")
                    
                current_run_buffer.clear()

                last_gene = elements[1]
                last_seq = elements[12]


        # awk END clause: keep ONLY the last line ($0), not the leftover buffer
        if line is not None:
            temp_fin.write(line)

    numbered = tempfile.NamedTemporaryFile("w", dir=tmpdir, delete=False).name

    counter = 0
    with open(numbered, "w") as out:

        with open(tmpf) as f:
            f.readline()
            line2 = f.readline()
        if line2:
            counter += 1
            out.write(f"{counter}\t{line2}")

        with open(tmpf) as f:
            for line in f:
                counter += 1
                out.write(f"{counter}\t{line}")

        with open(in_path) as f:
            for line in f:
                counter += 1
                out.write(f"{counter}\t{line}")

    sorted1 = tempfile.NamedTemporaryFile("w", dir=tmpdir, delete=False).name

    with open(sorted1, "w") as fout:
        subprocess.run(
            [sort_bin, "-S", sort_mem, "-T", tmpdir, "-t", "\t",
             "-k3,3", "-k14,14", "-k1,1n", numbered],
             stdout=fout, env=env, check=True
        )

    deduped = tempfile.NamedTemporaryFile("w", dir=tmpdir, delete=False).name
    prev_key = None

    with open(sorted1) as fin, open(deduped, "w") as fout:
        for line in fin:
            elements = line.split("\t")
            key = (elements[2], elements[13])

            if key != prev_key:
                fout.write(line)
            prev_key = key

    resorted = tempfile.NamedTemporaryFile("w", dir=tmpdir, delete=False).name

    with open(resorted, "w") as fout:
        subprocess.run(
        [sort_bin, "-S", sort_mem, "-T", tmpdir, "-t", "\t", "-k1,1n", deduped],
        stdout=fout, env=env, check=True,
        )

    with open(resorted) as fin, open(out_path, "w") as fout:
        for line in fin:
            fout.write(line.split("\t", 1)[1])

    for f in (tmpf, numbered, sorted1, deduped, resorted):
        os.remove(f)







"""
def collapse_blast(lines: Iterable[str]) -> Iterator[str]:
    # Load every row up front: this needs two passes over the data (the sort
    # below, then `rows` is reused whole in stage 3), so it can't stay lazy.
    rows = [l.rstrip("\n") for l in lines]

    # STAGE 1 - sort by column 13 (the mapped sequence) so identical sequences
    # sit next to each other; duplicates can then be found by comparing each row
    # to the one before it. Ties break on the whole line (the `, L`).
    srt = sorted(rows, key=lambda L: (L.split("\t")[12], L))

    # STAGE 2 - the legacy awk's roundabout way of collecting rows into `temp`.
    # Buffer each run of identical (gene, seq) in `a`; dump it at every boundary.
    # This block has quirks (drops some rows, odd counting) that DO NOT change the
    # final result - stage 3 re-adds every original row. It exists only to
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
            # The current row is deliberately NOT added here - a legacy quirk
            # that drops "singleton" rows; stage 3 recovers them from `rows`.
            for k in sorted(a):
                temp.append(a[k])
            n = 0
            a = {}

        l2, l13 = g, s

    # awk END clause: always emit the very last row.
    if srt:
        temp.append(srt[-1])

    # STAGE 3 - where the real dedup happens. Concatenate: one prepended row
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
"""

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print("usage: collapse_blast <in.blast> <out.blast>", file=sys.stderr)
        return 1
    collapse_blast(argv[0], argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
