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
from hyb2.tools import config


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

def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print("usage: collapse_blast <in.blast> <out.blast>", file=sys.stderr)
        return 1
    collapse_blast(argv[0], argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
