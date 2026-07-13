"""Port of bin/create_reference_file.pl (Tier 1 #1).

Reads a *_mtophits.blast, emits the ranking reference table (gene name +
criteria) consumed by remove_duplicate_hybrids. No dependencies.

Parity: fixtures/sam_composition_run/test_mtophits.ref
"""

import sys
from typing import Iterable, Iterator
import re


def create_reference(lines: Iterable[str]) -> Iterator[str]:
    """Transform mtophits.blast lines into reference-table lines."""
    
    total_cnt = 0
    seen = set()
    cnt = {}

    skip = re.compile(r"_pr-tr|pseudo|different|all_hits")

    for raw in lines:
        line = raw.rstrip("\n")
        if not line:
            continue
        
        columns = line.split("\t")

        readID_int = columns[0].split("_")


        if columns[0] not in seen:
            seen.add(columns[0])
            total_cnt += int(readID_int[-1])

        if int(columns[9]) > int(columns[8]) and not skip.search(line):
            cnt[columns[1]] = cnt.get(columns[1], 0) + int(readID_int[-1])

    for gene, c in sorted(cnt.items(), key=lambda kv: -kv[1]):
        yield f"{gene}\t{c}\t{100 * c / total_cnt:.15g}\n"



def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: create_reference_file <mtophits.blast>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        for out in create_reference(f):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
