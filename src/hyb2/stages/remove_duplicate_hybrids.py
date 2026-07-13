"""Port of bin/remove_duplicate_hybrids_hOH5_2.pl 

Takes the ranking reference (from create_reference_file) plus a .hyb, keeps
exactly one best hybrid per read ID. Ranking criteria (in order): sum of the
two bits' e-values; then miRNA-mRNA hybrids (when PREFER_MIM); then rank of the
higher- then lower-ranked bit in the reference; else keep the first seen.

Uses the shared Hybrid data class (models.Hybrid), the one cross-file
dependency in the Tier 1 set.

Parity: fixtures/sam_composition_run/test.ua.hyb
"""

import sys
import re

from hyb2.models import Hybrid  


def get_ordered_ranks(name1, name2, rank):
    r1 = rank.get(name1)
    r2 = rank.get(name2)

    if r1 is None and r2 is None:
        return (None, None)
    if r2 is None:
        return (r1, None)
    if r1 is None:
        return (r2, None)
    if r1 <= r2:
        return (r1, r2)
    return (r2, r1)

def rank_two_hybrids(x, y, rank, prefer_mim=True):
    x_mRNA, x_miRNA = x.match_bit_name("_mRNA"), x.match_bit_name("_microRNA")
    y_mRNA, y_miRNA = y.match_bit_name("_mRNA"), y.match_bit_name("_microRNA")

    x_top, x_bottom = get_ordered_ranks(*x.get_bit_names(), rank)
    y_top, y_bottom = get_ordered_ranks(*y.get_bit_names(), rank)

    # sum of e-values: lower wins
    xs, ys = x.sum_e_values(), y.sum_e_values()
    if xs != ys:
        return x if xs < ys else y
    
    # prefer miRNA-mRNA hybrids
    if prefer_mim:
        x_mim = x_mRNA is not None and x_miRNA is not None
        y_mim = y_mRNA is not None and y_miRNA is not None
        if x_mim and not y_mim:
            return x
        if y_mim and not x_mim:
            return y

    # rank of the higer-ranked bit
    if x_top is not None and y_top is not None and x_top != y_top:
        return x if x_top < y_top else y
    if x_top is not None and y_top is None:
        return x
    if x_top is None and y_top is not None:
        return y
    
    # rank of the lower-ranked bit
    if x_bottom is not None and y_bottom is not None and x_bottom != y_bottom:
        return x if x_bottom < y_bottom else y
    if x_bottom is not None and y_bottom is None:
        return x
    if x_bottom is None and y_bottom is not None:
        return y
    
    # otherwise keep the first
    return x

def remove_duplicate_hybrids(ref_lines ,hyb_lines, prefer_mim = True):
    """Yield the single best .hyb line per read ID."""
    rank = {}
    for i, line in enumerate(ref_lines, start=1):
        rank[line.split("\t")[0]] = i

    group, previous_id = [], None
    for line in hyb_lines:
        if line.startswith("#") or not re.search(r"[a-zA-Z0-9]", line):
            continue
        h = Hybrid()
        h.initialize_hyb(line, "default", ".", ".")
        if h.seq_ID() != previous_id:
            if group:
                yield _emit_best(group, rank, prefer_mim)
            group, previous_id = [], h.seq_ID()
        group.append(h)
    if group:
        yield _emit_best(group, rank, prefer_mim)

def _emit_best(group, rank, prefer_mim):
    best = group[0]
    for h in group[1:]:
        best = rank_two_hybrids(best, h, rank, prefer_mim)
    return best.line() + "\n"


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    prefer_mim = True
    while argv and re.match(r"^[A-Za-z_0-9]+=", argv[0]):   # legacy FOO=bar switches
        key, _, val = argv[0].partition("=")
        if key == "PREFER_MIM":
            prefer_mim = bool(int(val))
        argv = argv[1:]
    if len(argv) != 2:
        print("usage: remove_duplicate_hybrids [PREFER_MIM=1] <ref> <hyb>", file=sys.stderr)
        return 1
    with open(argv[0]) as ref, open(argv[1]) as hyb:
        for out in remove_duplicate_hybrids(ref, hyb, prefer_mim):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
