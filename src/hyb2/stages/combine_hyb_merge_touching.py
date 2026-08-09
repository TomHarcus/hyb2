""" Port of combine_hyb_merge_touching.pl

NOTE on main(): the legacy Perl's FOO=bar switch-parsing line (the one every
other ported script uses to read EXP/TARGET/GUIDE/etc. from argv) is commented
out in combine_hyb_merge_touching.pl - verified directly, not assumed. So in
the actual legacy script today, TARGET/GUIDE/TWO_WAY_MERGE/PRINT_SEQ_IDS can
never be set from the command line; they're permanently stuck at their
defaults. main() below re-enables that switch parsing anyway, since the
underlying function already cleanly supports it via kwargs - this is a
deliberate, documented enhancement over current legacy behaviour, not a
translation of something the Perl actually does today.
"""

import re
import sys
from typing import Iterator

from hyb2.models import Hybrid2

def combine_hyb_merge_touching(
        hyb_lines,
        *,
        exp="default",
        target=".",
        guide=".",
        two_way_merge=False,
        print_seq_ids=False
    ) -> Iterator[str]:

    all_hits = []
    bucket = {}

    for line in hyb_lines:
        curr_hit = Hybrid2()
        curr_hit.initialize_hyb(line, exp)

        if not curr_hit.check_hit_names(target, guide):
            continue

        reversed_hit = None
        if two_way_merge:
            reversed_hit = Hybrid2()
            reversed_hit.initialize_hyb(line, exp)
            reversed_hit.reverse_bit_order()

        for old_hit in bucket.get(curr_hit.get_sorted_bit_names(), []):
            if old_hit.touches(curr_hit):
                curr_hit.found_overlap(2)
                old_hit.merge_with(curr_hit)
                break
            elif two_way_merge and old_hit.touches(reversed_hit):
                reversed_hit.found_overlap(2)
                curr_hit.found_overlap(2)
                old_hit.merge_with(reversed_hit)
                old_hit.twoway_overlap(1)
                break

        if curr_hit.found_overlap() <= 1:
            all_hits.append(curr_hit)
            bucket.setdefault(curr_hit.get_sorted_bit_names(), []).append(curr_hit)

    for hit in all_hits:
        seq_ids = ",".join(hit.seq_ID_list())
        out = hit.print_hyb_15_columns()
        out += f"\tcount_total={hit.count()};count_last_clustering={hit.found_overlap()};two_way_merged={hit.twoway_overlap()};"

        if print_seq_ids:
            out += f"seq_IDs_in_cluster={seq_ids};"

        yield out + "\n"


_SWITCHES = {
    "EXP": ("exp", str),
    "TARGET": ("target", str),
    "GUIDE": ("guide", str),
    "TWO_WAY_MERGE": ("two_way_merge", lambda v: bool(int(v))),
    "PRINT_SEQ_IDS": ("print_seq_ids", lambda v: bool(int(v))),
}


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    kwargs = {}
    while argv and re.match(r"^[A-Za-z_0-9]+=", argv[0]):   # legacy-style FOO=bar switches
        key, _, val = argv[0].partition("=")
        if key in _SWITCHES:
            name, conv = _SWITCHES[key]
            kwargs[name] = conv(val)
        argv = argv[1:]
    if len(argv) != 1:
        print(
            "usage: combine_hyb_merge_touching [EXP=..] [TARGET=..] [GUIDE=..] "
            "[TWO_WAY_MERGE=1] [PRINT_SEQ_IDS=1] <in.hyb>",
            file=sys.stderr,
        )
        return 1
    with open(argv[0]) as f:
        for out in combine_hyb_merge_touching(f, **kwargs):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())