""" Similarity contact-map merge/min (shared by hyb2_compare and plot_similarity_map)

DELIBERATE DEVIATION from the legacy similarity chain (bin/hyb2_compare:106,
identical in bin/plot_similarity_map:38). The legacy awk pipeline
    awk 'n=x[$1,$2]{print n"\n"$0} {x[$1,$2]=$0}' | awk 'ORS=NR%2?FS:RS' | ...
is buggy: it (a) emits a malformed 4-field first row from its FNR==1 branch, and
(b) pairs each key with only its immediately-previous occurrence, so a key seen k
times yields k-1 rows of pairwise-consecutive minima (with non-min values) instead
of one true minimum. Both feed garbage into similarity_heatmap.R. Instead take the
clean global minimum count per (x,y) key, one row each.
NOTE: changes results vs legacy -> diffing this against the legacy chain won't match.
"""

def similarity_contact(merged_lines):
    groups = {}
    for line in merged_lines:
        if any(c.isalpha() for c in line):
            continue

        cols = line.split()
        if len(cols) < 3:
            continue

        key = (cols[0], cols[1])

        groups.setdefault(key, []).append(cols[2])

    out_rows = [f"{a}\t{b}\t{min(counts, key=float)}" for (a,b), counts in groups.items() if len(counts) > 1]

    return out_rows
