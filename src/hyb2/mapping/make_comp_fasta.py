""" Port of make_comp_fasta.pl

Collapses identical reads and counts occurrences, emitting one FASTA record per
unique sequence

Ordering fragility: the legacy iterates a perl hash and sorts with GNU sort, so
sequences with an equal count get an random tie order. This port uses the sorted function,
(ties keep first seen order), so on real data the output read IDs might differ
from the legacy. This is cosmetic though and shouldn't break it.

"""

import re

def make_comp_fasta(lines):

    cnt = {}
    cnt_barcodes = {}
    fnd = set()
    barcodes = 0

    for i, line in enumerate(lines):
        cols = line.rstrip("\n").split("\t")
        fld1 = cols[0]
        fld2 = cols[1] if len(cols) > 1 else ""

        if i == 0:
            barcodes = 1 if re.match(r"[^\t]*_[A-Z]*\t", line) else 0

        if barcodes == 1:
            barcode = fld1.split("_")[-1]
            cnt[fld2] = cnt.get(fld2, 0) + 1

            if (fld2, barcode) not in fnd:
                fnd.add((fld2, barcode))
                cnt_barcodes[fld2] = cnt_barcodes.get(fld2, 0) + 1

        else:
            cnt[fld2] = cnt.get(fld2, 0) + 1

    out = []

    if barcodes == 1:
        order = sorted(cnt, key=lambda s: (cnt[s], cnt_barcodes[s]), reverse=True)

        for nr, seq in enumerate(order, 1):
            out.append(f">{nr}-{cnt_barcodes[seq]}_{cnt[seq]}")
            out.append(seq)

    else:
        order = sorted(cnt, key=lambda s: cnt[s], reverse=True)

        for nr, seq in enumerate(order, 1):
            out.append(f">{nr}_{cnt[seq]}")
            out.append(seq)

    return "\n".join(out) + "\n"