"""Score one RNA pair without merging the coordinate systems of distinct RNAs."""

from pathlib import Path
from urllib.parse import quote

from hyb2.common.fasta2tab import fasta_to_tab
from hyb2.common.histogram import histogram


def paired_coordinates(ct_text, records):
    """Resolve opaque CT identifiers against original HYB arm metadata.

    Only pairs crossing between the two cofolded arms are evidence for this
    interaction. Original names are never reconstructed by splitting delimiters.
    """
    for line in ct_text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if "dG" in fields or "ENERGY" in fields:
            start1, end1, start2, end2, reverse = records[fields[-1]]
            length1 = end1 - start1 + 1
            length2 = end2 - start2 + 1
            continue
        position, partner = int(fields[0]), int(fields[4])
        if not (1 <= position <= length1 < partner <= length1 + length2):
            continue
        x, y = start1 + position - 1, start2 + partner - length1 - 1
        if reverse:
            x, y = y, x
        yield f"{x}\t{y}\n"


def score_gene_pair(in_hyb, ref_fasta, gene1, gene2, *, fold="vienna",
                    vienna_bin=None, output=None):
    """Write a labelled score table in each RNA's original coordinates.

    This is a scoring output, not a single-sequence folding constraint file.
    All cross-arm pairs for the selected genes are scored, without a window cap.
    """
    from hyb2.folding.comrades_make_constraints import _fold_vienna, _fold_unafold

    if fold not in ("vienna", "cplfold", "unafold"):
        raise ValueError(f"unknown folder: {fold}")
    sequences = dict(line.split("\t", 1) for line in
                     fasta_to_tab(Path(ref_fasta).read_text()).splitlines())
    for gene in (gene1, gene2):
        if gene not in sequences:
            raise ValueError(f"Gene not found in reference: {gene}")

    if output is None:
        pair = f"{quote(gene1, safe='')}.{quote(gene2, safe='')}"
        output = str(Path(in_hyb).with_suffix("")) + f".{pair}.basepair_scores_intergenic.txt"
    bit1, bit2, ct = (str(output) + suffix for suffix in (".bit_1.fasta", ".bit_2.fasta", ".ct"))
    records = {}
    with open(in_hyb) as source, open(bit1, "w") as first, open(bit2, "w") as second:
        for number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                raise ValueError(f"Invalid HYB row {number}: expected at least 14 columns")
            pair = fields[3], fields[9]
            if pair == (gene1, gene2):
                reverse = False
            elif pair == (gene2, gene1):
                reverse = True
            else:
                continue
            start1, end1, start2, end2 = map(int, (fields[6], fields[7], fields[12], fields[13]))
            for gene, start, end in ((pair[0], start1, end1), (pair[1], start2, end2)):
                if not 1 <= start <= end <= len(sequences[gene]):
                    raise ValueError(f"Invalid coordinates for {gene} on HYB row {number}")
            token = f"hyb2pair{number}"
            records[f"{token}a-{token}b"] = start1, end1, start2, end2, reverse
            first.write(f">{token}a\n{sequences[pair[0]][start1-1:end1]}\n")
            second.write(f">{token}b\n{sequences[pair[1]][start2-1:end2]}\n")

    if records:
        if fold == "unafold":
            _fold_unafold(bit1, bit2, ct)
        else:
            _fold_vienna(bit1, bit2, ct, vienna_bin)
        scores = histogram(paired_coordinates(Path(ct).read_text(), records))
    else:
        scores = ()
    with open(output, "w") as result:
        result.write(f"{gene1}\t{gene2}\tscore\n")
        result.writelines(scores)
    return str(output)
