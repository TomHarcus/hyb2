""" Port of hyb2fasta_bits_allRNAs.awk

Inputs a .tab file with reference sequences and a .hyb file (order matters -
tab first, hyb second, matching the legacy awk's NR==FNR file-1/file-2 idiom).
Outputs two paired FASTA files: bit_1.fasta (every chimera's first arm) and
bit_2.fasta (every chimera's second arm), record N in one corresponding to
record N in the other.

Output filenames are derived from the .hyb path's basename (before its
extension), same as the legacy script - that naming lives in main(), not in
the generator itself, since the generator only needs iterables of lines.
"""

import sys


def hyb2fasta_bits_allRNAs(tab_file, hyb_file):
    seq = {}
    for line in tab_file:
        columns = line.rstrip("\n").split("\t")

        seq[columns[0]] = columns[1]

    for line in hyb_file:
        columns = line.split("\t")

        miRNA_name, miRNA_start, miRNA_end = columns[3], int(columns[6]), int(columns[7])
        mRNA_name, mRNA_start, mRNA_end = columns[9], int(columns[12]), int(columns[13])
        
        if (miRNA_end < miRNA_start or mRNA_end < mRNA_start):
            continue

        miRNA_seq = seq[miRNA_name][miRNA_start-1:miRNA_end]
        mRNA_seq = seq[mRNA_name][mRNA_start-1:mRNA_end]

        miRNA_name_extended = columns[0] + "_" + miRNA_name + "_" + str(miRNA_start) + "_" + str(miRNA_end)
        mRNA_name_extended = columns[0] + "_" + mRNA_name + "_" + str(mRNA_start) + "_" + str(mRNA_end)

        yield (f">{miRNA_name_extended}\n{miRNA_seq}\n", f">{mRNA_name_extended}\n{mRNA_seq}\n")


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print("usage: hyb2fasta_bits_allRNAs <in.tab> <in.hyb>", file=sys.stderr)
        return 1
    tab_path, hyb_path = argv
    out_name = hyb_path[:-4] if hyb_path.endswith(".hyb") else hyb_path
    out1_path = out_name + ".bit_1.fasta"
    out2_path = out_name + ".bit_2.fasta"
    with open(tab_path) as tab_file, open(hyb_path) as hyb_file, \
         open(out1_path, "w") as out1, open(out2_path, "w") as out2:
        for bit1_rec, bit2_rec in hyb2fasta_bits_allRNAs(tab_file, hyb_file):
            out1.write(bit1_rec)
            out2.write(bit2_rec)
    return 0


if __name__ == "__main__":
    sys.exit(main())

