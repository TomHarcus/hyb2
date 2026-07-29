""" Port of hyb2_compare

"""

import argparse, sys
from hyb2 import config

def print_help():
    print("No Options Specified!")
    print("Usage:")
    print("To compare between datasets:")
    print("hyb2_compare -i <Input.table> -o <Outpud_ID>")
    print("Input.table format:")
    print("Output1.hyb   Output1.GENE.contact.txt   condition_one")
    print("Output2.hyb   Output2.GENE.contact.txt   condition_one")
    print("Output3.hyb   Output3.GENE.contact.txt   condition_two")
    print("Output4.hyb   Output4.GENE.contact.txt   condition_two")
    print("You can prepare however many rows you like (at least 2 for each condition), categorized by condition_one and condition_two")
    print(" ")
    print("  -i   Input table You MANUALLY GENERATED")
    print("  -d   Fasta file used for mapping")
    print("  -a   Gene ID of interest")
    print("  -j   Directory of VARNAcmd.jar (default set when installing hyb2)")
    print("  -0   Folding option: 0 to disable, 1 to activate automatic folding of enriched interactions (default=0)")
    print("  -m   DESeq2 chimera count filtering threshold (default=2)")
    print("  -r   Range to define an interaction (default=100). 100nt from midpoint in both directions, meaning 200nt long for each defined interaction")
    print("  -q   Upper limit for heatmap chimeric count (default=0.95)")
    print(" ") 
    print("Currently Not Compatible with Mac-OS-ARM due to several packages not being implemented yet.")
    print(" ")
    print("Run hyb2 without options for details about hyb2, hyb2_coverage, and hyb2_app")
    print(" ")
    print("Any other queries, email me at laujianyou@live.com")

    return 0

def hyb2_compare(input_table, out, min_reads, interaction_range, LIMIT, GENE, FASTA, 
                 VARNA, FOLDING):

    

    pass

def build_parser() -> argparse.ArgumentParser:
    
    p = argparse.ArgumentParser(
        prog="hyb2-compare",
        description="finds significant differences between two sets of data",
        add_help=False,
    )
    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-i", dest="input_table", required=True, metavar="INPUT.HYB", help="Input table You MANUALLY GENERATED")
    p.add_argument("-o", dest="out", required=True, help="out destination")
    p.add_argument("-m", dest="min_reads", type=int, default=2, help="DESeq2 chimera count filtering threshold (default=2)")
    p.add_argument("-r", dest="interaction_range", type=int, default=100, help="Range to define an interaction (default=100). 100nt from midpoint in both directions, meaning 200nt long for each defined interaction")
    p.add_argument("-q", dest="LIMIT", type=float, default=0.95, help="Upper limit for heatmap chimeric count (default=0.95)")
    p.add_argument("-a", dest="GENE", help="Gene ID of interest")
    p.add_argument("-d", dest="FASTA", help="Fasta file used for mapping")
    p.add_argument("-j", dest="VARNA", default=None, help="Directory of VARNAcmd.jar")
    p.add_argument("-0", dest="FOLDING", type=int, default=0, help="Folding option: 0 to disable, 1 to activate automatic folding of enriched interactions (default=0)")
    return p


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv

    if not argv:
        print_help()
        return 0

    args = build_parser().parse_args(argv)
    hyb2_compare(
        args.input_table, args.out, args.min_reads, args.interaction_range,
        args.LIMIT, args.GENE, args.FASTA, args.VARNA, args.FOLDING
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())