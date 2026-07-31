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

import subprocess, math, re, glob
from pathlib import Path
from hyb2.stages.make_hybrid_annotation_table import make_hybrid_annotation_table
from hyb2.stages.DESeq_interaction_split_select import split_select

def _num(s):
    try: return float(s)
    except ValueError: return 0.0

def hyb2_compare(input_table, out, min_reads, interaction_range, LIMIT, GENE, FASTA, 
                 VARNA, FOLDING):

    rows = [tuple(l.split()) for l in open(input_table) if l.strip()] 

    for x, y, z in rows:
        
        name = y.replace(".contact.txt", "")

        lines = [l.rstrip("\n").split("\t") for l in open(y) if not l.startswith("#")]

        lines.sort(key=lambda c: _num(c[2]), reverse=True)

        content = "\n".join(f"{c[0]}_{c[1]}\t{name}={c[2]}" for c in lines) + "\n"
        Path(f"{name}.forTable.txt").write_text(content)

    fortable_files = [y.replace(".contact.txt", ".forTable.txt")
                        for (x, y, z) in rows]

    Path(f"{out}.forTable.list.txt").write_text("\n".join(fortable_files) + "\n")

    names = [y.replace(".contact.txt", "") for (x, y, z) in rows]
    header = ("#seq_ID\t" + "\t".join(names)).replace("-", "_")
    Path(f"{out}.table.txt").write_text(header + "\n")

    body = [r.replace("NA", "0") for r in make_hybrid_annotation_table(fortable_files).splitlines()[1:]]

    with open(f"{out}.table.txt", "a") as f:
        f.write("\n".join(body) + "\n")

    names_table = "\n".join(
        f"{y.replace(".contact.txt", "")}\t{z}".replace("-", "_")
        for (x, y, z) in rows
    ) + "\n"

    Path(f"{out}_names.table").write_text(names_table)

    subprocess.run(["Rscript", config.rscript("DESeq_run.R"),
                    f"{out}.table.txt", f"{out}_names.table",
                    str(min_reads)],
                    check=True)

    deseq_output = []
    for line in open(f"DESeq_{out}.txx").read().splitlines()[1:]:
        line = line.replace('"', "")
        columns = line.split()

        deseq_output.append(f"{columns[0]}\tbaseMean={columns[1]};log2FoldChange={columns[2]};lfcSE={columns[3]};stat={columns[4]};pvalue={columns[5]};padj={columns[6]}")
        
    Path(f"{out}.tmp1").write_text("\n".join(deseq_output) + "\n")

    lines = open(f"{out}.table.txt").read().splitlines()
    val = lines[0].split("\t")[1:]
    transpose_output = []
    for line in lines[1:]:
        columns = line.split("\t")
        counts = columns[1:]
        transpose_output.append(columns[0] + "\t" + "".join(f"{n}={c};" for n, c in zip(val, counts)))

    Path(f"{out}.tmp2").write_text("\n".join(transpose_output) + "\n")

    merged = make_hybrid_annotation_table([f"{out}.tmp2", f"{out}.tmp1"]).splitlines()
    header, data = merged[0], merged[1:]

    kept = []
    for row in data:
        columns = row.split("\t")
        if columns[-1] == "NA" or float(columns[-1]) >= 0.05:
            continue
        kept.append(row)

    kept.sort(key=lambda r: float(r.split("\t")[-1]))

    Path(f"DESeq_{out}_significant.txx").write_text("\n".join([header] + kept) + "\n")

  
    emitted = []
    # skip the header (splitlines()[1:]): the legacy runs it through the awk too,
    # where "log2FoldChange">"0" is a true STRING comparison so the header hits the
    # >0 branch and is then dropped by `awk 'NR>1'` -- i.e. NR>1 removes the header,
    # not a data row. Skipping it here is the equivalent, so there is NO extra [1:].
    for row in open(f"DESeq_{out}_significant.txx").read().splitlines()[1:]:
        cols = row.split()

        if float(cols[-5]) < 0:
            emitted.append(f"{cols[0]}\t{math.log(float(cols[-1]))/math.log(10):.6g}\t{cols[-5]}\t{cols[-1]}")

        elif float(cols[-5]) > 0:
            emitted.append(f"{cols[0]}\t{-1*math.log(float(cols[-1]))/math.log(10):f}\t{cols[-5]}\t{cols[-1]}")

    tail = [r.replace("_", "\t") for r in emitted]
    tail = ["x\ty\tlogpadj\tlog2FoldChange\tpadj"] + tail

    Path(f"DESeq_{out}_significant.padj_heatmap.txt").write_text("\n".join(tail) + "\n")

    subprocess.run(["Rscript", config.rscript("differential_coverage_map.R"),
                    f"DESeq_{out}_significant.padj_heatmap.txt", out],
                    check=True)

    merged = "".join(open(y).read() for (x, y, z) in rows)
    Path(f"{out}.merge.txt").write_text(merged)

  
    groups = {}
    for line in open(f"{out}.merge.txt").read().splitlines():
        if any(c.isalpha() for c in line):
            continue

        cols = line.split()
        if len(cols) < 3:
            continue

        key = (cols[0], cols[1])

        groups.setdefault(key, []).append(cols[2])

    out_rows = [f"{a}\t{b}\t{min(counts, key=float)}" for (a,b), counts in groups.items() if len(counts) > 1]

    Path(f"{out}.contact.txt").write_text("\n".join(out_rows) + "\n")

    subprocess.run(["Rscript", config.rscript("similarity_heatmap.R"),
                    f"{out}.contact.txt", str(LIMIT)], check=True)

    split_select(out, interaction_range)


    for sign in ("pos", "neg"):
        f = f"{out}_{interaction_range}range_{sign}_enrichment.txt"
        Path(f).write_text(
            Path(f).read_text().replace("0\t0\t0.0\t0.0", "x\ty\tlog2FoldChange\tpadj")
        )

    for sign in ("pos", "neg"):
        src = f"{out}_{interaction_range}range_{sign}_enrichment.txt"
        windows = []
        for row in Path(src).read_text().splitlines()[1:]:       
            c = row.split("\t")
            x, y = int(c[0]), int(c[1])
            windows.append(f"{x-150}\t{x+150}\t{y-150}\t{y+150}\t{c[2]}\t{c[3]}")
        body = ["x1\tx2\ty1\ty2\tlog2FoldChange\tpadj"] + windows   
        body = [re.sub(r"-.[0-9]+\t", "1\t", ln, count=1) for ln in body]  
        body = [ln.replace("\t\t", "\t-") for ln in body]        
        Path(f"{out}_{interaction_range}range_{sign}_enrichment.heatmap.txt").write_text(
            "\n".join(body) + "\n"
        )

    value = [contact.replace(".contact.txt", "") for (hyb, contact, cond) in rows]

    for sign in ("pos", "neg"):
        hm = f"{out}_{interaction_range}range_{sign}_enrichment.heatmap.txt"
        for row in Path(hm).read_text().splitlines()[1:11]:       # awk NR>1 && NR<12
            c = row.split("\t")
            subprocess.run(
                ["Rscript", config.rscript("contact_density_map_zoom.R"),
                 c[0], c[1], c[2], c[3], str(LIMIT), *value],
                check=True,
            )

    condition_one_files = [hyb for (hyb, contact, cond) in rows if cond == "condition_one"]
    condition_two_files = [hyb for (hyb, contact, cond) in rows if cond == "condition_two"]

    if FOLDING == 1:
        _fold_enriched(condition_one_files, "pos", out, interaction_range, GENE, FASTA, VARNA)
        _fold_enriched(condition_two_files, "neg", out, interaction_range, GENE, FASTA, VARNA)
    else:
        print("Use Options -0 to -9, and -j to Plot RNA Secondary Structures for Enriched Interactions")

    print("Comparison Completed")


def _fold_enriched(condition_files, sign, out, rng, GENE, FASTA, VARNA):
    """Port of bin/hyb2_compare lines 137-141 (the FOLDING==1 branch, one condition).
    Builds a <IN>.<sign>.hyb of GENE-GENE chimeras for each condition file, then folds
    the top-10 enriched interactions of every heatmap. NOTE (faithful legacy quirk): the
    fold loop uses `in_hyb` = the LAST condition file, exactly as the legacy $IN leaks out
    of the build loop above it."""
    in_hyb = None
    for src in condition_files:                                 
        filtered = []
        for path in glob.glob(src + "*"):                      
            for line in open(path):
                c = line.split("\t")
                if len(c) >= 10 and re.search(GENE, c[3]) and re.search(GENE, c[9]):
                    filtered.append(line.rstrip("\n"))
        in_hyb = src.replace("hyb", f"{sign}.hyb", 1)          
        Path(in_hyb).write_text("\n".join(filtered) + "\n")

    for hm in glob.glob(f"{out}_{rng}*enrichment.heatmap.txt"):  
        lines = Path(hm).read_text().splitlines()
        for i, row in enumerate(lines):
            if not (1 < i + 1 < 12):                             
                continue
            c = row.split("\t")
            if int(c[0]) > int(c[2]):                           
                x1, y1 = int(c[2]), int(c[0])
            elif int(c[0]) < int(c[2]):
                x1, y1 = int(c[0]), int(c[2])
            else:
                continue
            res = y1 - x1
            cmd = [sys.executable, "-m", "hyb2.pipelines.hyb2_fold",
                   "-i", in_hyb, "-a", GENE, "-d", FASTA, "-x", str(x1 + 1)]
            if 300 > res:                                        
                cmd += ["-l", str(res + 300)]
            else:                                                
                cmd += ["-y", str(y1 + 1), "-l", "300"]
            if VARNA:
                cmd += ["-j", VARNA]
            subprocess.run(cmd, check=True)


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