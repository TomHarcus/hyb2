"""Port of bin/hyb2_sam_composition.sh -- the pipeline "spine".

Chimera-calling from an already-mapped SAM, then single-read vs chimera
composition stats. Parity oracle: fixtures/sam_composition_run/ (regenerate
with scripts/generate_sam_composition_baseline.sh).

Three run modes, mirroring the shell script:
  full     -i <sam>              steps 1-9 on the whole SAM
  subset   -i <sam> -n <nlines>  diagnostic run on the first <nlines> lines
  stats    -i <sam> -u <ua.hyb>  steps 7-9 only, from an existing dedup hyb

Fixed params (as in the shell): BLAST_THRESHOLD=0.1, MODE=2, MAX_OVERLAP=4.
MODE is exposed so the _antisense variant (MODE=1) folds in here as a flag
rather than a separate script.

Wiring note carried over from the shell: get_mtop_hybrids reads the RAW .blast,
while the collapse -> mtophits -> reference branch only builds the ranking
.ref consumed by remove_duplicate. Keep those two data paths separate.
"""

import argparse
import sys

from hyb2.stages.sam2blast import sam2blast
from hyb2.stages.collapse_blast import collapse_blast
from hyb2.stages.mtophits_blast import deduplicate_by_second_fragement_start
from hyb2.stages.create_reference_file import create_reference
from hyb2.stages.get_mtop_hybrids import get_mtop_hybrids
from hyb2.stages.remove_duplicate_hybrids import remove_duplicate_hybrids
from hyb2.stages.histogram import histogram


def run(
    in_sam: str,
    *,
    out=None,
    nlines: int | None = None,
    ua_hyb: str | None = None,
    hmax: int = 10,
    blast_threshold: float = 0.1,
    mode: int = 2,
    max_overlap: int = 4,
) -> None:
    """Drive the pipeline. See module docstring for the three modes."""

    print(f"input={in_sam} (whole file, no copy) h={hmax}", file=sys.stderr)

    if out is None:
        out = in_sam[:-4] if in_sam.endswith(".sam") else in_sam

    blast = out + ".blast"
    with open(in_sam) as fin, open(blast, "w") as fout:
        fout.writelines(sam2blast(fin))

    collapse = out + ".collapse.blast"
    with open(blast) as fin, open(collapse, "w") as fout:
        fout.writelines(collapse_blast(fin))

    mtophits = out + "_mtophits.blast"
    with open(collapse) as f:
        text = deduplicate_by_second_fragement_start(f.read())

    with open(mtophits, "w") as fout:
        fout.write(text)

    ref = out + "_mtophits.ref"
    with open(mtophits) as fin, open(ref, "w") as fout:
        fout.writelines(create_reference(fin))

    hyb = out + ".hyb"
    with open(blast) as fin, open(hyb, "w") as fout:
        fout.writelines(get_mtop_hybrids (
            fin, blast_threshold=blast_threshold, mode=mode, max_overlap=max_overlap, max_hits=hmax
        ))
    
    ua = out + ".ua.hyb"
    with open(ref) as r, open(hyb) as h, open(ua, "w") as fout:
        fout.writelines(remove_duplicate_hybrids(r, h, prefer_mim=True))

    gene_stats = out + ".ua.hyb_stats_by_gene.txt"
    with open(ua) as fin, open(gene_stats, "w") as fout:
        cut = (f"{c[3]}\t{c[9]}\n" for c in (line.rstrip("\n").split("\t") for line in fin))
        fout.writelines(histogram(cut))

    tophit = out + "_tophit_by_gene.txt"
    def _single_reads(sam):
        for line in sam:
            if line.startswith("@"):
                continue
            c = line.split("\t")
            flag = int(c[1])
            if (flag // 4) % 2 == 0 and (flag // 256) % 2 == 0:
                yield c[2] + "\n"
    
    with open(in_sam) as fin, open(tophit, "w") as fout:
        fout.writelines(histogram(_single_reads(fin)))

    print(f"Done. Outputs prefixed with: {out}", file=sys.stderr)

    """hyb2_composition_pies.py not here yet"""


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hyb2-sam-composition",
        description="Chimera calling + composition stats from a mapped SAM "
        "(port of bin/hyb2_sam_composition.sh).",
        add_help=False,
    )
    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-i", dest="in_sam", required=True, metavar="INPUT.SAM", help="Input SAM (required)")
    p.add_argument("-n", dest="nlines", type=int, default=None, help="Diagnostic: use only the first N SAM lines")
    p.add_argument("-h", dest="hmax", type=int, default=10, help="MAX_HITS_PER_SEQUENCE (default 10)")
    p.add_argument("-u", dest="ua_hyb", default=None, metavar="DEDUP.UA.HYB", help="Existing dedup hyb: run steps 7-9 only")
    p.add_argument("--mode", type=int, default=2, choices=[0, 1, 2, 3], help="get_mtop_hybrids MODE (2=sense only; 1=antisense allowed)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run(
        args.in_sam,
        nlines=args.nlines,
        ua_hyb=args.ua_hyb,
        hmax=args.hmax,
        mode=args.mode,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
