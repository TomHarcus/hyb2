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

from __future__ import annotations

import argparse
import sys

# Stage imports (filled in as each stage is ported):
# from hyb2.stages import (
#     sam2blast, collapse_blast, mtophits_blast, create_reference_file,
#     get_mtop_hybrids, remove_duplicate_hybrids, histogram,
# )


def run(
    in_sam: str,
    *,
    nlines: int | None = None,
    ua_hyb: str | None = None,
    hmax: int = 10,
    blast_threshold: float = 0.1,
    mode: int = 2,
    max_overlap: int = 4,
) -> None:
    """Drive the pipeline. See module docstring for the three modes."""
    # TODO(tier1): once stages land, wire:
    #   1. sam2blast(in_sam)                       -> <out>.blast
    #   2. collapse_blast | mtophits_blast          -> <out>_mtophits.blast
    #   3. create_reference_file                    -> <out>_mtophits.ref
    #   4. get_mtop_hybrids(<out>.blast, ...)       -> <out>.hyb
    #   5. remove_duplicate_hybrids(ref, hyb)       -> <out>.ua.hyb
    #   6. histogram(cut -f4,10 ua.hyb)             -> <out>.ua.hyb_stats_by_gene.txt
    #   7. histogram(single-read filter of sam)     -> <out>_tophit_by_gene.txt
    #   8. hyb2_composition_pies (Python; get from Grzegorz)
    raise NotImplementedError("sam_composition pipeline not yet wired -- Tier 1 in progress")


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
