#!/bin/bash -e
# Regenerates fixtures/coverage_run/ -- the parity oracle for stages/plot_hybrids_3.py
# (the deterministic contact-binning that hyb2_coverage does before shelling out to R).
#
# Captures the legacy bin/plot_hybrids_3.awk output on the real test.ua.hyb for the
# two cases hyb2_coverage uses:
#   single-gene  (GENE_1 == GENE_2)          -> plot_hybrids_3.single.golden
#   two-gene     (with the line-36 swap prep)-> plot_hybrids_3.twogene.golden
#
# The R plotting itself is not golden-tested (same as VARNA rendering isn't).
#
# Usage:  scripts/generate_coverage_baseline.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TIER2="$REPO_ROOT/fixtures/tier2_run"
OUT="$REPO_ROOT/fixtures/coverage_run"
HYB="$TIER2/test.ua.hyb"
AWK="$REPO_ROOT/legacy_bin/plot_hybrids_3.awk"

[ -f "$HYB" ] || { echo "Error: $HYB not found (run scripts/generate_tier2_baseline.sh first)" >&2; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"

echo "[1/2] single-gene: plot_hybrids_3.awk (GENE_1=GENE_2=Zika_virusRNA)"
awk -f "$AWK" USE_ENTIRE_HYBRIDS=1 BIN_SIZE=10 \
    GENE_1=Zika_virusRNA GENE_2=Zika_virusRNA "$HYB" \
    > "$OUT/plot_hybrids_3.single.golden"

echo "[2/2] two-gene: swap prep (legacy_bin/hyb2_coverage line 36) | plot_hybrids_3.awk"
awk -v GENE_1=Zika_virusRNA \
    '{if ($4==GENE_1) print $0; if ($10==GENE_1) print $1"\t"$2"\t"$3"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' \
    "$HYB" \
  | awk -f "$AWK" USE_ENTIRE_HYBRIDS=1 BIN_SIZE=10 \
        GENE_1=Zika_virusRNA GENE_2=18S_rRNA \
    > "$OUT/plot_hybrids_3.twogene.golden"

echo "done -> $OUT"
ls -la "$OUT"
