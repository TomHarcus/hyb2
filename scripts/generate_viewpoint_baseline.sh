#!/bin/bash -e
# Regenerates fixtures/viewpoint_run/ -- parity oracles for the two deterministic
# leaves of plot_viewpoint: hyb2blast.awk and blast2gplot.pl. (viewpoint_graph.R
# itself isn't golden-tested, same as VARNA/CDM rendering.)
#
# Usage:  scripts/generate_viewpoint_baseline.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TIER2="$REPO_ROOT/fixtures/tier2_run"
OUT="$REPO_ROOT/fixtures/viewpoint_run"
HYB="$TIER2/test.ua.hyb"
REF_FASTA="$REPO_ROOT/data/Zika_18S_formatted.fasta"
GENE=Zika_virusRNA
HB="$REPO_ROOT/legacy_bin/hyb2blast.awk"

[ -f "$HYB" ] || { echo "Error: $HYB not found (run scripts/generate_tier2_baseline.sh first)" >&2; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"

echo "[1/2] hyb2blast.awk on test.ua.hyb -> hyb2blast.golden"
awk -f "$HB" "$HYB" > "$OUT/hyb2blast.golden"

echo "[2/2] blast2gplot.pl inputs + golden (single-gene $GENE, EXP=exp)"
# reproduce plot_viewpoint's single-gene prep to get faithful blast2gplot inputs
LEN=$(expr $(grep -A1 "$GENE" "$REF_FASTA" | awk 'FNR==2' | wc -c) - 1)
printf '%s\t%s\n' "$GENE" "$LEN" > "$OUT/blast2gplot.lengths.txt"
awk -f "$HB" "$HYB" | grep "$GENE" | head -n1 > "$OUT/blast2gplot.ref.blast"
awk -v G="$GENE" '$4~G && $10~G' "$HYB" | awk -f "$HB" > "$OUT/blast2gplot.blast"

# blast2gplot.pl writes <EXP>_<gene>.gplot into cwd, so run it in a scratch dir
d="$(mktemp -d)"
cp "$OUT/blast2gplot.ref.blast"   "$d/ref.blast"
cp "$OUT/blast2gplot.blast"       "$d/in.blast"
cp "$OUT/blast2gplot.lengths.txt" "$d/lengths.txt"
( cd "$d" && perl "$REPO_ROOT/legacy_bin/blast2gplot.pl" \
    EXP=exp N_GENES=1 REF_BLAST_FILE=ref.blast BLAST_FILE=in.blast GENE_LENGTHS_FILE=lengths.txt )
cp "$d/exp_${GENE}.gplot" "$OUT/blast2gplot.gplot.golden"
rm -rf "$d"

echo "done -> $OUT"
ls -la "$OUT"
