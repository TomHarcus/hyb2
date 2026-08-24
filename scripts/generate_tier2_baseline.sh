#!/bin/bash -e
# Regenerates fixtures/tier2_run/ -- golden outputs for the Tier 2 (folding-
# adjacent) stages that were ported and validated interactively but had no
# committed parity fixture yet:
#
#   hyb2constraints.pl            in.hyb            -> test.constraints
#   hyb2fasta_bits_allRNAs.awk    ref.tab in.hyb   -> test.bit_1.fasta / test.bit_2.fasta
#   combine_hyb_merge_touching.pl in.hyb           -> test.combine.hyb
#
# Inputs are the real .hyb fixtures already captured by the sam_composition
# baseline (fixtures/sam_composition_run/), plus a .tab reference built from
# the underscore-formatted Zika/18S fasta (data/Zika_18S_formatted.fasta -- the
# pipe-delimited one does NOT work here, same known quirk as hyb2_fold).
#
# Deliberate local adaptations (mirroring generate_sam_composition_baseline.sh):
#   - Legacy perl scripts carry a hardcoded conda shebang, so they are invoked
#     through an explicit interpreter (PERL/AWK/PYTHON overridable) to run on
#     any machine without the hyb2 conda env active.
#
# hyb2constraints input note: test.hyb (from the sam_composition baseline) is
# used deliberately because it carries a '#'-comment header -- that header was
# the edge case that crashed an early version of the port, so it belongs in the
# golden.
#
# Usage:  scripts/generate_tier2_baseline.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$REPO_ROOT/legacy_bin"
SRC="$REPO_ROOT/fixtures/sam_composition_run"
OUT="$REPO_ROOT/fixtures/tier2_run"

PERL="${PERL:-perl}"
AWK="${AWK:-awk}"
PYTHON="${PYTHON:-python3}"

for f in test.hyb test.ua.hyb; do
    if [ ! -f "$SRC/$f" ]; then
        echo "Error: $SRC/$f not found; run scripts/generate_sam_composition_baseline.sh first" >&2
        exit 1
    fi
done

rm -rf "$OUT"; mkdir -p "$OUT"
cp "$SRC/test.hyb" "$OUT/test.hyb"
cp "$SRC/test.ua.hyb" "$OUT/test.ua.hyb"
cd "$OUT"

echo "[1/4] ref.tab: Zika_18S_formatted.fasta -> tab (name<TAB>seq)"
"$PYTHON" -c "
import sys
sys.path.insert(0, '$REPO_ROOT/src')
from hyb2.stages.fasta2tab import fasta_to_tab
sys.stdout.write(fasta_to_tab(open('$REPO_ROOT/data/Zika_18S_formatted.fasta').read()))
" > ref.tab

echo "[2/4] hyb2constraints.pl: test.hyb -> test.constraints"
"$PERL" "$BIN/hyb2constraints.pl" test.hyb > test.constraints

echo "[3/4] hyb2fasta_bits_allRNAs.awk: ref.tab + test.ua.hyb -> test.ua.bit_{1,2}.fasta"
# awk derives output basenames from the .hyb FILENAME -> test.ua.bit_1.fasta etc.
"$AWK" -f "$BIN/hyb2fasta_bits_allRNAs.awk" ref.tab test.ua.hyb

echo "[4/4] combine_hyb_merge_touching.pl: test.ua.hyb -> test.combine.hyb"
"$PERL" "$BIN/combine_hyb_merge_touching.pl" test.ua.hyb > test.combine.hyb

echo "done -> $OUT"
ls -la "$OUT"
