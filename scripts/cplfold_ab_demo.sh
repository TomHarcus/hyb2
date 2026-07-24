#!/bin/bash -e
# CPLfold bonus-normalization A/B demo (for the CoupleFold walkthrough).
#
# Folds ONE fragment three ways and prints each structure + energy so the
# raw-vs-log normalization decision can be made side by side:
#
#   1. no bonus   -- plain pseudoknot fold (the control / baseline)
#   2. raw bonus  -- --normalize raw --alpha 0.3  (raw counts dominate unless alpha is tiny)
#   3. log bonus  -- --normalize log --alpha 0.5  (log1p behaves well at alpha 0.5-1.0)
#
# The support matrix + fragment are built ONCE (comrades_make_constraints), then
# reused for all three folds so the ONLY thing that changes is the bonus handling.
#
# Defaults to the dense Zika_virusRNA [3900,4199] window (where the bonus pulls in
# an experimentally-supported pseudoknot at better energy than no-bonus).
#
# Usage:
#   scripts/cplfold_ab_demo.sh [HYB] [GENE] [X_START] [LENGTH]
# e.g.
#   scripts/cplfold_ab_demo.sh                       # Zika 3900 300 on the test data
#   scripts/cplfold_ab_demo.sh ~/data/mouse.hyb 18S_rRNA 1200 250
#
# Requires the hyb2 env active (RNAcofold on PATH) + a runnable HotKnots binary
# (CPLfold's computeEnergy, built for this machine -- see CLAUDE.md).

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VIENNA_BIN="${VIENNA_BIN:-$HOME/miniconda3/envs/hyb2/bin}"
PY="$VIENNA_BIN/python"

HYB="${1:-$REPO_ROOT/fixtures/tier2_run/test.ua.hyb}"
GENE="${2:-Zika_virusRNA}"
X="${3:-3900}"
LEN="${4:-300}"
REF="$REPO_ROOT/data/Zika_18S_formatted.fasta"
END=$((X + LEN - 1))

for f in "$HYB" "$REF"; do
    [ -f "$f" ] || { echo "Error: not found: $f" >&2; exit 1; }
done
[ -x "$VIENNA_BIN/RNAcofold" ] || { echo "Error: RNAcofold not in $VIENNA_BIN (set VIENNA_BIN=)" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp "$HYB" "$WORK/in.hyb"
cp "$REF" "$WORK/ref.fasta"
cd "$WORK"

echo "=== CPLfold bonus A/B: $GENE [$X,$END] (len $LEN) on $(basename "$HYB") ==="
echo "[setup] building support matrix + fragment (comrades_make_constraints)..."
PYTHONPATH="$REPO_ROOT/src" PATH="$VIENNA_BIN:$PATH" "$PY" -m hyb2.pipelines.comrades_make_constraints \
    -i in.hyb -f ref.fasta -b "$X" -e "$END" >/dev/null 2>&1

FRAG="ref_${X}-${END}.fasta"          # fragment fasta (comrades_make_constraints output)
SCORES="in.basepair_scores.txt"       # support matrix (genome coords)
[ -f "$FRAG" ]   || { echo "Error: fragment $FRAG not produced" >&2; exit 1; }
[ -f "$SCORES" ] || { echo "Error: scores $SCORES not produced" >&2; exit 1; }
echo "[setup] fragment len $(awk 'NR==2{print length($0)}' "$FRAG"), $(grep -c . "$SCORES") score rows"
echo

# fold(label, extra-flags...) -> runs comrades_fold cplfold, prints structure + energy.
# On failure it prints the error and continues to the next fold (does not abort).
fold () {
    local label="$1"; shift
    cp "$FRAG" "frag.fasta"
    rm -f frag.fasta.vienna
    printf '%-22s folding... ' "$label"; local t0=$(date +%s)
    if ! PYTHONPATH="$REPO_ROOT/src" PATH="$VIENNA_BIN:$PATH" "$PY" -m hyb2.pipelines.comrades_fold \
            -c none -i frag.fasta -r cplfold -b "$X" -e "$END" "$@" >fold.log 2>&1; then
        printf 'FAILED -- last lines of the error:\n'
        tail -6 fold.log | sed 's/^/   | /'
        echo
        return 0
    fi
    local struct energy pk
    struct=$(sed -n '3p' frag.fasta.vienna | awk '{print $1}')
    energy=$(sed -n '3p' frag.fasta.vienna | grep -oE '\([-0-9.]+\)$' | tr -d '()')
    pk=$(printf '%s' "$struct" | grep -o '\[' | wc -l | tr -d ' ')
    printf 'done in %ss  energy=%-10s pseudoknot-pairs=%s\n' "$(( $(date +%s) - t0 ))" "${energy:-N/A}" "$pk"
    printf '   %s\n' "$struct"
    if [ -n "$OUT" ]; then cp frag.fasta.vienna "$OUT/${label}.vienna"; fi   # OUT=dir to keep them
    return 0
}

echo "=== results (only the bonus handling differs) ==="
fold "1_no-bonus"
fold "2_raw_a0.3"  -p "$SCORES" --normalize raw --alpha 0.3
fold "3_log_a0.5"  -p "$SCORES" --normalize log --alpha 0.5

