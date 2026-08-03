#!/bin/bash -e
# Regenerates fixtures/unafold_run/ -- the UNAFold (OligoArrayAux) parity oracle
# for the `-r unafold` folding backend. Captures two goldens:
#
#   1. ct2b_gk3  : a REAL hybrid-ss-min .ct (a folded Zika fragment) + the legacy
#      Ct2B_GK_3.pl HYBRID_SS_MIN=1 output on it. Pins stages/ct2b_gk3.py.
#   2. make_constraints -r 0 : basepair_scores from the legacy
#      comradesMakeConstraints_2 -r 0 (hybrid-min) on the same 20-chimera
#      comrades_mini window. Pins the port's fold="unafold" constraints path.
#
# Requires oligoarrayaux (hybrid-ss-min + hybrid-min) and the UNAFOLDDAT energy
# tables. Both default to the hyb2 conda env; override with UNAFOLD_BIN=... /
# UNAFOLDDAT=...  (fixtures/ is gitignored, regenerated per-machine.)
#
# Usage:  scripts/generate_unafold_baseline.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TIER2="$REPO_ROOT/fixtures/tier2_run"
OUT="$REPO_ROOT/fixtures/unafold_run"

UNAFOLD_BIN="${UNAFOLD_BIN:-$HOME/miniconda3/envs/hyb2/bin}"
export UNAFOLDDAT="${UNAFOLDDAT:-$HOME/miniconda3/envs/hyb2/share/oligoarrayaux}"
REF="$REPO_ROOT/data/Zika_18S_formatted.fasta"

for tool in hybrid-ss-min hybrid-min; do
    if [ ! -x "$UNAFOLD_BIN/$tool" ]; then
        echo "Error: $UNAFOLD_BIN/$tool not found." >&2
        echo "       conda install -n hyb2 -c bioconda oligoarrayaux" >&2
        echo "       or set UNAFOLD_BIN=/path/to/env/bin" >&2
        exit 1
    fi
done
[ -d "$UNAFOLDDAT" ] || { echo "Error: UNAFOLDDAT dir $UNAFOLDDAT not found." >&2; exit 1; }
[ -f "$REF" ] || { echo "Error: reference $REF not found." >&2; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"

echo "[1/2] ct2b_gk3: real hybrid-ss-min .ct + Ct2B_GK_3.pl HYBRID_SS_MIN=1 golden"
# fold a real 200-nt Zika fragment with hybrid-ss-min -> genuine .ct
tmp="$(mktemp -d)"
PATH="$REPO_ROOT/bin:$UNAFOLD_BIN:$PATH" fasta2tab.awk "$REF" > "$tmp/ref.tab"
awk 'BEGIN{FS="\t"} $1 ~ /Zika/ {print ">"$1"\n"substr($2,3900,200)}' "$tmp/ref.tab" > "$tmp/frag.fasta"
( cd "$tmp" && PATH="$UNAFOLD_BIN:$PATH" hybrid-ss-min frag.fasta >/dev/null 2>&1 )
cp "$tmp/frag.fasta.ct" "$OUT/ct2b.ct"                                   # input fixture
perl "$REPO_ROOT/bin/Ct2B_GK_3.pl" HYBRID_SS_MIN=1 "$tmp/frag.fasta.ct" \
    > "$OUT/ct2b_gk3.vienna.golden"                                      # golden
rm -rf "$tmp"

echo "[2/2] comradesMakeConstraints_2 -r 0 (hybrid-min) -> basepair_scores golden"
mini_dir="$(mktemp -d)"
head -20 "$TIER2/test.ua.hyb" > "$mini_dir/mini.hyb"
cp "$REF" "$mini_dir/ref.fasta"
( cd "$mini_dir" && PATH="$REPO_ROOT/bin:$UNAFOLD_BIN:$PATH" LC_ALL=C UNAFOLDDAT="$UNAFOLDDAT" \
    bash "$REPO_ROOT/bin/comradesMakeConstraints_2" -i mini.hyb -f ref.fasta -b 1 -e 10298 -r 0 >/dev/null 2>&1 )
cp "$mini_dir/mini.hyb"                  "$OUT/mini.hyb"                  # port re-runs on this
cp "$mini_dir/mini.basepair_scores.txt"  "$OUT/mini.basepair_scores.golden"
rm -rf "$mini_dir"

echo "done -> $OUT"
ls -la "$OUT"
echo "--- ct2b_gk3.vienna.golden ---"; cat "$OUT/ct2b_gk3.vienna.golden"
echo "--- basepair_scores (head) ---"; head -5 "$OUT/mini.basepair_scores.golden"; wc -l "$OUT/mini.basepair_scores.golden"