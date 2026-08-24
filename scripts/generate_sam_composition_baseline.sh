#!/bin/bash -e
# Regenerates fixtures/sam_composition_run/ -- the golden-output regression
# harness for the *new* pipeline spine that the team actually runs:
# hyb2_sam_composition.sh (sent by Grzegorz; not part of the original bin/).
#
# This is the parity oracle the Python port must reproduce stage-by-stage.
# It replaces the old hyb2-orchestrator oracle (scripts/generate_fixtures.sh,
# fixtures/legacy_run/) as the primary target: hyb2_sam_composition.sh is the
# manually-verified path used for real analyses, whereas bin/hyb2 has known
# dead branches (see fixtures/README.md).
#
# The steps below are copied verbatim from hyb2_sam_composition.sh (full-run
# mode, steps 2-9), with three deliberate local adaptations, all noted inline:
#   1. TMPDIR is set to a local scratch dir (the script hardcodes a cluster path).
#   2. Tools are invoked with explicit interpreters (bin/ perl scripts carry a
#      hardcoded conda shebang; sam2blast_3 wants `python`), so this runs on any
#      machine regardless of whether the hyb2 conda env is active.
#   3. Step 9 (hyb2_composition_pies.py) is SKIPPED -- that script is not in the
#      repo or hyb2_scripts/; Grzegorz has it on the cluster. Ask him for it to
#      capture the final pie-chart output. Steps 2-8 produce every intermediate.
#
# Fixed params, exactly as in hyb2_sam_composition.sh:
#   BLAST_THRESHOLD=0.1  MODE=2  MAX_OVERLAP=4  MAX_HITS_PER_SEQUENCE=10
#
# Usage:  scripts/generate_sam_composition_baseline.sh [input.sam]
#   input.sam defaults to ../hyb2_baseline/test.sam (the SAM Tom was given).

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$REPO_ROOT/legacy_bin"
# Grzegorz's newer scripts live in a sibling checkout; only collapse_blast_2.sh
# is taken from here (memory-efficient rewrite, byte-identical output). Override
# with NEWBIN=... if your copy lives elsewhere.
NEWBIN="${NEWBIN:-$REPO_ROOT/../hyb2_scripts}"

IN_SAM="${1:-$REPO_ROOT/../hyb2_baseline/test.sam}"
OUT="$REPO_ROOT/fixtures/sam_composition_run"

# explicit interpreters (see note 2 above)
PERL="${PERL:-perl}"
PYTHON="${PYTHON:-python3}"
AWK="${AWK:-awk}"

if [ ! -f "$IN_SAM" ]; then echo "Error: input SAM '$IN_SAM' not found" >&2; exit 1; fi
if [ ! -f "$NEWBIN/collapse_blast_2.sh" ]; then
    echo "Error: collapse_blast_2.sh not found in NEWBIN=$NEWBIN" >&2
    echo "       set NEWBIN=/path/to/hyb2_scripts" >&2; exit 1
fi

rm -rf "$OUT"; mkdir -p "$OUT"
export TMPDIR="$OUT/tmp"; mkdir -p "$TMPDIR"   # note 1: local scratch, not cluster path
cp "$IN_SAM" "$OUT/test.sam"
cd "$OUT"

echo "[1/8] sam2blast_3: SAM -> blast"
"$PYTHON" "$BIN/sam2blast_3" test.sam > test.blast

echo "[2/8] collapse_blast_2.sh -> collapse.blast ; mtophits_blast -> mtophits.blast"
# Split the shell's `collapse | mtophits` pipe so the collapse-only output is
# captured as its own golden (enables an independent collapse_blast parity test).
# mtophits reads the file exactly as it would read the pipe, so its output is
# unchanged vs the original one-liner.
bash "$NEWBIN/collapse_blast_2.sh" test.blast > test.collapse.blast
"$AWK" -f "$BIN/mtophits_blast" test.collapse.blast > test_mtophits.blast

echo "[3/8] create_reference_file.pl -> mtophits.ref"
"$PERL" "$BIN/create_reference_file.pl" test_mtophits.blast > test_mtophits.ref

echo "[4/8] get_mtop_hybrids.pl (call chimeras) -> test.hyb"
"$PERL" "$BIN/get_mtop_hybrids.pl" \
    BLAST_THRESHOLD=0.1 MODE=2 MAX_OVERLAP=4 MAX_HITS_PER_SEQUENCE=10 OUTPUT_FORMAT=HYB \
    test.blast > test.hyb

echo "[5/8] remove_duplicate_hybrids_hOH5_2.pl (dedup) -> test.ua.hyb"
"$PERL" -I"$BIN" "$BIN/remove_duplicate_hybrids_hOH5_2.pl" \
    PREFER_MIM=1 test_mtophits.ref test.hyb > test.ua.hyb

echo "[6/8] chimera gene-pair counts -> test.ua.hyb_stats_by_gene.txt"
cut -f4,10 test.ua.hyb | "$PERL" "$BIN/histogram.pl" > test.ua.hyb_stats_by_gene.txt

echo "[7/8] single-read top-hit counts -> test_tophit_by_gene.txt"
grep -v '^@' test.sam \
    | "$AWK" -F'\t' 'int($2/4)%2==0 && int($2/256)%2==0' \
    | cut -f3 | "$PERL" "$BIN/histogram.pl" > test_tophit_by_gene.txt

echo "[8/8] composition pies -> SKIPPED (hyb2_composition_pies.py not available; ask Grzegorz)"

rm -rf "$TMPDIR"
echo "Baseline regenerated under $OUT"
