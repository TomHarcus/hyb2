#!/bin/bash -e
# Regenerates fixtures/legacy_run/ -- the golden-output regression harness used
# to validate each Python-ported stage against the current Perl/Bash/AWK
# pipeline in bin/. Run this from a clean checkout whenever you need a fresh
# baseline (e.g. after the legacy scripts change, or on a different machine).
#
# Requires the "hyb2" conda environment (conda env create -f bin/hyb2.yml).

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURES="$REPO_ROOT/fixtures/legacy_run"
DATA="$REPO_ROOT/data"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hyb2

rm -rf "$FIXTURES"
mkdir -p "$FIXTURES"

# --- Stage 1+2: full hyb2 orchestrator, SAM -> hyb -> contact density maps ---
# NOTE: as of this writing, bin/hyb2 passes -1/-2 to plot_viewpoint and
# hyb2_fold, but those scripts' getopts strings only accept -d/-b -- so the
# folding and viewpoint steps silently no-op ("incorrect option", exit 0)
# when run through this top-level entry point. This run captures that
# real, current behavior; see fixtures/README.md.
mkdir -p "$FIXTURES/main_short_range"
(
    cd "$FIXTURES/main_short_range"
    hyb2 -i "$DATA/testData.sam" -d "$DATA/Zika_18S.fasta" -o test -a Zika_virusRNA -x 1001 -l 500 \
        > run.stdout.log 2> run.stderr.log || true
)

# --- Stage 3: folding, invoked directly (this is the CoupleFold integration point) ---
# NOTE: must pass the underscore-formatted fasta (Zika_18S_formatted.fasta),
# not the pipe-delimited Zika_18S.fasta the README shows -- the gene-ID grep
# in hyb2_fold's fasta extraction won't match "Zika|virusRNA" against
# "Zika_virusRNA". See fixtures/README.md.
mkdir -p "$FIXTURES/direct_fold"
(
    cd "$FIXTURES/direct_fold"
    cp "$FIXTURES/main_short_range/test.hyb" .
    hyb2_fold -i test.hyb -d "$DATA/Zika_18S_formatted.fasta" -a Zika_virusRNA -x 1001 -l 500 \
        > run.stdout.log 2> run.stderr.log || true
)

# --- Stage 2b: viewpoint graph, invoked directly ---
mkdir -p "$FIXTURES/direct_viewpoint"
(
    cd "$FIXTURES/direct_viewpoint"
    cp "$FIXTURES/main_short_range/test.hyb" .
    plot_viewpoint -i test.hyb -d "$DATA/Zika_18S_formatted.fasta" -a Zika_virusRNA \
        > run.stdout.log 2> run.stderr.log || true
)

echo "Fixtures regenerated under $FIXTURES"
