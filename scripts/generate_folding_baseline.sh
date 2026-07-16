#!/bin/bash -e
# Regenerates fixtures/folding_run/ -- the real RNA-folding oracle that the
# remaining Tier 2 parser stages need (ct2bps_2, Ct2B_GK_3, make_nicer_vienna,
# add_dG_hyb). Until this existed there was no captured .ct/.vienna output to
# diff those ports against ("No folding oracle yet" in CLAUDE.md's open items).
#
# It reproduces the folding step of comradesMakeConstraints_2 VERBATIM (that
# script's line 57, the -r 1 / ViennaRNA branch):
#
#   paste -d "&" bit_1.fasta bit_2.fasta \
#     | RNAcofold --noconv --noPS \
#     | sed 's/&>/-/g;s/&//g' \                 <- .vienna  (make_nicer_vienna, add_dG input)
#     | b2ct | sed 's/ENERGY =/dG =/g' > .ct    <- .ct      (ct2bps_2, Ct2B_GK_3 input)
#
# The paired-bit FASTAs (bit_1/bit_2) come from fixtures/tier2_run/, produced by
# hyb2fasta_bits_allRNAs on the real test.ua.hyb -- so this is genuine ViennaRNA
# output on the real Zika/18S test data, not a synthetic structure.
#
# Requires ViennaRNA (RNAfold/RNAcofold/b2ct) on PATH, or point at the hyb2
# conda env's bin with VIENNA_BIN=... (default below).
#
# Usage:  scripts/generate_folding_baseline.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TIER2="$REPO_ROOT/fixtures/tier2_run"
OUT="$REPO_ROOT/fixtures/folding_run"

VIENNA_BIN="${VIENNA_BIN:-$HOME/miniconda3/envs/hyb2/bin}"
RNACOFOLD="$VIENNA_BIN/RNAcofold"
B2CT="$VIENNA_BIN/b2ct"

for tool in "$RNACOFOLD" "$B2CT"; do
    if [ ! -x "$tool" ]; then
        echo "Error: $tool not found/executable." >&2
        echo "       Install ViennaRNA (conda install -n hyb2 -c bioconda viennarna)" >&2
        echo "       or set VIENNA_BIN=/path/to/env/bin" >&2
        exit 1
    fi
done
for f in test.ua.bit_1.fasta test.ua.bit_2.fasta; do
    if [ ! -f "$TIER2/$f" ]; then
        echo "Error: $TIER2/$f not found; run scripts/generate_tier2_baseline.sh first" >&2
        exit 1
    fi
done

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$OUT"

echo "[1/2] RNAcofold: $(grep -c '^>' "$TIER2/test.ua.bit_1.fasta") bit pairs -> test.ua.cofold.vienna"
paste -d "&" "$TIER2/test.ua.bit_1.fasta" "$TIER2/test.ua.bit_2.fasta" \
    | "$RNACOFOLD" --noconv --noPS \
    | sed 's/&>/-/g;s/&//g' > test.ua.cofold.vienna

echo "[2/3] b2ct: test.ua.cofold.vienna -> test.ua.cofold.ct"
"$B2CT" < test.ua.cofold.vienna | sed 's/ENERGY =/dG =/g' > test.ua.cofold.ct

echo "[3/5] ct2bps_2.awk: test.ua.cofold.ct -> test.ua.cofold.bps (parser golden)"
# Golden for the ct2bps_2 port: legacy awk output on the real .ct above.
awk -f "$REPO_ROOT/bin/ct2bps_2.awk" test.ua.cofold.ct > test.ua.cofold.bps

# Everything below reproduces the bp2hyb feed from comradesMakeConstraints_2:
#   ct2bps | histogram > basepair_scores           (line 61)
#   awk 'printed<=1000 && in-window' > fragment     (line 64)  <- caps the input
#   bp2hyb.sh < fragment > ranked_interactions      (line 65)
# The cap matters: bp2hyb reformats each base pair into an all-"RNA" chimera, so
# every record lands in a single combine bucket -> combine is O(n^2). The real
# pipeline never feeds it more than ~1000 lines; feeding the full ~172k basepair
# scores makes the legacy Perl combine effectively never finish. We also pre-sort
# the scores so the frozen fragment fixture is reproducible (bp2hyb re-sorts
# internally, so this does not change its output).
PERL_BIN="${PERL:-perl}"
export PATH="$REPO_ROOT/bin:$PATH" LC_ALL=C   # bp2hyb.sh shells out to combine_hyb_merge_touching.pl

echo "[4/5] histogram + fragment cap -> test.fragment_scores.txt (bp2hyb input, ~1000 lines)"
awk -f "$REPO_ROOT/bin/ct2bps_2.awk" test.ua.cofold.ct \
    | "$PERL_BIN" "$REPO_ROOT/bin/histogram.pl" \
    | sort -k1,1n -k2,2n \
    | awk 'printed<=1000 && $1>=1 && $1<=10298 && $2>=1 && $2<=10298{print;printed++}' \
    > test.fragment_scores.txt

echo "[5/5] bp2hyb.sh: test.fragment_scores.txt -> test.ranked_interactions.txt (bp2hyb golden)"
bash "$REPO_ROOT/bin/bp2hyb.sh" < test.fragment_scores.txt > test.ranked_interactions.txt

echo "done -> $OUT"
ls -la "$OUT"
