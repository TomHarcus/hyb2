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

echo "[5/6] bp2hyb.sh: test.fragment_scores.txt -> test.ranked_interactions.txt (bp2hyb golden)"
bash "$REPO_ROOT/bin/bp2hyb.sh" < test.fragment_scores.txt > test.ranked_interactions.txt

echo "[6/6] comradesMakeConstraints_2 on a small window -> comrades_mini golden (end-to-end)"
# End-to-end golden for the comrades_make_constraints pipeline. Deliberately a
# SMALL input (20 chimeras): with so few base pairs the printed<=1000 fragment
# cap never truncates, so the histogram tie-order is irrelevant and the whole
# pipeline is deterministic -- giving a stable byte-exact golden. (On a dense
# window the cap would select 1000 rows in histogram order, which differs
# between the Perl and Python histograms -- see the ordering-fragility note.)
# The legacy script rm's its intermediates and writes them into cwd, so run it
# in a scratch dir and copy out only mini.hyb + the folding_constraints golden.
CMC_REF="$REPO_ROOT/data/Zika_18S_formatted.fasta"
mini_dir="$(mktemp -d)"
head -20 "$TIER2/test.ua.hyb" > "$mini_dir/mini.hyb"
cp "$CMC_REF" "$mini_dir/ref.fasta"
( cd "$mini_dir" && PATH="$REPO_ROOT/bin:$VIENNA_BIN:$PATH" LC_ALL=C \
    bash "$REPO_ROOT/bin/comradesMakeConstraints_2" -i mini.hyb -f ref.fasta -b 1 -e 10298 -r 1 >/dev/null 2>&1 )
cp "$mini_dir/mini.hyb" comrades_mini.hyb
cp "$mini_dir/mini.1-10298_folding_constraints.txt" comrades_mini.folding_constraints
rm -rf "$mini_dir"

echo "[7/7] make_VARNA_scores_2.sh -> make_varna.* (per-position colour-score golden)"
# make_VARNA_scores colours a folded structure's base pairs by their experimental
# support. Two inputs: basepair_scores (-i) and a .bps of the folded pairs (-b).
# The .bps must come from a CONSTRAINT-HONOURING fold, so it needs the Python
# comrades_fold (greedy loop) -- a one-shot RNAfold with all constraints at once
# is infeasible and folds nothing, and the legacy comradesFold2 is the buggy
# no-op (unconstrained), either of which gives a trivial all-zero golden.
# Window begin=10000 e=10300 so both arms of a real chimera fall in-fragment.
mv_dir="$(mktemp -d)"
head -20 "$TIER2/test.ua.hyb" > "$mv_dir/mini.hyb"
cp "$CMC_REF" "$mv_dir/ref.fasta"
(
  cd "$mv_dir"
  export PYTHONPATH="$REPO_ROOT/src" PATH="$REPO_ROOT/bin:$VIENNA_BIN:$PATH" LC_ALL=C
  "$VIENNA_BIN/python" -c "from hyb2.pipelines.comrades_make_constraints import run; run('mini.hyb','ref.fasta',10000,10300, vienna_bin='$VIENNA_BIN')" >/dev/null 2>&1
  "$VIENNA_BIN/python" -c "from hyb2.pipelines.comrades_fold import run; run('mini.10000-10300_folding_constraints.txt','ref_10000-10300.fasta', fold='vienna', vienna_bin='$VIENNA_BIN')" >/dev/null 2>&1
  # .bps from the constrained fold, genome coords (OFFSET = begin-1 = 9999)
  awk -v OFFSET=9999 'NR==1{print $1 "\t" $5}$1<$5 && NR>1{print $1+OFFSET "\t" $5+OFFSET}' ref_10000-10300.fasta.ct > mini.bps
  bash "$REPO_ROOT/bin/make_VARNA_scores_2.sh" -m 1000000 -l 10300 -t 301 -i mini.basepair_scores.txt -b mini.bps
)
cp "$mv_dir/mini.basepair_scores.txt" make_varna.basepair_scores.txt
cp "$mv_dir/mini.bps" make_varna.bps
cp "$mv_dir/mini__mini.VARNA_scores.txt" make_varna.VARNA_scores.golden
# the constrained-fold .ct is also the bp_score test input (bp_score extracts
# the same .bps from it via manifest offset=9999, len=10300, tail=301).
cp "$mv_dir/ref_10000-10300.fasta.ct" make_varna.ct
# comrades_fold golden: the same constrained fold, captured as a standalone
# regression snapshot (Option B, constraint-honouring != buggy legacy
# comradesFold2, so baselined FROM THE PORT -- Grzegorz-approved). Inputs =
# constraints + fragment fasta; goldens = the port's .vienna + .ct.
cp "$mv_dir/mini.10000-10300_folding_constraints.txt" comrades_fold.constraints
cp "$mv_dir/ref_10000-10300.fasta"                    comrades_fold.frag.fasta
cp "$mv_dir/ref_10000-10300.fasta.vienna"             comrades_fold.vienna
cp "$mv_dir/ref_10000-10300.fasta.ct"                 comrades_fold.ct
rm -rf "$mv_dir"

# svg_mod_coord oracle: VARNA renders make_varna.ct -> SVG, then the legacy
# svg_mod_coord.sh relabels it (branch 1 = single strand, branch 3 = two strand).
# Needs Java + the VARNA jar; skipped with a warning if either is missing.
VARNA_JAR="${VARNA_JAR:-$REPO_ROOT/VARNA/build/jar/VARNAcmd.jar}"
if command -v java >/dev/null 2>&1 && [ -f "$VARNA_JAR" ]; then
    echo "[8/8] VARNA -> svg ; svg_mod_coord.sh -> svg_mod.* goldens (branch 1 + 3)"
    svg_dir="$(mktemp -d)"
    cp make_varna.ct "$svg_dir/frag.ct"
    ( cd "$svg_dir"
      java -jar "$VARNA_JAR" -i frag.ct -bpStyle simple -spaceBetweenBases "0.6" -o frag.svg >/dev/null 2>&1
      bash "$REPO_ROOT/bin/svg_mod_coord.sh" -i frag.svg -x 100 && cp frag_plot.svg branch1
      rm -f frag_plot.svg
      bash "$REPO_ROOT/bin/svg_mod_coord.sh" -i frag.svg -x 100 -y 5000 -l 150 && cp frag_plot.svg branch3 )
    cp "$svg_dir/frag.svg" svg_mod.input.svg
    cp "$svg_dir/branch1" svg_mod.branch1.golden
    cp "$svg_dir/branch3" svg_mod.branch3.golden
    rm -rf "$svg_dir"

    echo "[9/9] plot_VARNA -> plot_varna.golden (headless: .ct -> VARNA -> svg_mod_coord)"
    pv_dir="$(mktemp -d)"
    cp make_varna.ct "$pv_dir/frag.ct"
    cp make_varna.VARNA_scores.txt "$pv_dir/s__frag.VARNA_scores.txt" 2>/dev/null \
        || cp make_varna.VARNA_scores.golden "$pv_dir/s__frag.VARNA_scores.txt"
    ( cd "$pv_dir" && PATH="$REPO_ROOT/bin:$VIENNA_BIN:$PATH" \
        bash "$REPO_ROOT/bin/plot_VARNA" -i frag.ct -s s__frag.VARNA_scores.txt \
             -j "$VARNA_JAR" -x 100 -l 150 >/dev/null 2>&1 )
    cp make_varna.ct plot_varna.ct
    cp make_varna.VARNA_scores.golden plot_varna.scores.txt 2>/dev/null || true
    cp "$pv_dir/s.frag_plot.svg" plot_varna.golden
    rm -rf "$pv_dir"
else
    echo "[8-9/9] SKIPPED svg_mod_coord + plot_VARNA oracles (java or VARNA jar not found)"
fi

echo "[10] hyb2_fold mode-1 (port) -> hyb2_fold.mode1.* goldens"
# hyb2_fold's transform (bin/hyb2_fold line 73) and fragment extraction (line 76)
# are legacy-faithful awk, so they are captured as PARITY goldens straight from
# the legacy awk. The fold is comrades_fold Option B (constraint-honouring,
# deliberately != the buggy legacy comradesFold2), so .vienna/.ct/log2scores are
# REGRESSION snapshots captured from the port (Grzegorz-approved re-baseline).
# Window: Zika_virusRNA 3900, len 150 (4 in-window chimeras -> full pipeline, fast fold).
HF_GENE=Zika_virusRNA; HF_X=3900; HF_L=150; HF_X2=$((HF_X + HF_L - 1))
# parity golden: legacy transform (line 73)
awk -v GENE=$HF_GENE -v X1=$HF_X -v X2=$HF_X2 -v X_COORD=$HF_X -v LEN=$HF_L \
  '{if ($4~GENE && $10~GENE && $7>=X1 && $8<=X2 && $13>=X1 && $14<=X2) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7-X_COORD+1"\t"$8-X_COORD+1"\t"$9"\t"$10"\t"$11"\t"12"\t"$13-X_COORD+1"\t"$14-X_COORD+1}' \
  "$TIER2/test.ua.hyb" > hyb2_fold.mode1.transform.golden
# parity golden: legacy fragment (line 76)
awk '{if(NR==1){print $0} else {if($0 ~ /^>/){print "\n"$0} else {printf $0}}}' "$CMC_REF" \
  | grep -A1 $HF_GENE \
  | awk -v b=$HF_X -v e=$HF_L '{if(/^>/){print $1}else{print substr($0,b,e)}}' > hyb2_fold.mode1.frag.golden
# port run -> Option B fold snapshots (+ plot svg if VARNA present)
hf_dir="$(mktemp -d)"
cp "$TIER2/test.ua.hyb" "$hf_dir/hf.hyb"
( cd "$hf_dir"
  export PYTHONPATH="$REPO_ROOT/src"
  "$VIENNA_BIN/python" -c "from hyb2.pipelines.hyb2_fold import run; run('hf.hyb','$HF_GENE',None,'$CMC_REF',$HF_X,None,$HF_L,'$VARNA_JAR',False,1,vienna_bin='$VIENNA_BIN')" >/dev/null 2>&1 || true )
cp "$hf_dir/$HF_GENE"_$HF_X-$HF_X2.fasta.vienna hyb2_fold.mode1.vienna
cp "$hf_dir/$HF_GENE"_$HF_X-$HF_X2.fasta.ct hyb2_fold.mode1.ct
cp "$hf_dir/hf__$HF_GENE"_$HF_X-$HF_X2.fasta.VARNA_log2scores.txt hyb2_fold.mode1.log2scores
cp "$hf_dir/hf.$HF_GENE"_$HF_X-$HF_X2.fasta_plot.svg hyb2_fold.mode1.plot.svg 2>/dev/null \
  || echo "  (no VARNA -> hyb2_fold.mode1.plot.svg skipped)"
rm -rf "$hf_dir"

echo "done -> $OUT"
ls -la "$OUT"
