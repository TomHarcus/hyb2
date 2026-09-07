#!/bin/sh
set -e
# GNU sort, parity-critical for collapse_blast
sort --version | grep -q "GNU coreutils" || { echo "MISSING: GNU sort"; exit 1; }
# every runtime tool resolves on PATH
for t in hyb2 RNAfold bowtie2 hybrid-ss-min hybrid-min java Rscript; do
    command -v "$t" >/dev/null || { echo "MISSING: $t"; exit 1; }
done
# DESeq2 loads
Rscript -e 'suppressMessages(library(DESeq2))'

echo "preflight OK"