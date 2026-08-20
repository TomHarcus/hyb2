#!/bin/sh
set -e
# GNU sort, parity-critical for collapse_blast
sort --version | grep -q "GNU coreutils" || { echo "MISSING: GNU sort"; exit 1; }
# every runtime tool resolves on PATH
for t in hyb2-py RNAfold bowtie2 hybrid-ss-min hybrid-min java Rscript; do
    command -v "$t" >/dev/null || { echo "MISSING: $t"; exit 1; }
done
# DESeq2 loads
Rscript -e 'suppressMessages(library(DESeq2))'
# HotKnots binary actually executes on this arch (126=not executable, 127=not found)
rc=0
/opt/hyb2/CPLfold/Utils/HotKnots_v2.0/bin/computeEnergy >/dev/null 2>&1 || rc=$?
if [ "$rc" -eq 126 ] || [ "$rc" -eq 127 ]; then
    echo "HotKnots computeEnergy did not execute (rc=$rc): arch/build problem"; exit 1
fi
echo "preflight OK"