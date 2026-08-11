# HYB2 Python port: running guide

This is the Python reimplementation of the HYB2 RNA proximity-ligation pipeline
(CLASH/PARIS/COMRADES/KARR-seq). It takes fastq/SAM input -> a chimeric-interaction
`.hyb` file -> contact-density maps, viewpoint graphs, and RNA secondary-structure
folding (ViennaRNA / UNAFold / **CPLfold**) rendered in VARNA.

It reproduces the legacy `bin/` pipeline stage-by-stage, validated by golden-diff tests.
R scripts (plotting, DESeq2) are unchanged and called as subprocesses.

---

## 1. Prerequisites

Everything runs inside a single conda env (called `hyb2` here). Tools needed:

| Tool | Provides | Install |
|---|---|---|
| Python ≥ 3.10 | the port (3.12 is what's tested) | (in the env) |
| ViennaRNA | `RNAfold`/`RNAcofold` (vienna backend, fixtures) | `conda install -c bioconda viennarna` |
| oligoarrayaux | `hybrid-ss-min`/`hybrid-min` (unafold backend) | `conda install -c bioconda oligoarrayaux` |
| bowtie2 | mapping fastq/fasta → SAM | `conda install -c bioconda bowtie2` |
| R + DESeq2 + ggplot2 + data.table | all plotting + `hyb2_compare` | `conda install -c bioconda -c conda-forge bioconductor-deseq2 r-ggplot2 r-data.table` |
| Java (JRE) + VARNA jar | structure rendering | `conda install -c conda-forge openjdk`; VARNA jar ships in `VARNA/` |
| numpy, pandas, pyyaml, tqdm, numba, scipy | the port + `hyb2_compare` + `--config` + progress bars + CPLfold | pulled by `pip install -e ".[cplfold]"` (below) |
| CPLfold + HotKnots | the CPLfold backend (local clone) | see §2 |

### Quick env build
```bash
conda create -n hyb2 python=3.12          # 3.10+ works; 3.12 is what's tested
conda activate hyb2
conda install -c bioconda -c conda-forge viennarna oligoarrayaux bowtie2 \
    bioconductor-deseq2 r-ggplot2 r-data.table openjdk
# from the repo root: installs the hyb2 package + console scripts + the
# Python deps (numpy, pandas, numba, scipy). Drop [cplfold] for numpy+pandas only.
pip install -e ".[cplfold]"
```

---

## 2. Per-machine setup (IMPORTANT: do these on each machine)

Three things are machine-specific and **must** be set up locally:

1. **CPLfold path.** `src/hyb2/config.py` has `CPLFOLD_DIR` hardcoded. Edit it to point at
   your CPLfold clone:
   ```python
   CPLFOLD_DIR = "/path/to/your/CPLfold"
   ```

2. **HotKnots is architecture-specific and must be compiled per machine.** CPLfold's energy
   binaries are C++; a binary built on one arch (e.g. Apple Silicon) will not run on another
   (e.g. x86-64 cluster). Rebuild after cloning / on a new machine:
   ```bash
   cd $CPLFOLD_DIR/Utils/HotKnots_v2.0
   make            # purge stale .o files first if switching arch
   ```
   If this is wrong, CPLfold folds return `energy = None` and the port raises a clear
   "rebuild HotKnots" error (it will not silently write a 0.0 energy).

3. **VARNA jar path (REQUIRED - machine-specific).** The port reads the jar location from
   the committed file `bin/VARNA.dir`, whose contents are an **absolute path from the
   original machine**, so after cloning it points somewhere that doesn't exist for you.
   Either edit `bin/VARNA.dir` so its single line is your jar path:
   ```
   /path/to/your/hyb2/VARNA/build/jar/VARNAcmd.jar
   ```
   or override it per-shell with `export HYB2_VARNA_JAR=/path/to/VARNAcmd.jar` (this takes
   precedence over the file). The jar itself ships in the repo at
   `VARNA/build/jar/VARNAcmd.jar`.

---

## 3. Environment gotchas (read before your first run)

- **Activate the env** (or put its `bin/` on PATH) before running. Every R step calls
  `Rscript` **bare**, so without the env active it grabs the *system* R (no DESeq2) and
  `hyb2_compare` fails. Same for bowtie2 / hybrid-* / java.
- **`UNAFOLDDAT`** (for the unafold backend) is set automatically by the port to
  `$CONDA_PREFIX/share/oligoarrayaux`. If you run `hybrid-ss-min` by hand it must be set, or
  it dies with "stack file is corrupt".
- **`hyb2_compare` input table:** use **non-numeric** dataset stems (e.g. `ctrl_rep1`, not
  `1`). R's `read.table` mangles numeric column headers and breaks DESeq2.

---

## 4. Commands

Console scripts (after `pip install -e .`): `hyb2-py`, `hyb2-coverage`, `plot-cdm`,
`hyb2-sam-composition`. The rest run via `python -m`:

| Command | Console script | or `python -m …` |
|---|---|---|
| main orchestrator | `hyb2-py` | `hyb2.cli` |
| coverage / CDM | `hyb2-coverage` | `hyb2.pipelines.hyb2_coverage` |
| folding | - | `hyb2.pipelines.hyb2_fold` |
| dataset comparison | - | `hyb2.pipelines.hyb2_compare` |
| mapping only | - | `hyb2.pipelines.bowtie2_map` |

Run `hyb2-py` with no args for full flag help.

---

## 5. Quick start: full pipeline from a SAM

```bash
conda activate hyb2
hyb2-py -i reads.sam -d reference.fasta -o myrun -a MyRNA -x 3900 -l 300 -r cplfold
```
Produces (prefix `myrun`): `myrun.hyb` (all chimeras) -> contact-density map PDF ->
viewpoint graph PDF -> CPLfold structure (`.ct`/`.vienna`) -> VARNA `_plot.svg`.

**Input types** (auto-detected from the extension): `.sam` (skips mapping), `.fastq` /
`.fastq.gz` / `.fasta` (mapped with bowtie2 first), or `.hyb` (skips straight to plotting/
folding). All five work.

Key flags: `-a` gene of interest, `-b` second gene, `-x`/`-y` fragment start coords,
`-l` fragment length, `-r` folding backend (`cplfold` default | `vienna` | `unafold`),
`-h` max hits/read, `-v` BLAST threshold, `-m` max overlap.

### The four analysis modes (as in the original README)
```bash
# short-range intramolecular
hyb2-py -i reads.sam -d ref.fasta -o test -a MyRNA -x 1001 -l 500
# long-range intramolecular
hyb2-py -i reads.sam -d ref.fasta -o test -a MyRNA -x 1001 -y 5001 -l 500
# intermolecular (two genes)
hyb2-py -i reads.sam -d ref.fasta -o test -a RNA_A -b RNA_B -x 7501 -y 501 -l 500
# homodimer
hyb2-py -i reads.sam -d ref.fasta -o test -a MyRNA -b MyRNA -x 3501 -y 3501 -l 200
```

### Reproducible runs with a YAML config (`--config`)

Both `hyb2-py` and `hyb2-fold` accept a YAML config of arguments, handy for
reproducible/shareable runs and CPLfold parameter sweeps. Templates ship in the repo:
`full_pipeline.yml` (full pipeline) and `standalone_fold.yml` (fold only).

```bash
cp full_pipeline.yml run.yml           # then edit paths/params
hyb2-py --config run.yml                  # run entirely from the config
hyb2-py --config run.yml --alpha 0.3      # CLI flags OVERRIDE the config
```
Precedence is **CLI flag > config value > built-in default**. Config keys are the
long-flag names (`input`, `reference`, `blast_threshold`, `x_start`, `alpha`, …); an
unknown key fails loudly. `hyb2-fold` has its own template with the same key names for
the overlapping args, so you can lift the fold section of a run straight across.

### Output verbosity (`-V`)

By default a run prints clean step-by-step status: `[1] Calling chimeras`, the folded
`change in G`, output files: with progress bars on the slow front-of-pipeline stages (interactive
terminals only; silent when piped/`tee`'d/on the cluster). Add **`-V/--verbose`** for the
full detail: the per-stage legacy messages plus the R/VARNA subprocess output. Useful for
debugging; leave it off for normal runs.

### Readable long flags

Every short flag has a descriptive long alias (`-i/--input`, `-d/--reference`,
`-a/--gene-1`, `-x/--x-start`, `-r/--fold-backend`, …). The short flags still work
(drop-in with the legacy), and the long names double as the config-file keys.

---

## 6. Reuse the `.hyb`: display any RNA's structure on demand

Chimeric calling is the expensive part, and the `.hyb` it produces contains chimeras for
**every** RNA in the reference at once. **Call chimeras once, then fold/plot any RNA
repeatedly without re-mapping or re-calling**, the orchestrator skips the whole
chimera-calling spine when the `.hyb` already exists.

```bash
# 1. Call chimeras once (all RNAs):
hyb2-py -i reads.sam -d ref.fasta -o test          # -> test.hyb

# 2. Then pick any RNA and fold/plot it - spine is skipped:
hyb2-py -i test.hyb -a RNA_A -x 100  -l 300
hyb2-py -i test.hyb -a RNA_B -x 500  -l 300
python -m hyb2.pipelines.hyb2_fold -i test.hyb -d ref.fasta -a RNA_C -x 900 -l 300 -r cplfold
```
Only coverage/viewpoint (CDM) needs no coords: `hyb2-py -i test.hyb -a RNA_A`.

---

## 7. Folding backends & CPLfold tuning

Choose the backend with `-r`: `cplfold` (default, finds pseudoknots), `vienna` (ViennaRNA),
`unafold` (UNAFold). **Both `hyb2-py` and `hyb2-fold` expose the full CPLfold parameter
surface** (`--alpha` etc.), so you can tune it from the one-shot pipeline or the standalone
fold command:

```bash
python -m hyb2.pipelines.hyb2_fold -i test.hyb -d ref.fasta -a MyRNA -x 3900 -l 300 \
    -r cplfold -p test_MyRNA_3900-4199.basepair_scores.txt \
    --alpha 0.5 --beta 0.0 --normalize log --beam-size 100
```
- `-p <basepair_scores>` turns the experimental bonus **on**; omit it for the no-bonus
  baseline (the fair A/B control).
- `--normalize raw|log`, `--alpha`, `--beta`, `--beam-size`, `--energy-delta`,
  `--max-phase1`, `--max-phase2`, `--energy-model` - all default from `config.CPL_DEFAULTS`.
- `-0 1` launches the interactive VARNA GUI (needs a display); omit for headless SVG output.

---

## 8. Comparing datasets (differential + similarity)

Compares ≥2 replicates per condition, runs DESeq2, produces a differential coverage map,
similarity heatmap, and enrichment tables. Needs a **manually-made** tab-delimited table
(non-numeric stems: see §3):
```
ctrl_rep1.hyb   ctrl_rep1.MyRNA.contact.txt   condition_one
ctrl_rep2.hyb   ctrl_rep2.MyRNA.contact.txt   condition_one
expt_rep1.hyb   expt_rep1.MyRNA.contact.txt   condition_two
expt_rep2.hyb   expt_rep2.MyRNA.contact.txt   condition_two
```
```bash
python -m hyb2.pipelines.hyb2_compare -i input.table -o cmp -a MyRNA -d ref.fasta
```

---

## 9. Testing on large files (the important bit for scale testing)

The pipeline has been validated for **correctness** on the Zika `testData.sam`, but **not
yet for scale on real large data**. The front-of-pipeline stages (where the data is biggest)
now stream, **except `collapse_blast`, the one remaining bottleneck:**

| Stage | Status |
|---|---|
| `sam2blast` | **streams (O(1) memory), ~3.3× faster**: byte-identical output |
| `bowtie2_map` fastq/fastq.gz prep | **streams (O(unique reads))**: byte-identical |
| `mtophits_blast` | **streams (O(unique read IDs))**: byte-identical |
| `collapse_blast` | **still materializes the blast ~2.7× in RAM (≈1.5× the SAM), linear**: the bottleneck |

**Empirically: a real ~50 GB SAM OOM-kills at `collapse_blast`.** `sam2blast` completes
(26 GB blast, ~11 min), then `collapse` is killed loading that blast into RAM (needs ~75 GB;
died at ~9 GB on a 15 GB-RAM box). This is a regression vs the legacy `collapse_blast_2.sh`,
whose `sort -k13` was disk-backed/memory-safe, so the faithful fix restores that. **The
pipeline cannot complete a 50 GB SAM until `collapse_blast` is fixed** (external-sort vs
clean-rewrite: decision pending).

Everything **downstream of the `.hyb`** (folding, coverage, compare) operates on
already-reduced data and is not a large-file concern.

To profile a run's peak memory:
```bash
/usr/bin/time -v hyb2-py -i big.sam -d ref.fasta -o big -a MyRNA 2>&1 | grep "Maximum resident"
```

---

## 10. Run the test suite

```bash
conda activate hyb2
PYTHONPATH=src python -m pytest tests/ -q
```
Some tests are **gated**, they skip (not fail) if a tool or a regenerated fixture is
missing. Fixtures are gitignored and rebuilt per-machine:
```bash
scripts/generate_sam_composition_baseline.sh   # primary spine oracle
scripts/generate_coverage_baseline.sh
scripts/generate_viewpoint_baseline.sh
scripts/generate_folding_baseline.sh           # needs ViennaRNA + java
scripts/generate_unafold_baseline.sh           # needs oligoarrayaux
```

---

## 11. Known limitations

- **Scale untested on real large data** (see §9): the reason for this hand-off.
- **`comradesScore`** (randomized 1000× parallel folding, significance scoring): not
  ported; needs a `qsub` cluster.
- **`hyb2_app`** (Shiny GUI): the R app is unchanged; only its thin launcher wrapper is
  not yet ported.
- **No per-stage checkpointing yet:** the `.hyb` is reused if present (§6), but if the
  `.hyb` is missing the spine regenerates the blast intermediates from scratch (no
  skip-if-`test.blast`-exists yet).
