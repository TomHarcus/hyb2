# HYB2 Python port: running guide

This is the Python reimplementation of the HYB2 RNA proximity-ligation pipeline
(CLASH/PARIS/COMRADES/KARR-seq). It takes fastq/SAM input -> a chimeric-interaction
`.hyb` file -> contact-density maps, viewpoint graphs, and RNA secondary-structure
folding (ViennaRNA / UNAFold / **CPLfold**) rendered in VARNA.

It reproduces the legacy `legacy_bin/` pipeline stage-by-stage, validated by golden-diff tests.
R scripts (plotting, DESeq2) are unchanged and called as subprocesses.

---

## Prerequisites

You need **one** container runtime (matched to where you will run) and **nothing else** (no conda, R, viennaRNA, or per-tool installs. Its 
all within the image):

- **Cluster (Eddie / HPC):** **Apptainer**: usually available as a module (`module load apptainer`), so nothing to install (find it with `module avail apptainer`).
- **Local Linux / WSL / macOS:** **Docker** install via the official guide: 
<https://docs.docker.com/get-docker/>.
- **Apptainer on your own Linux machine** (not a cluster): see
  <https://apptainer.org/docs/admin/main/installation.html>
  (recent Debian/Ubuntu: `sudo apt install -y apptainer`).

You do *not* need both: pick the runtime for your machine. (Building the image yourself,
rather than pulling it, additionally needs Docker: see *Developing*.)

## 1. Get the container (no install)

The pipeline ships as a container image on GHCR with everything needed already there: Python, the conda env, R + DESeq2
, ViennaRNA, CPLfold + HotKnots, VARNA, bowtie2, Java. **No conda, no compiling, no per-machine install.**

**Linux cluster (e.g. Eddie): Use Apptainer:**

```bash
module load apptainer   # use apptainer/1.4.4 on Eddie
apptainer pull docker://ghcr.io/tomharcus/hyb2:latest   # produces hyb2_latest.sif
apptainer run --bind "$TMPDIR" hyb2_latest.sif hyb2-py --config run.yml
```

**Local Linux / WSL / macOS: Use Docker:**

```bash
docker pull ghcr.io/tomharcus/hyb2:latest

# -t shows progress bars when supported, --user keeps outputs owned by you
docker run --rm -t \
      --user "$(id -u):$(id -g)" -e HOME=/tmp \
      -v "$PWD:/data" -w /data \
      ghcr.io/tomharcus/hyb2:latest hyb2-py --config run.yml
```
> On macOS (Docker Desktop), omit `--user ... -e HOME=/tmp`. Ownership is mapped to you automatically and `--user`
> can cause permission errors there. Drop `-t` when piping output to a file (it errors without a terminal).

(Apple silicon should emulate running the image. This is fine for testing but slow for big mapping. The cluster is recommended for real runs)

`:latest` tracks the most up to date build. If you want a specific version: instead of `:latest` add `:<git-sha>`.
To pull a new build, re-pull (`apptainer pull --force... / docker pull...`)

**To run the pipeline with a simple command see the convenience part of section 2.**

## For developing / modifying the code

You only need this if you are **changing the pipeline itself**. To just run it, use just the container, no clone or setup required.

```bash
git clone -b python-migration https://github.com/TomHarcus/hyb2.git
cd hyb2
```

### Quick edits (recommended)

The image already contains the full environment, so bind-mount your source over it and the edits are live with no rebuild. The
package is installed with `pip install -e`, so it reads straight from the mounted source. Mount `src/` for Python and `rscripts/` for 
the R plotting scripts:

```bash
docker run --rm \
    -v "$PWD/src:/opt/hyb2/src" \
    -v "$PWD/rscripts:/opt/hyb2/rscripts" \
    -v "$PWD/mydata:/data" -w /data \
    ghcr.io/tomharcus/hyb2:latest hyb2-py --config run.yml
```

Make the changes, then re-run and the pipeline changes right away.

### Full native env for the test suite or bigger work

Recreate the environment the image is built from (the `Dockerfile` contains the setup directions). Follow its `conda env create`, 
CPLfold clone, and HotKnots `make` steps.

```bash
conda env create -f environment.yml && conda activate hyb2
# then CPLfold + HotKnots exactly as the Dockerfile does (clone, purge *.o/*.a, make),
# and export UNAFOLDDAT / HYB2_VARNA_JAR / HYB2_CPLFOLD_DIR as its ENV block sets them
PYTHONPATH=src python -m pytest tests/ -q   # run the tests
```

### How a change works

edit files -> test (bind-mount or native) -> git commit + push
           -> CI rebuilds the image -> publishes ghcr.io/tomharcus/hyb2:latest
           -> re-pull on the cluster / local machine to run the new version

Any changes made to `src/`, `rscripts/`, the `Dockerfile`, `environment.yml`, etc and pushed trigger the CI automatically.
So `:latest` points to that newest build.

To build the full image locally before pushing (catches build breaks without waiting for the CI):

```bash
docker_scripts/build.sh build   # full image build
docker_scripts/build.sh smoke   # entrypoint + tools resolve
```

---

## 2. Commands

The pipeline exposes four commands. Examples throughout the `README` show the bare command (e.g. `hyb2-py` etc). How you actually
run them depends on how you are running:

- **Container: cluster (Apptainer):**
  `apptainer run --bind "$TMPDIR" --bind /path/to/scratch hyb2_latest.sif <command>`

- **Container: local (Docker):**
  `docker run --rm -t  --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD:/data" -w /data ghcr.io/tomharcus/hyb2:latest <command>`

- **Dev clone (native env):** the commands are console scripts: run `<command>` directly, or use the 
`python -m ...` module form.

| Command | Console script | or `python -m …` |
|---|---|---|
| main orchestrator | `hyb2-py` | `hyb2.cli` |
| folding | `hyb2-fold` | `hyb2.folding.hyb2_fold` |
| coverage / CDM | `hyb2-coverage` | `hyb2.coverage.hyb2_coverage` |
| dataset comparison | `hyb2-compare` | `hyb2.compare.hyb2_compare` |

Run any command with no args (or `--help`) for full flag help — e.g. `hyb2-py`.

### Convenience: a shell wrapper

The full `docker run ...` / `apptainer run ...` line is long. You can define a wrapper **once** in your shell's startup file: `~/.bashrc` on Linux/WSL (default bash) or
`~/.zshrc` on macOS (default zsh). Check your shell with `echo $SHELL`. This makes every command short.
Pick the one for your platform:

**Linux / WSL (Docker):**
```bash
hyb2() {
  local tty=""; [ -t 1 ] && tty="-t"    # progress bars when interactive, safe when piped
  docker run --rm $tty --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$PWD:/data" -w /data ghcr.io/tomharcus/hyb2:latest "$@"
}
```

**macOS (Docker):** omit `--user/HOME` (Docker desktop maps ownership automatically):
```bash
hyb2() {
  local t=""; [ -t 1 ] && t="-t" 
  docker run --rm $t -v "$PWD:/data" -w /data ghcr.io/tomharcus/hyb2:latest "$@"; 
}
```

**Cluster (Apptainer):**
```bash
hyb2() { 
  apptainer run --bind "$TMPDIR" --bind /exports/eddie/scratch/$USER ~/hyb2_latest.sif "$@"; 
}
```

Afterwards, run `source ~/.bashrc` or `source ~/.zshrc` or open a new shell and run from your data directory:
```bash
cd my_data_dir
hyb2 hyb2-py --config run.yml
hyb2 hyb2-fold -i myrun.hyb -d ref.fasta -a MyRNA -x 108 -l 1168 -r cplfold
```

The wrapper mounts your current directory as `/data`, so `cd` into your data folder first
and reference files by plain name (or `/data/...`). It handles the progress bars, file ownership,
and required binds for you.

**IMPORTANT: The rest of the commands assume the wrapper has been set up.**

---

## 3. Running Hyb2 on a cluster (e.g. Eddie)

The container runs on any cluster with Apptainer, but after testing there are a few Eddie specific
rules that make it run.

1. **Run from scratch, not home.** Home's ~10 GB quota can't hold the pipeline's intermediate files
(SAM/blast can be 10s of GB). Work in `/exports/eddie/scratch/$USER/...`.

2. **Request cores**: when using a small amount of cores bowtie2 crawls and seems to hang. Use
`#$ -pe sharedmem 16` (batch) or `qlogin -pe sharedmem N` (interactive). The pipeline maps with
exactly your allocated cores.

3. **`h_vmem >= 16G`**: needed both for mapping and for the `apptainer pull` squashfs conversion
(it OOMs on 8G).

4. **Keep Apptainer's cache/tmp off home:**
```bash
export APPTAINER_CACHEDIR=/exports/eddie/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/exports/eddie/scratch/$USER/apptainer_tmp
```

5. **Run each job in a clean/empty dir**: a leftover `.tab`/index from a prior run makes the pipeline
skip the database build (bowtie2 then errors "index does not exist").

6. **Bind scratch + `$TMPDIR`**: Apptainer auto mounts `$HOME` but not scratch or the node-local
`$TMPDIR` (where the big sort spills).

Batch script (run_hyb2.sh, submit from your scratch run dir):
```bash
#!/bin/bash
#$ -cwd
#$ -N hyb2
#$ -pe sharedmem 16
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00

. /etc/profile.d/modules.sh
module load apptainer/1.4.4
hyb2 hyb2-py --config run.yml
```

```bash
qsub run_hyb2.sh    # queues + runs when cores free up. Check status with qstat
```

For a quick interactive test: `qlogin -pe sharedmem 4 -l h_vmem=16G`, then
`module load apptainer/1.4.4` and run directly.

---

## 4. Quick start: full pipeline from a SAM

Full pipeline in one command (cluster/Apptainer shown. See section 2 for Docker/native):

```bash
hyb2 hyb2-py -i reads.sam -d reference.fasta -o myrun -a MyRNA -x 3900 -l 300 -r cplfold
```

Produces (prefix `myrun`): `myrun.hyb` (all chimeras) -> contact-density map PDF ->
viewpoint graph PDF -> CPLfold structure (`.ct`/`.vienna`) -> VARNA `_plot.svg`.

**Input types** (auto-detected from the extension): `.sam` (skips mapping), `.fastq` /
`.fastq.gz` / `.fasta` (mapped with bowtie2 first), or `.hyb` (skips straight to plotting/
folding). 

Key flags: `-a` gene of interest, `-b` second gene, `-x`/`-y` fragment start coords,
`-l` fragment length, `-r` folding backend (`cplfold` default | `vienna` | `unafold`),
`-h` max hits/read, `-v` BLAST threshold, `-m` max overlap.

### The four analysis modes (as in the original README)
```bash
# short-range intramolecular
hyb2 hyb2-py -i reads.sam -d ref.fasta -o test -a MyRNA -x 1001 -l 500
# long-range intramolecular
hyb2 hyb2-py -i reads.sam -d ref.fasta -o test -a MyRNA -x 1001 -y 5001 -l 500
# intermolecular (two genes)
hyb2 hyb2-py -i reads.sam -d ref.fasta -o test -a RNA_A -b RNA_B -x 7501 -y 501 -l 500
# homodimer
hyb2 hyb2-py -i reads.sam -d ref.fasta -o test -a MyRNA -b MyRNA -x 3501 -y 3501 -l 200
```

### Reproducible runs with a YAML config (`--config`)

All `hyb2-py`, `hyb2-fold`, `hyb2-coverage`, and `hyb2-compare` accept a YAML config of arguments, handy for
reproducible/shareable runs and CPLfold parameter sweeps. To get the templates run:
```bash
# container user: fetch a template from GitHub:
curl -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/pipeline_templates/full_pipeline.yml
mv full_pipeline.yml run.yml        # then edit paths/params

# dev clone: copy it locally instead:
cp pipeline_templates/full_pipeline.yml run.yml
```

In `pipeline_templates/` you will find:
`full_pipeline.yml` (full pipeline), `folding.yml` (fold only), `compare.yml` (comparison only), and `coverage.yml` (coverage only).

```bash
cp pipeline_templates/full_pipeline.yml run.yml           # then edit paths/params
hyb2 hyb2-py --config run.yml                  # run entirely from the config
hyb2 hyb2-py --config run.yml --alpha 0.3      # CLI flags OVERRIDE the config
```
Precedence is **CLI flag > config value > built-in default**. Config keys are the
long-flag names (`input`, `reference`, `blast_threshold`, `x_start`, `alpha`, …); an
unknown key fails loudly. Leave `varna_jar` **unset/commented** to auto-resolve the bundled
jar, only set it to point at a VARNA jar in a non-default location.

### Deterministic runs (`--reproducible`)

By default bowtie2 maps with multiple threads and emits reads in **thread-completion order**, which varies run-to-run. The alignments are identical, only their *order* differs, but the order-sensitive `collapse` / `mtophits` stages turn that into a different `.hyb`, so two runs of the same data produce byte-different intermediates. The resulting folded structures and contact maps seem to be unaffected, only the bytes differ.

Pass **`--reproducible`** (or `reproducible: true` in the config) to make bowtie2 emit reads in input order (`--reorder`), so the whole pipeline is byte-deterministic run-to-run. This slows down the mapping stage slightly, so by default it is off.

```bash
hyb2 hyb2-py --config run.yml --reproducible
```

PDF's still byte differ, even with `--reproducible`. The PDF outputs (contact maps, viewpoint graphs) won't match byte-for-byte between runs, because R's `pdf()` device embeds a creation timestamp. The plots are identical, only the metadata differs. To verify two runs match, compare without the PDFs:

```bash
diff -rq run_1 run_2 --exclude='*.pdf'
```

### Output verbosity (`-V`)

By default a run prints clean step-by-step status. On an interactive terminal you get live progress bars/spinners on the slow stages of the pipeline. When piped or on the cluster (no terminal)
those same steps print as `stage ... / done (Xs)` lines with elapsed time, so cluster logs stay readable and show per-stage timing. Add `-V/--verbose` for full details (per-stage legacy messages + R/VARNA subprocess output).



### Readable long flags

Every short flag has a descriptive long alias (`-i/--input`, `-d/--reference`,
`-a/--gene-1`, `-x/--x-start`, `-r/--fold-backend`, …). The short flags still work
(drop-in with the legacy), and the long names double as the config-file keys.

---

## 5. Reuse the `.hyb`: display any RNA's structure on demand

Chimeric calling is the expensive part, and the `.hyb` it produces contains chimeras for
**every** RNA in the reference at once. **Call chimeras once, then fold/plot any RNA
repeatedly without re-mapping or re-calling**, when you pass the `.hyb` as input (`-i test.hyb`)

```bash
# 1. Call chimeras once (all RNAs):
hyb2 hyb2-py -i reads.sam -d ref.fasta -o test          # -> test.hyb

# 2. Then pick any RNA and fold/plot it - spine is skipped:
hyb2 hyb2-py -i test.hyb -a RNA_A -x 100  -l 300
hyb2 hyb2-py -i test.hyb -a RNA_B -x 500  -l 300
hyb2 hyb2-fold -i test.hyb -d ref.fasta -a RNA_C -x 900 -l 300 -r cplfold
```
Only coverage/viewpoint (CDM) needs no coords: `hyb2 hyb2-py -i test.hyb -a RNA_A`.

---

## 6. Folding backends & CPLfold tuning

Choose the backend with `-r`: `cplfold` (default, finds pseudoknots), `vienna` (ViennaRNA),
`unafold` (UNAFold). **Both `hyb2-py` and `hyb2-fold` expose the full CPLfold parameter
surface** (`--alpha` etc.), so you can tune it from the one-shot pipeline or the standalone
fold command:

> **`-r unafold` and constraints.** The pipeline bundles **OligoArrayAux** (`hybrid-ss-min`),
> the freely-redistributable subset of UNAFold (the *full* UNAFold is licensed and cannot be
> shipped in the conda env). OligoArrayAux's `hybrid-ss-min` **silently ignores `--force`
> constraints**, so `-r unafold` folds the fragment **unconstrained** (it does not incorporate the
> experimental base-pair support) and no error is raised. To get constraint-honouring UNAFold folds
> you must install the full licensed UNAFold-3.8 and point the tools at its `hybrid-ss-min`. For
> constraint-guided folding, use **`cplfold` (default)** or **`vienna`**, which both honour the
> experimental constraints.

```bash
hyb2 hyb2-fold -i test.hyb -d ref.fasta -a MyRNA -x 3900 -l 300 \
    -r cplfold -p test_MyRNA_3900-4199.basepair_scores.txt \
    --alpha 0.5 --beta 0.0 --normalize log --beam-size 100
```
- `-p <basepair_scores>` turns the experimental bonus **on**; omit it for the no-bonus
  baseline (the fair A/B control).
- `--normalize raw|log`, `--alpha`, `--beta`, `--beam-size`, `--energy-delta`,
  `--max-phase1`, `--max-phase2`, `--energy-model` - all default from `config.CPL_DEFAULTS`.
- `-0 1` launches the interactive VARNA GUI (needs x11, else use the SVG); omit for headless SVG output.

---

## 7. Comparing datasets (differential + similarity)

Compares ≥2 replicates per condition, runs DESeq2, produces a differential coverage map,
similarity heatmap, and enrichment tables. Needs a **manually-made** tab-delimited table.

> **Gotcha:** use **non-numeric** dataset stems (e.g. `ctrl_rep1`, not `1`). R's
> `read.table` mangles numeric column headers and breaks DESeq2.

```
ctrl_rep1.hyb   ctrl_rep1.MyRNA.contact.txt   condition_one
ctrl_rep2.hyb   ctrl_rep2.MyRNA.contact.txt   condition_one
expt_rep1.hyb   expt_rep1.MyRNA.contact.txt   condition_two
expt_rep2.hyb   expt_rep2.MyRNA.contact.txt   condition_two
```
```bash
hyb2 hyb2-compare -i input.table -o cmp -a MyRNA -d ref.fasta
```

---

## 8. Large files & memory

The front-of-pipeline stages (where the data is biggest) all stream or are memory-safe:

| Stage | Status |
|---|---|
| `sam2blast` | streams (O(1) memory), ~3.3× faster; byte-identical output |
| `bowtie2_map` fastq/fastq.gz prep | streams (O(unique reads)); byte-identical |
| `mtophits_blast` | streams (O(unique read IDs)); byte-identical |
| `collapse_blast` | **disk-backed external sort** (memory bounded by the sort buffer); byte-identical |

`collapse_blast` uses a disk-backed `sort` (memory capped by `-S`, spilling to `TMPDIR`),
mirroring the legacy `collapse_blast_2.sh`. A real ~50 GB SAM that previously OOM-killed
now completes.

> **`TMPDIR` must be on real disk.** A tmpfs (RAM-backed) `/tmp` (common on Linux)
> re-introduces the OOM, because the sort spills into RAM. The **container's entrypoint**
> applies this guard automatically: if `TMPDIR` is unset or points at tmpfs it falls back
> to `$HOME/scratch_tmp` (real disk). On a cluster, set `TMPDIR` (and `APPTAINER_TMPDIR`)
> to scratch and bind it in (see section 3), so both the sort and Apptainer's own tmp land on
> large real disk.

Everything **downstream of the `.hyb`** (folding, coverage, compare) operates on
already-reduced data and is not a large-file concern.

To profile a run's peak memory:
```bash
/usr/bin/time -v hyb2 hyb2-py -i big.sam -d ref.fasta -o big -a MyRNA 2>&1 | grep "Maximum resident"
```

---

## 9. Run the test suite

The tests run against the source in a **dev clone**, not the container. `tests/` isn't
shipped in the image. Set up the native env or a bind-mount first (see *Developing /
modifying the code*), then:

```bash
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

## 10. Known limitations

- **`comradesScore`** (randomized 1000× parallel folding, significance scoring): not
  ported; needs a `qsub` cluster.
- **`hyb2_app`** (Shiny GUI): the R app is unchanged; only its thin launcher wrapper is
  not yet ported.

