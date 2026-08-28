# HYB2 Python port: running guide

This is the Python reimplementation of the HYB2 RNA pipeline. It is used for analyzing RNA proximity ligation experiments from mapped files in the fastq/SAM format to:
1. Generate a list of chimeric interactions with their coordinates, sequence, and folding energy
2. Plot contact density maps of selected genes and viewpoint graphs
3. Generate intra/intermolecular RNA structures of selected regions

## Table of Contents
- [Prerequisites](#prerequisites)
- [Local Linux/WSL and macOS Installation](#local-linuxwsl-and-macos-installation)
- [Eddie Installation](#eddie-installation)
- [Running hyb2](#running-hyb2)
- [Running The Pipeline on Eddie with Batch Jobs](#running-the-pipeline-on-eddie-with-batch-jobs)



## Prerequisites

### Install Docker (Ignore if wanting to run hyb2 on Eddie)

#### Linux/WSL:
Install Docker Engine directly by following the official guide for your distribution:
- Ubuntu: https://docs.docker.com/engine/install/ubuntu/
- Fedora: https://docs.docker.com/engine/install/fedora/
- Debian: https://docs.docker.com/engine/install/debian/
- Other distros: https://docs.docker.com/engine/install/

#### macOS:
On macOS install Docker Desktop:
1. Go to https://www.docker.com/products/docker-desktop/
2. Download the compatible build (Apple Silicon or Intel)
3. Run it and follow the prompts (no need to create an account or sign in)
4. Verify it worked by opening a terminal and running:
```bash
docker run hello-world
```

## Local Linux/WSL and macOS Installation

### If using Linux/WSL:
Start off by installing Docker Engine by following the Linux/WSL section in [Prerequisites](#prerequisites)

### If using macOS:
Start off by installing Docker Desktop by following the macOS section in [Prerequisites](#prerequisites)


Once Docker is installed you can now grab the hyb2 container by running:
```bash
docker pull ghcr.io/tomharcus/hyb2:latest
```
This command fetches the most up to date version of hyb2.

Next, to fetch the hyb2 setup script run:
```bash
curl -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/setup_hyb2.sh
```
This command downloads the script into your current directory. It is named: `setup_hyb2.sh`

To run the setup script, run these commands:
```bash
chmod +x setup_hyb2.sh
./setup_hyb2.sh
```
`chmod +x` makes the shell script runnable, and then `./setup_hyb2.sh` runs it.

The setup file creates the file `~/.hyb2.sh` which contains a simple function wrapper that simplifies the hyb2 commands (instead of rewriting the long docker run command each time):

### Linux/WSL:
```bash
hyb2() {
    local tty=""; [ -t 1 ] && tty="-t"
    docker run --rm $tty --user "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "$PWD:/data" -w /data ghcr.io/tomharcus/hyb2:latest hyb2 "$@"
}
```

### macOS:
```bash
hyb2() {
    local tty=""; [ -t 1 ] && tty="-t"
    docker run --rm $tty -v "$PWD:/data" -w /data ghcr.io/tomharcus/hyb2:latest hyb2 "$@"
}
```

It then writes one line to your shells startup configuration file (`~/.bashrc` or `~/.zshrc`):
```bash
[ -f ~/.hyb2.sh ] && source ~/.hyb2.sh
```
This makes sure that the `~/.hyb2.sh` file exists in your home directory and if it does, it executes it.

Then once you run the setup script close the terminal and reopen it again.

Afterwards you are ready to start.
Type `hyb2` into your terminal to test that installation worked and to see information on how to use the pipeline.

For help on a specific functionality type:
```bash
hyb2 coverage --help
hyb2 fold --help
hyb2 compare --help
```

You can now go to the next section on how to run the pipeline.

## Eddie Installation

Eddie uses **Apptainer** (not Docker) and is a shared batch cluster. The login node is limited, so compute intensive tasks should
be run on compute nodes via `qsub` (job batching) or `qlogin` (interactive).

To get the Apptainer image on Eddie run:
```bash
qlogin -l h_vmem=16G                    # needs a compute node, NOT A LOGIN NODE
cd /exports/eddie/scratch/$USER         # change to the scratch directory as the home quota is too small
module load apptainer
apptainer pull docker://ghcr.io/tomharcus/hyb2:latest
```
The pull converts to a `.sif` file.

For interactive use you can set up the wrapper. Fetch `setup_hyb2.sh` and run it with the `.sif` path:
```bash
curl -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/setup_hyb2.sh
chmod +x setup_hyb2.sh                      
./setup_hyb2.sh /exports/eddie/scratch/$USER/hyb2_latest.sif
source ~/.bashrc
```
>Important caveat: the `hyb2` wrapper works only in an interactive session (`qlogin`), but not inside a `qsub` batch job.
>Batch jobs run a non-interactive shell that doesn't source `~/.bashrc`, so the function is not available. For batch jobs
>use the full `qsub` batch template in [Running on Eddie](#running-the-pipeline-on-eddie-with-batch-jobs).

### Eddie Specific Quirks
- Work on the scratch directory, not home. The home directories quota is too small to run the pipeline.
- Request cores: the default is 1, and this makes bowtie2 very slow.
- Keep Apptainer's cache/tmp off of scratch by running:
  ```bash
  export APPTAINER_CACHEDIR=/exports/eddie/scratch/$USER/apptainer_cache
  export APPTAINER_TMPDIR=/exports/eddie/scratch/$USER/apptainer_tmp
  ```
- Bind scratch + `$TMPDIR` as Apptainer auto mounts to `$HOME` but not those (available in the template in [Running on Eddie](#running-the-pipeline-on-eddie-with-batch-jobs))

Afterwards you are ready to start.
To make sure installation worked and see information on how to use the pipeline start a new interactive session and load Apptainer:
```bash
qlogin -l h_vmem=16G    # request enough memory
module load apptainer
```

Then type `hyb2` into the terminal.

Within the interactive session, for help on a specific functionality type:
```bash
hyb2 coverage --help
hyb2 fold --help
hyb2 compare --help
```

You can now go to the next section on how to run the pipeline.


## Running hyb2

hyb2 is one command with four subcommands. You run it from your data directory. You can reference files that are in the directory by their plain name 
or if somewhere else on your machine, you can provide the full path. The following information is for running the pipeline locally or on the cluster interactively.

### The Different Commands

| Purpose | Command | 
|---|---|
| main orchestrator | `hyb2` |
| folding | `hyb2 fold` | 
| coverage / CDM | `hyb2 coverage` | 
| dataset comparison | `hyb2 compare` | 


Example full pipeline run:
```bash
cd data_dir
hyb2 -i reads.sam -d ref.fasta -o myrun -a MyRNA -x 3900 -l 300
```
- `-i` input, `-d` reference, `-o` output name, `-a` gene of interest, `-x`/`-l` fold window
- the example produces the file `myrun.hyb`, the contact density map, the viewpoint graph, and the folded structure files

The pipeline has different input types: `.sam` files skip the mapping stage, `.fastq`, `.fastq.gz`, `.fasta`, `.fasta.gz` get mapped first, and `.hyb` skips to plotting
and folding. The pipeline autodetects the input.

### The Four Analysis Modes
```bash
# short-range intramolecular
hyb2 -i reads.sam -d ref.fasta -o test -a MyRNA -x 1001 -l 500
# long-range intramolecular
hyb2 -i reads.sam -d ref.fasta -o test -a MyRNA -x 1001 -y 5001 -l 500
# intermolecular (two genes)
hyb2 -i reads.sam -d ref.fasta -o test -a RNA_A -b RNA_B -x 7501 -y 501 -l 500
# homodimer
hyb2 -i reads.sam -d ref.fasta -o test -a MyRNA -b MyRNA -x 3501 -y 3501 -l 200
```

The `.hyb` file and database files can be re-used to save time and test different inputs for example:
```bash
hyb2 -i reads.sam -d ref.fasta -o test      # slow: call chimeras once
hyb2 -i test.hyb -a RNA_A -x 100 -l 300     # fast: fold any RNA, no re-mapping
hyb2 -i test.hyb -a RNA_B -x 500 -l 300
```

Here is a list of the key flags:
| Flag | purpose |
|---|---|
| -b | second gene | 
| -r | folding backend: `cplfold` (default) / `vienna` / `unafold` | 
| -v | BLAST e-value threshold (default=0.1) | 
| -m | max overlap (default=4) |
| -h | max hits per read (default=10) |
| -q | heatmap upper quantile (default=0.95) |
| -V | verbose (shows per stage detail + subprocess errors) |
| --config | read args from a YAML file |
| --reproducible | byte-deterministic run | 


**For complete information about each flag type: `hyb2 --help`**

## Running The Pipeline on Eddie with Batch Jobs

The hyb2 wrapper works in an interactive `qlogin` session, but **not** inside a batch job. A batch job runs on a non-interactive
shell that doesn't load your `~/.bashrc`, so the `hyb2` function is not available. For batch jobs you write a small job script
that calls the container directly with the full `apptainer run` command, then submit it with `qsub`.

Create a file, e.g. `run_hyb2.sh`, in your scratch run directory and paste in this:
```bash
#!/bin/bash
#$ -cwd                      # run in the directory you submit from
#$ -N hyb2                   # job name (shows up in qstat)
#$ -pe sharedmem 32          # number of cores (bowtie2 uses these)
#$ -l h_vmem=8G              # memory PER CORE (so 32 x 8G = 256G total here)
#$ -l h_rt=12:00:00          # max runtime (hh:mm:ss)

. /etc/profile.d/modules.sh  # makes `module` available in the batch shell
module load apptainer/1.4.4

export APPTAINER_TMPDIR="$TMPDIR"

apptainer run --bind "$TMPDIR" --bind "$TMPDIR":/tmp --bind /exports/eddie/scratch/$USER \
    /exports/eddie/scratch/$USER/<your_dir>/hyb2_latest.sif \
    hyb2 -i reads.sam -d ref.fasta -o myrun -a MyRNA -x 3900 -l 300
```

The last line is the actual pipeline command. You can swap it for **any command from the sections above** (the full pipeline, `hyb2 fold`,
`hyb2 compare`, ...). Everything before it is the cluster wrapping: request resources, load Apptainer, and bind the directories as Apptainer
doesn't mount automatically (scratch and the node-local `$TMPDIR`, where the big sort spills).

Then you can submit the job and check its status:
```bash
qsub run_hyb2.sh      # queues the job
qstat                 # shows your queued/running jobs
```

### A few things to get right:
- `h_vmem` is per core, so total memory is `h_vmem * cores`. The `32 * 8G = 256G` above is more than enough. Do not set `h_vmem` too low or mapping runs out of memory.
- You can reuse the reference index. Keep `ref.fasta` and its generated `.bt2` / `.tab` files together in a persistent directory on scratch and you can point the `-d` flag at it,
skipping the database creation step.
- The sort buffer is 25% of the node's RAM, so keep `h_vmem * cores` above that. A much tighter allocation could run out of memory on a large sort.






