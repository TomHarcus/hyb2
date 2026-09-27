# HYB2 Quickstart Guide

This tool analyses RNA proximity ligation data. It takes in mapped reads (fastq/SAM) and tells you which RNA regions are physically interacting with each other, with
their sequence, coordinates, and predicted folded structure. Hyb2 outputs contact maps, viewpoint graphs, and the folded RNA structures.

This guide gets hyb2 setup and running as fast as possible using conda. For advanced use, see the full [README_python_port.md](README_python_port.md)

## 1. Install Miniforge

If you don't already have conda, install **Miniforge** (it ships the fast libmamba solver). Follow <https://github.com/conda-forge/miniforge#install>, or on Linux/WSL/macOS:
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"
```

Restart your shell, then check it worked:
```bash
conda --version
```

## 2. Get hyb2

Clone the repo and create the environment:
```bash
git clone -b python-migration https://github.com/TomHarcus/hyb2.git
cd hyb2

conda env create -f environment.yml     # env + all tools
conda activate hyb2

# CPLfold (pure python)
git clone https://github.com/Vicky-0256/CPLfold.git && git -C CPLfold checkout 24bab52
export HYB2_CPLFOLD_DIR="$PWD/CPLfold"
```

> **macOS (Apple Silicon):** `bioconductor-deseq2` has no arm64 build, so replace the create line with
> `CONDA_SUBDIR=osx-64 conda env create -f environment.yml` (builds the env as Intel, runs under Rosetta).

Check it worked:
```bash
hyb2
```

This should print the pipeline's information. In a new terminal, run `conda activate hyb2` first.

## 3. Verify it works with test data

```bash
mkdir hyb2_test && cd hyb2_test

curl -L -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/data/Zika_18S_formatted.fasta
curl -L -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/data/testData.sam

hyb2 -i testData.sam -d Zika_18S_formatted.fasta -o testrun -a Zika_virusRNA -x 3900 -l 300 -r cplfold
```

You should get a `testrun.hyb` file, a contact density map PDF, a viewpoint graph PDF, and a folded structure. If they are generated,
the pipeline is setup correctly.

## 4. Run it on your own data

```bash
hyb2 -i your_reads.sam -d your_reference.fasta -o my_run -a MyGeneName -x 3900 -l 300
```

- `-i`: your input reads (`.sam`, `.fastq`, `.fastq.gz`, `.fasta`, or `.fasta.gz` all work)
- `-d`: your reference genome/sequence file
- `-o`: a name for your output files
- `-a`: the gene/region you want to analyse
- `-x` / `-l`: the window (start position and length) to fold around

## 5. Next steps

This guide covers the basic pipeline. See the full [README_python_port.md](README_python_port.md) for running on Eddie, modifying the pipeline's code, tuning CPLfold, comparing datasets, 
or anything else.
