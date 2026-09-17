# HYB2 Quickstart Guide

This tool analyses RNA proximity ligation data. It takes in mapped reads (fastq/SAM) and tells you which RNA regions are physically interacting with each other, with
their sequence, coordinates, and predicted folded structure. Hyb2 outputs contact maps, viewpoint graphs, and the folded RNA structures.

This guide gets hyb2 setup and running as fast as possible using Docker. For advanced use, see the full [README_python_port.md](README_python_port.md)

## 1. Install Docker

- **Linux/WSL**: install Docker Engine (https://docs.docker.com/engine/install/)
- **macOS**: install Docker Desktop (https://www.docker.com/products/docker-desktop/)

Check if it worked:
```bash
docker run hello-world
```

## 2. Get hyb2

Pull the container:
```bash
docker pull ghcr.io/tomharcus/hyb2:latest
```

Download and run the setup script:
```bash
curl -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/setup_hyb2.sh
chmod +x setup_hyb2.sh
./setup_hyb2.sh
```

Close down the terminal and then reopen it. Check if it worked:
```bash
hyb2
```

This should print the pipelines information.

## 3. Verify it works with test data

```bash
mkdir hyb2_test && cd hyb2_test

curl -L -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/data/Zika_18S_formatted.fasta
curl -L -O https://raw.githubusercontent.com/TomHarcus/hyb2/python-migration/data/testData.sam

hyb2 -i testData.sam -d Zika_18S_formatted.fasta -o testrun -a Zika_virusRNA -x 3900 -l 300 -r cplfold
```

You should get a `testrun.hyb` file, a contact density map PDF, a viewpoint graph PDF, and a folded structure. If they are generated,
the pipeline is setup correctly.

## 4. Run it on you own data

```bash
hyb2 -i your_reads.sam -d your_reference.fasta -o my_run -a MyGeneName -x 3900 -l 300
```

- `-i`: your input reads (`.sam`, `.fastq`, `.fastq.gz`, `.fasta`, or `.fasta.gz` all work)
- `-d`: your reference genome/sequence file
- `-o`: a name for your output files
- `-a`: the gene/region you want to analyze
- `-x` / `-l`: the window (start position and length) to fold around

## 5. Next steps

This guide covers the basic pipeline. See the full [README_python_port.md](README_python_port.md) for running on Eddie, modifying the pipeline's code, tuning CPLfold, comparing datasets, 
or anything else.

