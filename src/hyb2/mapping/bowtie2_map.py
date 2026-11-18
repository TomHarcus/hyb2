""" Port of bowtie2_fasta2sam, bowtie2_fastq.gz2sam, and bowtie2_fastq2sam all condensed into one script

Changes from the legacy scripts (gk-dev):
  * reads are 3' adapter-trimmed with cutadapt before collapsing (MAPPING_DEFAULTS["adapter"];
    adapter=None skips trimming). Already-collapsed FASTA is not trimmed.
  * bowtie2 seed length and minimum score default to -L 10 / --score-min L,15,0 (legacy:
    -L 16 / L,18,0), and both can be set by the caller; extra bowtie2 arguments can be appended.
    On the Brr2 CLASH data (Kudla et al. 2011) trimming + these settings recovered 97% of the
    paper's U2 chimeras vs 58% with the legacy settings, mostly short (<= 20 nt) chimera arms.
    Provisional defaults: to be benchmarked on other datasets.

"""

from hyb2.mapping.make_hyb_db_2 import make_hyb_db_2
from hyb2.mapping.make_comp_fasta import make_comp_fasta
from hyb2.mapping.solexa2fasta import solexa_to_fasta_lines
from hyb2.common.fasta2tab import fasta_to_tab_lines
from hyb2.tools.ui import spinner
from hyb2.tools.logsetup import is_quiet

import subprocess, gzip, os, contextlib, re, shlex, shutil, logging

log = logging.getLogger(__name__)

MAPPING_DEFAULTS = {
    "adapter": "TGGAATTCTCGGGTGCCAAGG",   # Illumina small-RNA 3' adapter (CLASH / CRAC libraries)
    "min_length": 18,                       # reads shorter than this after trimming are dropped
    "seed_length": 10,                      # bowtie2 -L (legacy 16)
    "score_min": "L,15,0",                  # bowtie2 --score-min (legacy L,18,0)
}

# read IDs as written by make_comp_fasta: <nr>_<count> or <nr>-<barcodes>_<count>
_COLLAPSED_ID = re.compile(r"\d+(-\d+)?_\d+")

def bowtie2_map(in_file, db, out, reproducible=False, adapter=MAPPING_DEFAULTS["adapter"],
                min_length=MAPPING_DEFAULTS["min_length"], seed_length=MAPPING_DEFAULTS["seed_length"],
                score_min=MAPPING_DEFAULTS["score_min"], bowtie2_args=None):

    bt2 = dict(seed_length=seed_length, score_min=score_min, bowtie2_args=bowtie2_args)

    suffix = in_file.split(".")
    is_fasta = suffix[-1] == "fasta" or suffix[-2:] == ["fasta", "gz"]

    if adapter:
        if is_fasta and _is_collapsed(in_file, gzip.open if suffix[-1] == "gz" else open):
            log.info("input reads are already collapsed: adapter trimming skipped")
        else:
            in_file = _trim_adapters(in_file, out, adapter, min_length)
            suffix = in_file.split(".")

    # check if in_file is a .fasta or .fasta.gz file
    # downstream stages read the read count from the ID suffix (<nr>_<count>), so
    # raw reads are collapsed first, as for fastq; already-collapsed input is
    # mapped as-is so its counts are kept
    if suffix[-1] == "fasta" or suffix[-2:] == ["fasta", "gz"]:

        make_hyb_db_2(db)

        opener = gzip.open if suffix[-1] == "gz" else open

        if _is_collapsed(in_file, opener):
            reads = in_file
        else:
            reads = f"{out}_comp.fasta"
            with opener(in_file, "rt") as fin, open(reads, "w") as fout:
                fout.writelines(make_comp_fasta(fasta_to_tab_lines(fin)))

        _bowtie2(db, out, reads, reproducible=reproducible, **bt2)

        print("Mapping concluded")

    # check if in_file is a .fastq.gz file
    # uses generators + streaming so that when reading large input file ram limit doesnt shoot up
    elif suffix[-2] == "fastq" and suffix[-1] == "gz":


        make_hyb_db_2(db)

        comp_path = f"{out}_comp.fasta"

        with gzip.open(in_file, "rt") as fin, open(comp_path, "w") as fout:
            fasta = solexa_to_fasta_lines(fin)
            tab = fasta_to_tab_lines(fasta)
            fout.writelines(make_comp_fasta(tab))

        _bowtie2(db, out, comp_path, reproducible=reproducible, **bt2)

        print("Mapping concluded")

    # check if in_file is a .fastq file
    # uses generators + streaming so that when reading large input file ram limit doesnt shoot up
    elif suffix[-1] == "fastq":
    
        make_hyb_db_2(db)

        comp_path = f"{out}_comp.fasta"

        with open(in_file) as fin, open(comp_path, "w") as fout:
            fasta = solexa_to_fasta_lines(fin)
            tab = fasta_to_tab_lines(fasta)
            fout.writelines(make_comp_fasta(tab))

        _bowtie2(db, out, comp_path, reproducible=reproducible, **bt2)

        print("Mapping concluded")

    else:
        raise ValueError(f"unsupported input format: {in_file}")

# True if every read ID already has the collapsed <nr>_<count> form
def _is_collapsed(in_file, opener=open):
    with opener(in_file, "rt") as fin:
        return all(_COLLAPSED_ID.fullmatch((line[1:].split() or [""])[0])
                   for line in fin if line.startswith(">"))

def _trim_adapters(in_file, out, adapter, min_length):
    """3' adapter trimming with cutadapt; returns the trimmed reads file (same format, uncompressed)."""
    if shutil.which("cutadapt") is None:
        raise RuntimeError("cutadapt not found: install it (conda install -c bioconda cutadapt) "
                           "or switch trimming off with --adapter none")
    fmt = "fastq" if ".fastq" in in_file else "fasta"
    trimmed = f"{out}_trimmed.{fmt}"
    threads = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    with open(f"{out}.cutadapt.log", "w") as err:
        with spinner("adapter trimming ") if is_quiet() else contextlib.nullcontext():
            subprocess.run(["cutadapt", "-a", adapter, "-m", str(min_length), "-j", str(threads),
                            "-o", trimmed, in_file], stdout=err, stderr=subprocess.STDOUT, check=True)
    log.info(f"adapter-trimmed reads: {trimmed} (cutadapt report: {out}.cutadapt.log)")
    return trimmed

# run bowtie2
def _bowtie2(db, out, reads, reproducible, seed_length=MAPPING_DEFAULTS["seed_length"],
             score_min=MAPPING_DEFAULTS["score_min"], bowtie2_args=None):
    quiet = is_quiet()
    threads = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)

    cmd = ["bowtie2", "-D", "20", "-R", "3", "-N", "0", "-L", str(seed_length), "-k", "20", "--local",
                    "-i", "S,1,0.50", "--score-min", score_min, "--ma", "1", "--np", "0", "--mp", "2,2",
                    "--rdg", "5,1", "--rfg", "5,1", "-p", str(threads), "-x", db.replace(".fasta", "", 1),
                    "-f", reads]
    if bowtie2_args:
        cmd += shlex.split(bowtie2_args)
    log.info("bowtie2 command: " + " ".join(cmd))

    with open(f"{out}.sam", "w") as sam, open(f"{out}.blast.err", "w") as err:
        with spinner("bowtie2 mapping ") if quiet else contextlib.nullcontext():
            if not quiet:
                print("bowtie2 mapping:")

            if reproducible:
                cmd += ["--reorder"]

            subprocess.run(
                cmd,
                stdout=sam,
                stderr=err,
                check=True
            )