""" Port of bowtie2_fasta2sam, bowtie2_fastq.gz2sam, and bowtie2_fastq2sam all condensed into one script

"""

from hyb2.mapping.make_hyb_db_2 import make_hyb_db_2
from hyb2.mapping.make_comp_fasta import make_comp_fasta
from hyb2.mapping.solexa2fasta import solexa_to_fasta_lines
from hyb2.common.fasta2tab import fasta_to_tab_lines
from hyb2.tools.ui import spinner
from hyb2.tools.logsetup import is_quiet

import subprocess, gzip, os, contextlib

from pathlib import Path

def bowtie2_map(in_file, db, out, reproducible=False):

    suffix = in_file.split(".")

    # check if in_file is a .fasta file
    if suffix[-1] == "fasta":

        make_hyb_db_2(db)

        _bowtie2(db, out, in_file, reproducible=reproducible)
        
        print("Mapping concluded")

    # support for .fasta.gz files
    elif suffix[-2] == "fasta" and suffix[-1] == "gz":

        make_hyb_db_2(db)
        
        _bowtie2(db, out, in_file, reproducible=reproducible)
                
        print("Mapping concluded")

    # check if in_file is a .fastq.gz file
    # uses generators + streaming so that when reading large input file ram limit doesnt shoot up
    elif suffix[-2] == "fastq" and suffix[-1] == "gz":


        make_hyb_db_2(db)

        comp_path = f"{out}_comp.fasta"

        with gzip.open(in_file, "rt") as fin:
            fasta = solexa_to_fasta_lines(fin)
            tab = fasta_to_tab_lines(fasta)
            comp = make_comp_fasta(tab)
            Path(comp_path).write_text(comp)

        _bowtie2(db, out, comp_path, reproducible=reproducible)

        print("Mapping concluded")

    # check if in_file is a .fastq file
    # uses generators + streaming so that when reading large input file ram limit doesnt shoot up
    elif suffix[-1] == "fastq":
    
        make_hyb_db_2(db)

        comp_path = f"{out}_comp.fasta"

        with open(in_file) as fin:
            fasta = solexa_to_fasta_lines(fin)
            tab = fasta_to_tab_lines(fasta)
            comp = make_comp_fasta(tab)
            Path(comp_path).write_text(comp)

        _bowtie2(db, out, comp_path, reproducible=reproducible)

        print("Mapping concluded")

    else:
        raise ValueError(f"unsupported input format: {in_file}")

# run bowtie2
def _bowtie2(db, out, reads, reproducible):
    quiet = is_quiet()
    threads = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)

    cmd = ["bowtie2", "-D", "20", "-R", "3", "-N", "0", "-L", "16", "-k", "20", "--local",
                    "-i", "S,1,0.50", "--score-min", "L,18,0", "--ma", "1", "--np", "0", "--mp", "2,2",
                    "--rdg", "5,1", "--rfg", "5,1", "-p", str(threads), "-x", db.replace(".fasta", "", 1),
                    "-f", reads]

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