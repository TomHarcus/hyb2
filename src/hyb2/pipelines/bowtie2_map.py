""" Port of bowtie2_fasta2sam, bowtie2_fastq.gz2sam, and bowtie2_fastq2sam all condensed into one script

"""

from hyb2.stages.make_hyb_db_2 import make_hyb_db_2
from hyb2.stages.make_comp_fasta import make_comp_fasta
from hyb2.stages.solexa2fasta import solexa_to_fasta
from hyb2.stages.fasta2tab import fasta_to_tab

import subprocess, gzip

from pathlib import Path

def bowtie2_map(in_file, db, out):

    suffix = in_file.split(".")

    # check if in_file is a .fasta file
    if suffix[-1] == "fasta":

        if not Path(db.replace("fasta", "tab", 1)).is_file():
            print("Making database...")
            make_hyb_db_2(db)

        print("Bowtie2 mapping...")

        _bowtie2(db, out, in_file)
        
        print("Mapping concluded")

    # check if in_file is a .fastq.gz file
    elif suffix[-2] == "fastq" and suffix[-1] == "gz":

        if not Path(db.replace("fasta", "tab", 1)).is_file():
            print("Making database...")
            make_hyb_db_2(db)

        unzipped_in_file = gzip.open(in_file, "rt").read()

        fasta = solexa_to_fasta(unzipped_in_file)
        tab = fasta_to_tab(fasta)
        comp = make_comp_fasta(tab.splitlines())

        comp_path = f"{out}_comp.fasta"
        Path(comp_path).write_text(comp)

        print("Bowtie2 mapping...")

        _bowtie2(db, out, comp_path)

        print("Mapping concluded")

    # check if in_file is a .fastq file
    elif suffix[-1] == "fastq":
        if not Path(db.replace("fasta", "tab", 1)).is_file():
            print("Making database...")
            make_hyb_db_2(db)

        fasta = solexa_to_fasta(open(in_file).read())
        tab = fasta_to_tab(fasta)
        comp = make_comp_fasta(tab.splitlines())

        comp_path = f"{out}_comp.fasta"
        Path(comp_path).write_text(comp)

        print("Bowtie2 mapping...")

        _bowtie2(db, out, comp_path)

        print("Mapping concluded")

    else:
        raise ValueError(f"unsupported input format: {in_file}")

# run bowtie2
def _bowtie2(db, out, reads):
    with open(f"{out}.sam", "w") as sam, open(f"{out}.blast.err", "w") as err:
            
        subprocess.run(
            ["bowtie2", "-D", "20", "-R", "3", "-N", "0", "-L", "16", "-k", "20", "--local",
            "-i", "S,1,0.50", "--score-min", "L,18,0", "--ma", "1", "--np", "0", "--mp", "2,2",
            "--rdg", "5,1", "--rfg", "5,1", "-p", "64", "-x", db.replace(".fasta", "", 1),
            "-f", reads],
            stdout=sam,
            stderr=err,
            check=True
        )