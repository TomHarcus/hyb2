""" Port of bowtie2_fasta2sam, bowtie2_fastq.gz2sam, and bowtie2_fastq2sam all condensed into one script

"""

from hyb2.stages.make_hyb_db_2 import make_hyb_db_2
from hyb2.stages.make_comp_fasta import make_comp_fasta

import subprocess

from pathlib import Path

def bowtie2_map(in_file, db, out):

    suffix = in_file.split(".")

    if suffix[-1] == "fasta":

        if not Path(db.replace("fasta", "tab", 1)).is_file():
            print("Making database...")
            make_hyb_db_2(db)

        print("Bowtie2 mapping...")

        with open(f"{out}.sam", "w") as sam, open(f"{out}.blast.err") as err:

            subprocess.run(
                ["bowtie2", "-D", "20", "-R", "3", "-N", "0", "-L", "16", "-k", "20", "--local",
                "-i", "S,1,0.50", "--score-min", "L,18,0", "--ma", "1", "--np", "0", "--mp", "2,2",
                "--rdg", "5,1", "--rfg", "5,1", "-p", "64", "-x", db.replace(".fasta", "", 1),
                "-f", in_file],
                stdout=sam,
                stderr=err,
                check=True
            )
        
        print("Mapping concluded")
        


