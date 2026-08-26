""" Port of make_hyb_db_2

"""

from pathlib import Path
import subprocess
from hyb2.common.fasta2tab import fasta_to_tab
from hyb2.tools.ui import spinner
import logging
import contextlib
from hyb2.tools.logsetup import is_quiet

def make_hyb_db_2(in_fasta):

    quiet = is_quiet()

    if not Path(in_fasta).is_file():
        raise FileNotFoundError(f"error: {in_fasta} file not found")

    prefix = in_fasta.replace(".fasta", "", 1)
    tab = in_fasta.replace("fasta", "tab", 1)

    if not (Path(f"{prefix}.rev.2.bt2").is_file() or Path(f"{prefix}.rev.2.bt2l").is_file()):

        with spinner("making database ") if quiet else contextlib.nullcontext():
            if not quiet:
                print("making database:")
            subprocess.run(["bowtie2-build", in_fasta, prefix], check=True, 
                        stdout=subprocess.DEVNULL if quiet else None,
                        stderr=subprocess.DEVNULL if quiet else None
                        )
    if not Path(tab).is_file():
        with open(in_fasta) as fin:
            Path(tab).write_text(fasta_to_tab(fin.read()))

    else:
        print("Database already generate: skipping to bowtie2 mapping")