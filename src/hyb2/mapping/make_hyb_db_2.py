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

    with spinner("making database ") if quiet else contextlib.nullcontext():
        if not quiet:
            print("making database:")
        subprocess.run(["bowtie2-build", in_fasta, prefix], check=True, 
                       stdout=subprocess.DEVNULL if quiet else None,
                       )

    tab = in_fasta.replace("fasta", "tab", 1)
    with open(in_fasta) as fin:
        Path(tab).write_text(fasta_to_tab(fin.read()))