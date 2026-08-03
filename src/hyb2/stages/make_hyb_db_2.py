""" Port of make_hyb_db_2

"""

from pathlib import Path
import subprocess
from hyb2.stages.fasta2tab import fasta_to_tab

def make_hyb_db_2(in_fasta):

    if not Path(in_fasta).is_file():
        raise FileNotFoundError(f"error: {in_fasta} file not found")

    prefix = in_fasta.replace(".fasta", "", 1)
    subprocess.run(["bowtie2-build", in_fasta, prefix], check=True, stdout=subprocess.DEVNULL)

    tab = in_fasta.replace("fasta", "tab", 1)
    with open(in_fasta) as fin:
        Path(tab).write_text(fasta_to_tab(fin.read()))