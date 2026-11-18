"""

Shared paths/ external tool locations

"""

import os, sys, subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]

def rscript(name: str) -> str:  
    """
    Path to R scripts in bin/
    """
    return str(REPO_ROOT / "rscripts" / name)


def gnu_sort() -> str:
    """
    MacOS uses BSD sort by default so identical data would produce different results on
    MacOS and Linux. Prefer a verified GNU sort for byte-parity
    """
    for cand in ("sort", "gsort"):
        try:
            if "GNU coreutils" in subprocess.run(
                [cand, "--version"], capture_output=True, text=True).stdout:
                return cand
        except FileNotFoundError:
            continue
    raise RuntimeError("GNU sort not found: install coreutils (macOS: 'brew install coreutils' or conda 'coreutils')")

def str2bool(v):
    return str(v).strip().lower() in ("true", "1")

VARNA_JAR = os.environ.get("HYB2_VARNA_JAR") or str(REPO_ROOT / "VARNA" / "build" / "jar" / "VARNAcmd.jar") 

CPLFOLD_DIR = os.environ.get("HYB2_CPLFOLD_DIR") or str(REPO_ROOT / "CPLfold") # path to cplfold
# beam 200 / no pseudoknots: on full-length yeast U2 (1,175 nt) with chimera data, beam 100 lost
# the 3' domain while beam >= 200 recovered the 2011 structure; without pseudoknots hyb2 skips
# CPLfold's slow suboptimal-structure step (minutes -> seconds).
CPL_DEFAULTS = {"alpha": 0.5, "beta": 0.0, "normalize": "log", "beam_size": 200, "energy_delta": 5.0, "max_phase1": 10, "max_phase2": 5, "allow_pseudoknot": False, "energy_model": "DP09", "support": "cofold"}

UNAFOLD_DIR = Path(sys.prefix)/"share/oligoarrayaux" # path to unafold

# mfold 3.6 (bioconda "mfold"): replaces hybrid-ss-min for constrained folding.
# HYB2_MFOLD can point at the mfold script if it is not on PATH (e.g. a separate conda env).
import shutil as _shutil
MFOLD_BIN = os.environ.get("HYB2_MFOLD") or _shutil.which("mfold") or "mfold"
MFOLD_DAT_DIR = os.environ.get("HYB2_MFOLD_DAT") or str(Path(MFOLD_BIN).resolve().parents[1] / "share" / "mfold")