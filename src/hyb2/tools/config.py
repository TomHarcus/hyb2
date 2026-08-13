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
    return str(REPO_ROOT / "bin" / name)


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



VARNA_JAR = os.environ.get("HYB2_VARNA_JAR") or str(REPO_ROOT / "VARNA" / "build" / "jar") # Path to varna_jar

CPLFOLD_DIR = os.environ.get("HYB2_CPLFOLD_DIR") or str(REPO_ROOT / "CPLfold") # path to cplfold
CPL_DEFAULTS = {"alpha": 0.5, "beta": 0.0, "normalize": "log", "beam_size": 100, "energy_delta": 5.0, "max_phase1": 10, "max_phase2": 5, "energy_model": "DP09"}

UNAFOLD_DIR = Path(sys.prefix)/"share/oligoarrayaux" # path to unafold