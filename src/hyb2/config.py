"""

Shared paths/ external tool locations

"""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def varna_jar() -> str:
    """
    Path to VARNAcmd.jar. Port of VARNAcmd.sh
    """
    env = os.environ.get("HYB2_VARNA_JAR")

    if env:
        return env
    
    return (REPO_ROOT / "bin" / "VARNA.dir").read_text().splitlines()[0]