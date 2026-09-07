"""Shared test helpers. Path locators + a legacy-script runner.

These are methodology-agnostic (usable whether the suite diffs against the
captured golden fixtures or re-runs the legacy scripts). 
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
BIN = REPO_ROOT / "bin"
FIXTURES = REPO_ROOT / "fixtures" / "sam_composition_run"
TIER2_FIXTURES = REPO_ROOT / "fixtures" / "tier2_run"
FOLDING_FIXTURES = REPO_ROOT / "fixtures" / "folding_run"
COVERAGE_FIXTURES = REPO_ROOT / "fixtures" / "coverage_run"
VIEWPOINT_FIXTURES = REPO_ROOT / "fixtures" / "viewpoint_run"
UNAFOLD_FIXTURES = REPO_ROOT / "fixtures" / "unafold_run"


@pytest.fixture
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture
def bin_dir() -> Path:
    return BIN


@pytest.fixture
def fixtures_dir() -> Path:
    """Golden outputs from scripts/generate_sam_composition_baseline.sh.
    Skips the test if they haven't been generated on this machine."""
    if not FIXTURES.exists():
        pytest.skip(f"fixtures not generated; run scripts/generate_sam_composition_baseline.sh")
    return FIXTURES


@pytest.fixture
def tier2_fixtures_dir() -> Path:
    """Golden outputs from scripts/generate_tier2_baseline.sh (the folding-
    adjacent Tier 2 stages). Skips the test if they haven't been generated."""
    if not TIER2_FIXTURES.exists():
        pytest.skip("tier2 fixtures not generated; run scripts/generate_tier2_baseline.sh")
    return TIER2_FIXTURES


@pytest.fixture
def folding_fixtures_dir() -> Path:
    """Real RNA-folding oracle from scripts/generate_folding_baseline.sh
    (ViennaRNA .ct/.vienna on the test data, plus legacy parser goldens).
    Requires ViennaRNA installed to (re)generate; skips otherwise."""
    if not FOLDING_FIXTURES.exists():
        pytest.skip("folding fixtures not generated; run scripts/generate_folding_baseline.sh (needs ViennaRNA)")
    return FOLDING_FIXTURES


@pytest.fixture
def vienna_bin() -> str:
    """Directory containing ViennaRNA's RNAcofold/b2ct, for pipeline tests that
    actually fold. Prefers PATH, falls back to the hyb2 conda env; skips if
    neither has it."""
    import os
    exe = shutil.which("RNAcofold")
    if exe is not None:
        return str(Path(exe).parent)
    candidate = Path.home() / Path(os.environ["CONDA_PREFIX"]) / "bin" / "RNAcofold"
    if candidate.exists():
        return str(candidate.parent)
    pytest.skip("ViennaRNA (RNAcofold) not found on PATH or in the hyb2 env")


@pytest.fixture
def coverage_fixtures_dir() -> Path:
    """Legacy plot_hybrids_3.awk parity oracle from
    scripts/generate_coverage_baseline.sh. Skips if not generated."""
    if not COVERAGE_FIXTURES.exists():
        pytest.skip("coverage fixtures not generated; run scripts/generate_coverage_baseline.sh")
    return COVERAGE_FIXTURES


@pytest.fixture
def viewpoint_fixtures_dir() -> Path:
    """Legacy hyb2blast.awk / blast2gplot.pl parity oracles from
    scripts/generate_viewpoint_baseline.sh. Skips if not generated."""
    if not VIEWPOINT_FIXTURES.exists():
        pytest.skip("viewpoint fixtures not generated; run scripts/generate_viewpoint_baseline.sh")
    return VIEWPOINT_FIXTURES


@pytest.fixture
def cplfold_env() -> str:
    """Skip unless CPLfold is importable from config.CPLFOLD_DIR"""
    import os
    import sys as _sys
    from hyb2.tools import config

    cpldir = config.CPLFOLD_DIR
    if not os.path.isdir(cpldir):
        pytest.skip(f"CPLfold not found at {cpldir}")
    if cpldir not in _sys.path:
        _sys.path.insert(0, cpldir)
    try:
        import CPLfold  # noqa: F401
    except Exception as e:
        pytest.skip(f"CPLfold not importable: {e}")

    return cpldir


@pytest.fixture
def unafold_fixtures_dir() -> Path:
    """UNAFold (OligoArrayAux) parity oracle from
    scripts/generate_unafold_baseline.sh. Skips if not generated."""
    if not UNAFOLD_FIXTURES.exists():
        pytest.skip("unafold fixtures not generated; run scripts/generate_unafold_baseline.sh")
    return UNAFOLD_FIXTURES


@pytest.fixture
def unafold_env(monkeypatch) -> str:
    """Skip unless hybrid-ss-min AND hybrid-min (oligoarrayaux) are runnable.
    Prepends the hyb2 conda env's bin to PATH so the port's subprocesses find
    them when the env isn't activated (UNAFOLDDAT is set by the port itself from
    config.UNAFOLD_DIR). Returns the bin dir."""
    import os

    bindir = None
    if shutil.which("hybrid-ss-min") and shutil.which("hybrid-min"):
        bindir = Path(shutil.which("hybrid-ss-min")).parent
    else:
        cand = Path.home() / Path(os.environ["CONDA_PREFIX"]) / "bin"
        if (cand / "hybrid-ss-min").exists() and (cand / "hybrid-min").exists():
            monkeypatch.setenv("PATH", str(cand) + os.pathsep + os.environ.get("PATH", ""))
            bindir = cand
    if bindir is None:
        pytest.skip("hybrid-ss-min/hybrid-min (oligoarrayaux) not found")
    return str(bindir)


@pytest.fixture
def varna_jar() -> str:
    """Path to the VARNA jar, for tests that actually render (plot_VARNA).
    Skips if java or the jar is missing."""
    if shutil.which("java") is None:
        pytest.skip("java not found")
    jar = REPO_ROOT / "VARNA" / "build" / "jar" / "VARNAcmd.jar"
    if not jar.exists():
        pytest.skip("VARNA jar not found")
    return str(jar)


def run_legacy(cmd: list[str], stdin_text: str | None = None) -> str:
    """Run a legacy bin/ script (perl/awk/bash) and return its stdout. Kept for
    ad-hoc use; the Tier 1 suite diffs against captured golden fixtures rather
    than re-running legacy at test time."""
    result = subprocess.run(
        cmd, input=stdin_text, capture_output=True, text=True, check=True, cwd=REPO_ROOT
    )
    return result.stdout


@pytest.fixture
def strip_hyb_header():
    """Callable that drops the '#'-comment header from .hyb text (see data_lines)."""
    return data_lines


def data_lines(text: str) -> str:
    """Drop '#'-comment and blank lines. get_mtop_hybrids' .hyb output carries a
    comment header that embeds the invocation path/args (differs between the
    Perl and the Python port), so the parity contract is on data rows only."""
    return "".join(
        line + "\n"
        for line in text.splitlines()
        if line.strip() and not line.startswith("#")
    )
