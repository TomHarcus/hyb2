"""Shared test helpers. Path locators + a legacy-script runner.

These are methodology-agnostic (usable whether the suite diffs against the
captured golden fixtures or re-runs the legacy scripts). The exact test style
for Tier 1 is being confirmed with Grzegorz before test_*.py files are added.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
BIN = REPO_ROOT / "bin"
FIXTURES = REPO_ROOT / "fixtures" / "sam_composition_run"
TIER2_FIXTURES = REPO_ROOT / "fixtures" / "tier2_run"
FOLDING_FIXTURES = REPO_ROOT / "fixtures" / "folding_run"


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
