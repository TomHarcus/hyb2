#!/usr/bin/env python3

"""
HYB2 one-shot installer

Run once per machine: python3 install.py

Stdlib is the only necessity: run this BEFORE the hyb2 conda env exists, so it
cannot import numpy/yaml/etc. It orchestrates conda/git/make and then hands off
to the hyb2 launcher for actual runs

Platforms: Linux, Intel MacOS, Apple Silicon MacOS (osx-64 via Rosetta).
Windows is not supported natively, instead run within WSL2

"""

import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
ENV_NAME = "hyb2"
CPLFOLD_URL = "https://github.com/Vicky-0256/CPLfold.git"
CPLFOLD_COMMIT = "6f49167b25cf6312c45410a4a3b7f22ae1ca54b9"
CPLFOLD_DIR = REPO / "CPLfold"
HOTKNOTS_DIR = CPLFOLD_DIR / "Utils" / "HotKnots_v2.0"

# output helpers
def info(msg):  print(f"\033[1;34m==>\033[0m {msg}")
def ok(msg):    print(f"\033[1;32m  ok\033[0m {msg}")
def die(msg):   print(f"\033[1;31mERROR:\033[0m {msg}", file=sys.stderr); sys.exit(1)

def run(cmd, *, cwd=None, check=True, capture=False, env=None):
    # this subprocess wrapper fails loud with the command that broke

    print(f'\t${" ".join(map(str, cmd))}')
    try:
        r = subprocess.run(
            [str(c) for c in cmd], cwd=cwd, env=env,
            text=True, capture_output=capture
        )

    except FileNotFoundError:
        die(f"command not found: {cmd[0]}")

    if check and r.returncode != 0:
        if capture and r.stderr:
            print(r.stderr, file=sys.stderr)
        die(f'command failed ({r.returncode}): {" ".join(map(str, cmd))}')

    return r


# platform
def detect_platform():
    osname = platform.system()
    if osname == "Windows":
        die("Native Windows is not supported. Install WSL2, then re-run this installer inside WSL2")

    mac_arm = False
    if osname == "Darwin":
        # check if intel or arm mac
        r = run(["sysctl", "-n", "hw.optional.arm64"], check=False, capture=True)
        mac_arm = r.stdout.strip() == "1"

    return osname, mac_arm

def mac_prep(mac_arm):
    # compiler for hotknots
    if run(["xcode-select", "-p"], check=False, capture=True).returncode != 0:
        info("Installing Xcode command-line tools (GUI diaglog)...")
        run(["xcode-select", "--install"], check=False)
        die("Xcode CLT install was triggered. Complete the GUI dialog, then re-run this installer")
    ok("Xcode command-line tools present")

    # check if apple silicon: need rosetta to run osx-64 bioconda libraries
    if mac_arm:
        info("Apple silicon detected: ensuring Rosetta 2 is installed")
        run(["softwareupdate", "--install-rosetta", "--agree-to-license"], check=False)

    # build the env as osx-64 on both Mac types
    os.environ["CONDA_SUBDIR"] = "osx-64"
    ok("CONDA_SUBDIR=osx-64 (uniform MAC env)")

# conda
def find_conda():
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if not conda:
        die("conda not found on PATH. Install Miniconda/Miniforge first")

    base = run([conda, "info", "--base"], capture=True).stdout.strip()
    ok(f"conda base: {base}")
    return conda

def env_exists(conda):
    out = run([conda, "env", "list"], capture=True).stdout
    return any(line.split() and line.split()[0] == ENV_NAME for line in out.splitlines())

def create_env(conda):
    yml = REPO / "environment.yml"
    if not yml.is_file():
        die(f"environment.yml not found at {yml}")

    if env_exists(conda):
        info(f"env {ENV_NAME} exists, updating (prune)")
        run([conda, "env", "update", "-n", ENV_NAME, "-f", yml, "--prune"], cwd=REPO)
    else:
        info(f"Creating conda env {ENV_NAME} (this pulls R/DESeq2 so might take a few minutes)")
        run([conda, "env", "create", "-f", yml], cwd=REPO)
    ok(f"env {ENV_NAME} is ready")

def env_prefix(conda):
    return run([conda, "run", "-n", ENV_NAME, "python", "-c", "import sys;print(sys.prefix)"], capture=True).stdout.strip()

# cplfold + hotknots
def clone_cplfold():
    if not shutil.which("git"):
        die("git not found on PATH")

    if (CPLFOLD_DIR / ".git").is_dir():
        info(f"CPLfold present: fetching and pinning {CPLFOLD_COMMIT[:7]}")
        run(["git", "fetch", "--all"], cwd=CPLFOLD_DIR)

    else:
        info(f"Cloning CPLfold into {CPLFOLD_DIR}")
        run(["git", "clone", CPLFOLD_URL, str(CPLFOLD_DIR)])

    run(["git", "checkout", CPLFOLD_COMMIT], cwd=CPLFOLD_DIR)
    ok(f"CPLfold at {CPLFOLD_COMMIT[:7]}")

def build_hotknots():
    if not shutil.which("make"):
        die("make not found (Linux: install build-essential; MacOS: Xcode CLT)")

    if not HOTKNOTS_DIR.is_dir():
        die(f"HotKnots dir missing: {HOTKNOTS_DIR}")

    info("Building HotKnots (native arch)")

    # purge stale objects/archives so a cross-arch rebuild is clean
    for f in (*HOTKNOTS_DIR.rglob("*.o"), *HOTKNOTS_DIR.rglob("*.a")):
        f.unlink()

    # the graphics/ objects are NOT built by the top Makefile's `default` target,
    run(["make", "graphics.o", "PlotRna.o"], cwd=HOTKNOTS_DIR / "graphics")
    run(["make"], cwd=HOTKNOTS_DIR)
    ok("HotKnots built")

# activate.d
def write_activate_d(prefix):
    d = Path(prefix) / "etc" / "conda" / "activate.d"
    d.mkdir(parents=True, exist_ok=True)

    script = d / "hyb2.sh"
    script.write_text(
        "#!/bin/sh\n"
        "# Written by hyb2 install.py - runtime env for the pipeline\n"
        'export UNAFOLDDAT="$CONDA_PREFIX/share/oligoarrayaux"\n'
        f"export HYB2_VARNA_JAR='{REPO}/VARNA/build/jar/VARNAcmd.jar'\n"
        f"export HYB2_CPLFOLD_DIR='{CPLFOLD_DIR}'\n"
        # default TMPDIR to real disk if unset or tmpfs then collapse_blast spills
        # to TMPDIR and OOMs on a tmpfs (RAM-backed) dir. findmnt is Linux-only;
        # 2>/dev/null makes this a no-op on macOS (where /tmp isn't tmpfs anyway).
        'if [ -z "$TMPDIR" ] || '
        '[ "$(findmnt -no FSTYPE --target "${TMPDIR:-/tmp}" 2>/dev/null)" = "tmpfs" ]; then\n'
        '    export TMPDIR="$HOME/scratch_tmp"\n'
        '    mkdir -p "$TMPDIR"\n'
        "fi\n"
    )


    ok(f"wrote {script}")

# preflight
def preflight(conda):
    info("Preflight checks")

    def cr(*args, **kw):
        return run([conda, "run", "-n", ENV_NAME, *args], check=False, capture=True, **kw)

    # GNU sort: parity critical for collapse_blast (MacOS default is BSD sort)
    if "GNU coreutils" in cr("sort", "--version").stdout:
        ok("GNU sort")
    else:
        die("sort in the env is NOT GNU coreutils: collapse_blast parity would break. Check the coreutils package installed")

    # external tool binaries
    for tool in ("bowtie2", "RNAcofold", "hybrid-ss-min", "hybrid-min", "java", "Rscript"):
        if cr("which", tool).returncode == 0:
            ok(tool)
        else:
            die(f"required tool missing from env: {tool}")

    # DESeq2 loads
    if cr("Rscript", "-e", "suppressMessages(library(DESeq2))").returncode == 0:
        ok("R DESeq2 loads")
    else:
        die("DESeq2 failed to load in the env's R")

    # HotKnots arch smoke
    ce = HOTKNOTS_DIR / "bin" / "computeEnergy"
    if not ce.is_file():
        die(f"HotKnots binary not built: {ce}")

    try:
        subprocess.run([str(ce)], capture_output=True, timeout=10)
        ok("HotKnots computeEnergy executes (arch OK)")
    except OSError as e:
        die(f"HotKnots computeEnergy won't execute ({e}). Likely an architecture mismatch. Rebuild: make -C {HOTKNOTS_DIR}")
    except subprocess.TimeoutExpired:
        ok("HotKnots computeEnergy executes (arch OK)")

def ensure_tmpdir():
    tmp = os.environ.get("TMPDIR") or str(Path.home() / "scratch_tmp")
    os.makedirs(tmp, exist_ok=True)
    os.environ["TMPDIR"] = tmp
    ok(f"TMPDIR={tmp}")

def main():
    info(f"HYB2 installer (repo: {REPO})")
    ensure_tmpdir()
    osname, mac_arm = detect_platform()
    ok(f"platform {osname}{' (Apple Silicon)' if mac_arm else ''}")

    conda = find_conda()
    if osname == "Darwin":
        mac_prep(mac_arm)

    create_env(conda)
    clone_cplfold()
    build_hotknots()
    write_activate_d(env_prefix(conda))
    preflight(conda)

    print()
    info("Install complete")
    print("\tActivate:  conda activate hyb2")
    print("\tRun:       hyb2-py --config run.yml     (edit run.yml first)")

if __name__ == "__main__":
    sys.exit(main())