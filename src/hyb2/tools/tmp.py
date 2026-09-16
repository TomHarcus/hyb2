import os, tempfile, logging
log = logging.getLogger(__name__)

def _is_tmpfs(path: str) -> bool:
    try:
        with open("/proc/mounts") as f:
            mounts = [ln.split() for ln in f]
    except OSError:
        return False
    
    rp = os.path.realpath(path)
    best_mp, best_fs = "", ""
    for parts in mounts:
        if len(parts) < 3:
            continue
        mp, fstype = parts[1], parts[2]

        if (rp == mp or rp.startswith(mp.rstrip("/") + "/")) and len(mp) > len(best_mp):
            best_mp, best_fs = mp, fstype
    return best_fs in ("tmpfs", "ramfs")

def ensure_offtmp_tmpdir():
    current = os.environ.get("TMPDIR")

    if current and os.path.isdir(current) and not _is_tmpfs(current):
        chosen = current

    else:
        chosen = os.path.join(os.getcwd(), ".hyb2_tmp")
        log.warning("TMPDIR unset or tmpfs; using %s for temp files", chosen)

    os.makedirs(chosen, exist_ok=True)
    os.environ["TMPDIR"] = chosen
    tempfile.tempdir = chosen
    return chosen