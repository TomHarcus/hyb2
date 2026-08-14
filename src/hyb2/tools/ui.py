import os, sys
from tqdm import tqdm
import logging

import time, threading, itertools, contextlib

log = logging.getLogger(__name__)

class Steps:
    def __init__(self):
        self.n = 0

    def start(self, desc: str) -> None:
        self.n += 1
        prefix = f"[{self.n}]"
        log.info(f"{prefix} {desc}")

# progress bar for writing to files
def progress(fh, path: str, desc: str):

    bar = tqdm(total=os.path.getsize(path), unit="B", unit_scale=True,
               desc=f"  {desc}", disable=not sys.stderr.isatty(), leave=False)

    try:
        for line in fh:
            bar.update(len(line))
            yield line
    finally:
        bar.close()

# ellipsis loading for steps that have no way of knowing when finished
@contextlib.contextmanager
def spinner(desc: str):

    if not sys.stderr.isatty():
        log.info(f"  {desc}...")
        yield
        return

    bar = tqdm(total=None, bar_format="{desc}  {elapsed}", leave=False)
    loading = itertools.cycle(["|", "/", "-", "\\"])
    stop = threading.Event()

    def _load():
        while not stop.wait(0.4):
            bar.set_description_str(f"  {desc}{next(loading)}")

    t = threading.Thread(target=_load, daemon=True)
    t.start()

    try:
        yield
    finally:
        stop.set()
        t.join()
        bar.close()

# count bar for greedy constrained loop
def track(iterable, desc, total=None):

    return tqdm(iterable, total=total, desc=f"  {desc}",
                disable=not sys.stderr.isatty(), leave=False)