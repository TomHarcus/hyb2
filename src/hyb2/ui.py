import os, sys
from tqdm import tqdm
import logging

log = logging.getLogger(__name__)

class Steps:
    def __init__(self):
        self.n = 0

    def start(self, desc: str) -> None:
        self.n += 1
        prefix = f"[{self.n}]"
        log.info(f"{prefix} {desc}")

def progress(fh, path: str, desc: str):

    bar = tqdm(total=os.path.getsize(path), unit="B", unit_scale=True,
               desc=f"  {desc}", disable=not sys.stderr.isatty(), leave=False)

    try:
        for line in fh:
            bar.update(len(line))
            yield line
    finally:
        bar.close()