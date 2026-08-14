""" Port of plot_VARNA

Renders a folded .ct with VARNA (coloured by per-base support scores), then
relabels the SVG coordinates via svg_mod_coord. Interactive mode (-0 1) pops up
the VARNA GUI; headless mode writes <prefix>.<name>_plot.svg. The WSL/Windows
branch of the legacy is dropped (Linux-only). Intermediate SVG is preserved
(the legacy rm's it).
"""

import argparse
import subprocess
import sys
import contextlib

from hyb2.folding.svg_mod_coord import svg_mod_coord
from hyb2.tools import config
from hyb2.tools.logsetup import is_quiet
from hyb2.tools.ui import spinner

import logging

log = logging.getLogger(__name__)

def plot_VARNA(in_file, score, VARNA=None, *, x_coord, y_coord, length, interactive):

    quiet = is_quiet()

    if VARNA is None:
        VARNA = config.VARNA_JAR

    with open(in_file, "r") as f:
        ct_lines = f.readlines()

    header = ct_lines[0].split()[4]
    scores = [line.strip() for line in open(score) if line.rstrip() != ""]
    colorscores = ",".join(scores)
    smax = max(scores, key=float)
    sm = float(smax)

    colorstyle = f"0.00:#FFFFFF,{sm*0.2:g}:#4747FF,{sm*0.55:g}:#1CFF47,{sm*0.65:g}:#FFFF47,{sm*0.9:g}:#FF4747,{sm}:#B64747"
    
    varna_cmd = [
        "java", "-jar", VARNA,
        "-i", in_file,
        "-bpStyle", "simple",
        "-colorMap", colorscores,
        "-colorMapMin", "0",
        "-colorMapMax", smax,
        "-colorMapStyle", colorstyle,
        "-title", header,
        "-spaceBetweenBases", "0.6"
    ]

    if interactive:
        subprocess.run(varna_cmd, stdout=subprocess.DEVNULL if quiet else None, stderr=subprocess.DEVNULL if quiet else None)
        log.info("Modify base numbers: svg_mod_coord -i <svg> -x start_coord -y 2nd_strand_coord(only if it exists) -l length_of_fragment")
        return None
    
    else:
        svg = in_file[:-3] + ".svg" if in_file.endswith(".ct") else in_file + ".svg"
        with spinner("rendering structure (VARNA) ") if quiet else contextlib.nullcontext():
            if not quiet:
                print("rendering structure (VARNA)")
            subprocess.run(varna_cmd + ["-o", svg],
                        stdout=subprocess.DEVNULL if quiet else None, stderr=subprocess.DEVNULL if quiet else None)

        prefix = score.split("__")[0]

        out = svg_mod_coord(svg, x_coord, y_coord, length, prefix)

        return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plot-VARNA",
        description="Render a folded .ct with VARNA and relabel coordinates "
        "(port of plot_VARNA).",
    )
    p.add_argument("-i", dest="in_file", required=True, metavar="IN.ct", help="folded structure (.ct)")
    p.add_argument("-s", dest="score", required=True, metavar="VARNA_scores.txt", help="per-base colour scores")
    p.add_argument("-j", dest="varna", default=None, metavar="VARNAcmd.jar", help="VARNA jar (default: config.varna_jar())")
    p.add_argument("-x", dest="x_coord", type=int, default=None, help="start coordinate of 1st strand")
    p.add_argument("-y", dest="y_coord", type=int, default=None, help="start coordinate of 2nd strand")
    p.add_argument("-l", dest="length", type=int, default=None, help="length of 1st strand")
    p.add_argument("-0", dest="interactive", default=None, help="1 = VARNA interactive GUI mode")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    plot_VARNA(
        args.in_file,
        args.score,
        args.varna,
        x_coord=args.x_coord,
        y_coord=args.y_coord,
        length=args.length,
        interactive=(args.interactive == "1"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
