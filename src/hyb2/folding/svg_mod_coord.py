""" Port of svg_mod_coord.sh

Relabels a VARNA structure SVG's position numbers to genome coordinates
(shift by x_coord-1; two-strand case remaps across the length/gap boundaries)
and thins the stroke width. Legacy branch 2 (both X and Y unset) is unreachable
and not ported.
"""

import argparse
import re
import sys

def svg_mod_coord(in_file, x_coord, y_coord, length, prefix):

    X = x_coord - 1

    temp_file = in_file[:-4] if in_file.endswith(".svg") else in_file
    if prefix is None:
        out_file = f"{temp_file}_plot.svg"

    else:
        out_file = f"{prefix}.{temp_file}_plot.svg"

    with open(in_file) as f:
        lines = f.readlines()

    result = []

    for line in lines:
        fields = line.split()

        if is_position_number(fields):
            b2 = int(extract_number(fields[9]))

            if y_coord is None:
                new_tag = f">{b2 + X}</text>"
            else:
                Y = y_coord - 1
                L = length
                if b2 <= L:
                    new_tag = f">{b2 + X}</text>"
                elif b2 <= L + 100:
                    new_tag = fields[9]
                else:
                    new_tag = f">{b2 + Y - L - 100}</text>"

            fields.append(new_tag)
            fields[9] = " "
            result.append(" ".join(fields))
        
        else:
            result.append(line.rstrip("\n"))

    text = "\n".join(result) + "\n"
    text = text.replace('stroke-width="1.0"', 'stroke-width="0.25"')

    with open(out_file, "w") as fout:
        fout.write(text)

    return out_file


def is_position_number(fields):
    if len(fields) < 10:
        return False
    
    f6, f10 = fields[5], fields[9]
    return ("7.5" in f6) and ("text" in f10) and (re.search(r">[A-Z]<", f10) is None)

def extract_number(f10):
    return f10.split("<")[0].split(">")[1]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="svg-mod-coord",
        description="Relabel a VARNA SVG's position numbers to genome coords "
        "(port of svg_mod_coord.sh).",
    )
    p.add_argument("-i", dest="in_file", required=True, metavar="IN.svg", help="VARNA SVG (required)")
    p.add_argument("-x", dest="x_coord", type=int, required=True, help="start coordinate of 1st strand")
    p.add_argument("-y", dest="y_coord", type=int, default=None, help="start coordinate of 2nd strand (two-strand mode)")
    p.add_argument("-l", dest="length", type=int, default=None, help="length of 1st strand (two-strand mode)")
    p.add_argument("-p", dest="prefix", default=None, help="optional output-name prefix")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    svg_mod_coord(args.in_file, args.x_coord, args.y_coord, args.length, args.prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())