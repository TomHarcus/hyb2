""" Port of ct2bps_2.awk

Takes in a .ct file

Parses folding output where each chimeric reads header encodes the 
genomic coordinates of its two arms, then converts local fold
coordinates back to those absolute coordinates

"""

import sys

def ct2bps_2(ct_file):

    for line in ct_file.splitlines():
        columns = line.split()

        if "dG" in line:
            last_field = columns[-1]

            a = last_field.split("-")
            n = len(a)

            i = a[n // 2 - 1].split("_")
            j = a[-1].split("_")

            bit1_st = i[-2]
            bit1_en = i[-1]

            bit2_st = j[-2]
            #bit2_en = j[-1]

            bit1_len = int(bit1_en) - int(bit1_st) + 1

        else:
            pairing_partner = int(columns[4])
            current_position = int(columns[0])

            if pairing_partner == 0 or current_position > bit1_len:
                continue

            current_position += int(bit1_st) - 1
            pairing_partner = pairing_partner - bit1_len + int(bit2_st) - 1

            yield f"{current_position}\t{pairing_partner}\n"


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: ct2bps_2 <ct_file>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        sys.stdout.writelines(ct2bps_2(f.read()))
    return 0


if __name__ == "__main__":
    sys.exit(main())