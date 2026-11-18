""" Port of ct2bps_2.awk

Takes in a .ct file

Parses folding output where each chimeric reads header encodes the 
genomic coordinates of its two arms, then converts local fold
coordinates back to those absolute coordinates

Two deliberate departures from the legacy awk (pass legacy=True to reproduce it):
  * pairs formed within one arm (partner also in arm 1) are skipped. The awk mapped
    them through the arm-2 offset, producing spurious pairs (~1% of pairs on U2).
  * each pair is written as (lower, higher) coordinate. Chimeras are sequenced in
    both arm orders, so the awk reported the same base pair as "i j" or "j i" and the
    support for one pair was split between two keys downstream.

"""

import sys

def ct2bps_2(ct_file, legacy=False):

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

            if not legacy and pairing_partner <= bit1_len:
                continue   # intra-arm pair, not evidence for the chimeric interaction

            current_position += int(bit1_st) - 1
            pairing_partner = pairing_partner - bit1_len + int(bit2_st) - 1

            if not legacy and pairing_partner < current_position:
                current_position, pairing_partner = pairing_partner, current_position

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