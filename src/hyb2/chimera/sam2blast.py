"""Adopt bin/sam2blast_3 

bin/sam2blast_3 is already Python 3 - this is an ADOPT, not a translation

Original script already written in Python, just some shape changes
"""

import sys
from typing import Iterable, Iterator
import re
from math import exp, log

_CIGAR_RE = re.compile(r'(\d+)(\D)')


def sam2blast(lines: Iterable[str]) -> Iterator[str]:
    dblen = 0
    for line in lines:
        #li = line.strip()
        if not line.startswith("@"):
            arr = line.split()
            if (int(arr[1]) in [0,256]):
                yield print_line(arr, 'f', dblen)
            elif (int(arr[1]) in [16,272]):
                yield print_line(arr, 'r', dblen)
            elif (int(arr[1]) == 4):
                pass
            else:
                raise Exception("unknown flag %s" % arr[1])
        elif line.startswith("@SQ"):
            dblen += int(line.split("LN:")[1])

def print_line(arr, flag, dblen):
    r"""
    cigar, mismatches, gaps, identity = [], 0, 0, 0
    cigar_num = re.findall(r'\d+',arr[5])
    cigar_op = re.findall(r'\D',arr[5])
    ref_len = 0
    """
    ref_len = read_len = gaps = 0
    lead_clip = trail_clip = 0
    started = False

    
    # instead of building an array of all bases to find where alignment starts and ends
    # just use running counters in a single pass
    for num, op in _CIGAR_RE.findall(arr[5]):
        n = int(num)
        if op == "M":
            ref_len += n
            read_len += n
            started = True
            trail_clip = 0
        elif op == "I":
            read_len += n
            gaps += 1
            started = True
            trail_clip = 0
        elif op == "D":
            ref_len += n
            gaps += 1
        elif op == "S" or op == "H":
            read_len += n
            if not started:
                lead_clip += n
            trail_clip += n
        else:
            raise Exception("unknown cigar operator %s" % op)
    align_start = lead_clip
    align_end = read_len - trail_clip
    len_read = read_len
    len_align = align_end - align_start


    ref_start = int(arr[3])
    if (flag == 'f'):
        ref_end = ref_start + ref_len -1
    elif (flag == 'r'): 
        ref_end = ref_start
        ref_start = ref_start + ref_len -1
        align_start, align_end = len_read - align_end, len_read - align_start

    # Takes the first NM/AS tag, not the last using the (is None guards)
    # matches legacy behaviour as without the guards it takes the last occurrence
    mismatches = align_score = None
    for x in arr[11:]:
        if mismatches is None and x.startswith("NM:i:"):
            mismatches = int(x[5:])
        elif align_score is None and x.startswith("AS:i:"):
            align_score = int(x[5:])
        if align_score is not None and mismatches is not None:
            break

    identity = ((abs(align_end - align_start) - mismatches )/ 
                       float(len_align) ) * 100


    #ungapped l, k, h, n, m = 1.33, 0.621, 1.12, dblen, len_read
    l, k, h, n, m = 1.28, 0.46, 0.85, dblen, len_read
    np = n-log(k*n*m)/h
    mp = m-log(k*n*m)/h
    if mp < 0:
        mp = 0.1
    bit_score = ((l*float(align_score))-log(k))/log(2)
    # alternative evalue formula, gives 0's for small numbers
    # u = (log(k*mp*np))/l
    # evalue = 1-exp(-exp(-l*(float(align_score)-u)))
    evalue = mp*np*2**-bit_score
    return ('%s\t%s\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.2g\t%.1f\t%s'%(
           arr[0], arr[2], identity, len_align, 
           mismatches, gaps, align_start+1, align_end, ref_start, ref_end,
           evalue, bit_score, arr[9]) + "\n")

    


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: sam2blast <in.sam>", file=sys.stderr)
        return 1
    with open(argv[0]) as f:
        for out in sam2blast(f):
            sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
