""" Port of hyb2constraints.pl

Takes hyb file as input, predicte RNA stem coordinates
Output file in form "F 1 100 5", which can be used as a list of constraints for hybrid-ss-min

***UNTESTED RIGHT NOW***

"""

def hyb2constraints(file):
    out = []
    for line in file:
        columns = line.rstrip("\n").split("\t")
        arm1_start, arm1_end = int(columns[6]), int(columns[7])
        arm2_start, arm2_end = int(columns[12]), int(columns[13])

        if overlap(arm1_start, arm1_end, arm2_start, arm2_end) > 0:
            continue
        
        coords = [arm1_start, arm1_end, arm2_start, arm2_end]
        sorted_coords = sorted(coords)

        out.append(
            f"F\t{sorted_coords[0]}\t{sorted_coords[3]}\t{sorted_coords[1]-sorted_coords[0]+1}"
            )

    return "\n".join(out) + "\n"


def overlap(a, b, c, d):
    return 1 + min(b, d) - max(a, c)