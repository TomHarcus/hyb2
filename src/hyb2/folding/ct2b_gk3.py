""" Port of ct2b_gk3.pl

Needed for the UNAFold foler

"""

def ct2b_gk3(in_ct, hybrid_ss_min):

    n = 1
    curr_line = 1
    out = []

    for line in in_ct.splitlines():
        cols = line.split()

        if "dG" in line:
            energy = seq = bp = ""
            n = int(cols[0])
            a = line.split("=", 1)
            if len(a) > 1:
                e = a[1].split()
                if e:
                    energy = f"\t({e[0]})"

            if not hybrid_ss_min:
                out.append(cols[7])
            elif hybrid_ss_min:
                out.append(cols[4])

            curr_line = 1

        else:
            seq += cols[1]

            if int(cols[4]) == 0:
                bp += "."
            elif int(cols[0]) < int(cols[4]):
                bp += "("
            else:
                bp += ")"

        if curr_line % (n + 1) == 0:
            out.append(seq)
            out.append(f"{bp}{energy}")


        curr_line += 1

    return "\n".join(out) + "\n"

