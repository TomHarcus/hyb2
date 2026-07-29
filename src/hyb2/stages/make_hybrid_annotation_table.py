""" Port of make_hybrid_annotation_table.pl

"""

def make_hybrid_annotation_table(files, SORT_COLUMNS=0, OUTPUT_ZERO_FOR_NA=0):

    all_hyb_names = {}
    all_hyb_properties = {}
    data = {}

    for path in files:
        with open(path) as f:
            for line in f:
                elements = line.rstrip("\n").split("\t")
                if len(elements) < 2:
                    continue

                name, properties = elements[0], elements[1]
                curr_hyb_properties = properties.split(";")

                if name not in all_hyb_names:
                    all_hyb_names[name] = 1

                for prop in curr_hyb_properties:
                    if "=" not in prop:
                        break

                    key, val = prop.split("=", 1)
                    if key not in all_hyb_properties:
                        all_hyb_properties[key] = 1
                    data.setdefault(name, {})[key] = val
               
    columns = sorted(all_hyb_properties) if SORT_COLUMNS else list(all_hyb_properties)

    out = ["#seq_ID" + "".join(f"\t{p}" for p in columns)]

    for name in all_hyb_names:
        row = name
        for p in columns:
            val = data.get(name, {}).get(p)
            row += f"\t{val}" if val is not None else ("\t0" if OUTPUT_ZERO_FOR_NA else "\tNA")

        out.append(row)

    return "\n".join(out) + "\n"

 
