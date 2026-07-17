""" Port of comradesFold2

Greedily folds an RNA fragment, adding the ranked base-pair constraints one at a time
and keeping only those the fold actually honours

outputs: .vienna + .ct

"""

import subprocess, random, shutil, sys, argparse

def run(in_constraints, in_fasta, *, output_id=None, shuffling=False, fold="vienna",
        vienna_bin=None):
    
    """
    cluster array job code would live here
    """

    shuffled_constraints = f"{in_constraints}.shuf"
    current_constraints = f"{in_fasta}.aux"
    ct_output = f"{in_fasta}.ct"
    vienna_output = f"{in_fasta}.vienna"

    if current_constraints == in_constraints:
        raise ValueError("name your constraint file something else")

    with open(in_constraints, "r") as fin, open(shuffled_constraints, "w") as fout: 
        lines = fin.readlines()
        if shuffling:
            random.shuffle(lines)

        fout.writelines(lines)
    
    current = []
    accepted = []
    for constraint in lines:
        
        current.append(constraint)

        with open(current_constraints, "w") as fout:
            fout.writelines(current)

        if fold == "vienna":
            _fold_vienna_constrained(in_fasta, current_constraints, ct_output, vienna_bin)
        elif fold == "unafold":
            raise NotImplementedError("not implemented yet")
        elif fold == "cplfold":
            raise NotImplementedError("not implemented yet")
        else:
            raise ValueError("not a valid folding algorithm")
        
        _, i, j, _ = current[-1].split()
        satisfied = False

        with open(ct_output) as fin:
            for line in fin:
                cols = line.split()
                if len(cols) >= 5 and cols[0] == i and cols[4] == j:
                    satisfied = True
                    break
        
        if satisfied:
            accepted = current[:]
        else:
            current = accepted[:]
    
    with open(current_constraints, "w") as fout:
        fout.writelines(accepted)

    if fold == "vienna":
        _fold_vienna_constrained(in_fasta, current_constraints, ct_output, vienna_bin, vienna_output=vienna_output)
    elif fold == "unafold":
        raise NotImplementedError("not implemented yet")
    elif fold == "cplfold":
        raise NotImplementedError("not implemented yet")
    else:
        raise ValueError("not a valid folding algorithm")
    
    if output_id:
        for path in (ct_output, current_constraints, vienna_output):
            shutil.copy(path, f"{path}.{output_id}")

    return (vienna_output, ct_output)

def _fold_vienna_constrained(in_fasta, constraints_file, ct_output, vienna_bin, vienna_output=None):
    b = (vienna_bin.rstrip("/") + "/") if vienna_bin else ""

    with open(in_fasta) as fin:
        fasta_text = fin.read()

    fold = subprocess.run(
        [b + "RNAfold", "--noconv", "--noPS", f"--commands={constraints_file}"],
        input=fasta_text, capture_output=True, text=True
    ).stdout

    vienna = fold.replace("&>", "-").replace("&", "")

    if vienna_output is not None:
        with open(vienna_output, "w") as fout:
            fout.write(vienna)

    ct = subprocess.run([b + "b2ct"], input=vienna,
                        capture_output=True, text=True).stdout
    
    ct = ct.replace("ENERGY =", "dG =")

    with open(ct_output, "w") as fout:
        fout.write(ct)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="comrades-fold",
        description="constraint-guided RNA folding",
        add_help=False
    )

    p.add_argument("--help", action="help", help="Show this help message and exit")
    p.add_argument("-c", dest="in_constraints", required=True, metavar="INPUT.TXT", help="Input TXT (required)")
    p.add_argument("-i", dest="in_fasta", required=True, metavar="INPUT.FASTA", help="Input FASTA (required)")
    p.add_argument("-o", dest="output_id", default=None, help="output file name")
    p.add_argument("-s", dest="shuffling", help="shuffle constraints")
    p.add_argument("-r", dest="fold", type=str, default="vienna", help="folder to be used")

    return p

def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run(
        args.in_constraints,
        args.in_fasta,
        output_id=args.output_id,
        shuffling=(args.shuffling=="1"),
        fold=args.fold
    )

    return 0

if __name__ == "__main__":
    sys.exit(main())