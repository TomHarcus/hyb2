""" Port of comradesFold2

Greedily folds an RNA fragment, adding the ranked base-pair constraints one at a time
and keeping only those the fold actually honours

outputs: .vienna + .ct

"""

import subprocess, random, shutil, sys, argparse

def run(in_constraints, in_fasta, *, output_id=None, shuffling=False, fold="vienna",
        vienna_bin=None, basepair_scores=None, begin=None, end=None, alpha=0.5,
        beta=0.0, normalize="log", beam_size=100):
    
    """
    cluster array job code would live here
    """

    shuffled_constraints = f"{in_constraints}.shuf"
    current_constraints = f"{in_fasta}.aux"
    ct_output = f"{in_fasta}.ct"
    vienna_output = f"{in_fasta}.vienna"

    if fold == "cplfold":
        if basepair_scores is None or begin is None or end is None:
            raise ValueError("cplfold needs basepair_scores + begin/end")
        
        _fold_cplfold(in_fasta, basepair_scores, begin, end, ct_output, vienna_output,
                      alpha=alpha, beta=beta, normalize=normalize, beam_size=beam_size,
                      vienna_bin=vienna_bin)
        
        if output_id:
            for path in (ct_output, vienna_output):
                shutil.copy(path, f"{path}.{output_id}")
        
        return (vienna_output, ct_output)

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

def _fold_cplfold(in_fasta, basepair_scores, begin, end, ct_output, vienna_output,
                  alpha, beta, normalize, beam_size, vienna_bin):
    import numpy as np
    from hyb2 import config

    if config.CPLFOLD_DIR not in sys.path:
        sys.path.insert(0, config.CPLFOLD_DIR)
    
    from CPLfold import two_phase_pseudoknot_fold

    with open(in_fasta) as f:
        lines = f.readlines()
        gene_name = lines[0][1:].rstrip()
        seq = lines[1].rstrip().replace("T", "U")
        n = len(seq)

        if n != end - begin + 1:
                raise ValueError("wrong format")

    matrix = np.zeros((n, n))

    with open(basepair_scores) as f:
        lines = f.readlines()

        for line in lines:
            elements = line.split()

            if len(elements) < 3:
                continue

            i, j, count = int(elements[0]), int(elements[1]), float(elements[2])

            if not (begin <= i <= end and begin <= j <= end):
                continue
            
            v = np.log1p(count) if normalize == "log" else count

            matrix[i-begin, j-begin] = v
            matrix[j-begin, i-begin] = v

    results = two_phase_pseudoknot_fold(
        seq,
        bonus_matrix=matrix,
        alpha=alpha,
        beta=beta,
        beam_size=beam_size,
        verbose=False
    )

    if not results:
        raise ValueError("CPLfold returned no structures")
    
    best = results[0]
    structure = best["structure"]
    energy = best.get("energy")

    if energy is None:
        energy = 0.0

    vienna = f">{gene_name}\n{seq}\n{structure} ({energy})\n"

    if vienna_output is not None:
        with open(vienna_output, "w") as fout:
            fout.write(vienna)

    b = (vienna_bin.rstrip("/") + "/") if vienna_bin else ""

    # b2ct only understands nested () -- it silently emits nothing on pseudoknot
    # brackets. Strip [ ] -> . so the .ct carries the nested pairs (crossings are
    # kept in the full .vienna above)
    nested = structure.replace("[", ".").replace("]", ".")
    ct_input = f">{gene_name}\n{seq}\n{nested} ({energy})\n"

    ct = subprocess.run([b + "b2ct"], input=ct_input,
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