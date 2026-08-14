""" Port of comradesFold2

Greedily folds an RNA fragment, adding the ranked base-pair constraints one at a time
and keeping only those the fold actually honours

outputs: .vienna + .ct

"""

import subprocess, random, shutil, sys, argparse, os
from pathlib import Path

from hyb2.folding.ct2b_gk3 import ct2b_gk3
from hyb2.tools import config

import logging
from hyb2.tools import ui

log = logging.getLogger(__name__)

def run(in_constraints, in_fasta, *, output_id=None, shuffling=False, fold="vienna",
        vienna_bin=None, basepair_scores=None, begin=None, end=None, alpha=config.CPL_DEFAULTS["alpha"],
        beta=config.CPL_DEFAULTS["beta"], normalize=config.CPL_DEFAULTS["normalize"], beam_size=config.CPL_DEFAULTS["beam_size"],
        energy_delta=config.CPL_DEFAULTS["energy_delta"], max_phase1=config.CPL_DEFAULTS["max_phase1"], 
        max_phase2=config.CPL_DEFAULTS["max_phase2"], energy_model=config.CPL_DEFAULTS["energy_model"]):
    
    """
    cluster array job code would live here
    """

    if fold is None:
            fold = "cplfold"

    # input fold mapping to preserve legacy inputs: 0 and 1
    fold = {1: "vienna", "1": "vienna", 0: "unafold", "0": "unafold"}.get(fold, fold)

    shuffled_constraints = f"{in_constraints}.shuf"
    current_constraints = f"{in_fasta}.aux"
    ct_output = f"{in_fasta}.ct"
    vienna_output = f"{in_fasta}.vienna"

    log.debug(f"Input fasta file: {in_fasta}")

    if fold == "cplfold":
        if begin is None or end is None:
            raise ValueError("cplfold needs begin + end")
        
        _fold_cplfold(in_fasta, basepair_scores, begin, end, ct_output, vienna_output,
                      alpha=alpha, beta=beta, normalize=normalize, beam_size=beam_size,
                      energy_delta=energy_delta, max_phase1=max_phase1, max_phase2=max_phase2,
                      energy_model=energy_model, vienna_bin=vienna_bin)
        
        if output_id:
            for path in (ct_output, vienna_output):
                shutil.copy(path, f"{path}.{output_id}")
        
        return (vienna_output, ct_output)

    log.info(f"Input constraints file: {in_constraints}")
    log.info(f"Creating constraints file: {current_constraints}")
    log.info(f"Creating shuffled constraints file: {shuffled_constraints}")

    if current_constraints == in_constraints:
        raise ValueError("name your constraint file something else... Exiting")

    with open(in_constraints, "r") as fin, open(shuffled_constraints, "w") as fout: 
        lines = fin.readlines()
        if shuffling:
            random.shuffle(lines)

        fout.writelines(lines)
    
    current = []
    accepted = []
    for constraint in ui.track(lines, "fitting constraints", total=len(lines)):
        
        current.append(constraint)

        with open(current_constraints, "w") as fout:
            fout.writelines(current)

        if fold == "vienna":
            _fold_vienna_constrained(in_fasta, current_constraints, ct_output, vienna_bin)
        elif fold == "unafold":
            _fold_unafold(in_fasta, ct_output, vienna_output=None)
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
        _fold_unafold(in_fasta, ct_output, vienna_output=vienna_output)
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

    # call RNAfolder
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
                  alpha, beta, normalize, beam_size, energy_delta, max_phase1,
                  max_phase2, energy_model, vienna_bin):
    import numpy as np


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

    matrix = None

    # construct cplfold support matrix
    if basepair_scores is not None:
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

                # normalize or raw matrix
                v = np.log1p(count) if normalize == "log" else count

                matrix[i-begin, j-begin] = v
                matrix[j-begin, i-begin] = v

    # call cplfold
    with ui.spinner("folding (cplfold) "):
        results = two_phase_pseudoknot_fold(
            seq,
            bonus_matrix=matrix,
            alpha=alpha,
            beta=beta,
            beam_size=beam_size,
            energy_delta=energy_delta,
            max_phase1=max_phase1,
            max_phase2=max_phase2,
            energy_model=energy_model,
            verbose=False
        )

    if not results:
        raise ValueError("CPLfold returned no structures")
    
    best = results[0]
    structure = best["structure"]
    energy = best.get("energy")

    # CPLfold returns energy=None when HotKnots' computeEnergy fails (it catches
    # the error, warns, and leaves energy unset). The commonest cause is an
    # architecture mismatch in the compiled binary -- e.g. an aarch64 build run
    # on x86-64 gives "exec format error". Silently writing 0.0 here would emit
    # a structure with a meaningless energy and corrupt the COMRADES-score
    # ranking, so fail loudly with the fix instead.
    if energy is None:
        hk = os.path.join(config.CPLFOLD_DIR, "Utils", "HotKnots_v2.0")
        raise RuntimeError(
            f"CPLfold folded {gene_name!r} but HotKnots could not compute its "
            f"energy. This usually means the compiled binary {hk}/bin/computeEnergy "
            f"does not match this machine's architecture. Rebuild it:\n"
            f"    make -C {hk}\n"
            f"then verify: `uname -m` vs `file {hk}/bin/computeEnergy`."
        )

    vienna = f">{gene_name}\n{seq}\n{structure} ({energy})\n"

    if vienna_output is not None:
        with open(vienna_output, "w") as fout:
            fout.write(vienna)

    b = (vienna_bin.rstrip("/") + "/") if vienna_bin else ""

    ct = _dot_to_ct(gene_name, seq, structure, energy)

    with open(ct_output, "w") as fout:
        fout.write(ct)

def _fold_unafold(in_fasta, ct_output, vienna_output):

    # call unafold
    subprocess.run(
        ["hybrid-ss-min", "-c", in_fasta],
        capture_output=True, text=True, check=True,
        env={**os.environ, "UNAFOLDDAT": str(config.UNAFOLD_DIR)}
    )

    result = f"{in_fasta}.ct"
    if result != ct_output:
        shutil.move(result, ct_output)   

    if vienna_output is not None:
        Path(vienna_output).write_text(ct2b_gk3(Path(ct_output).read_text(),
                                                hybrid_ss_min=True))

    

def _dot_to_bpmap(dot):
    """Dot-bracket -> {position: partner} pseudoknot aware pair map. Crossing pairs
    resolve correctly (b2ct can't). Adapted from IPyRSSA's Structure.dot2bpmap
    Link: https://github.com/lipan6461188/IPyRSSA -- per Ke Wang's suggestion
    """
    stack, bpmap = [], {}

    for idx, sym in enumerate(dot):
        if sym in "([{<":
            stack.append((idx + 1, sym))
        elif sym in ".-_=:,":
            continue
        else:
            for i in range(len(stack) - 1, -1, -1):
                if stack[i][1] + sym in ("()", "[]", "{}", "<>"):
                    j = stack.pop(i)[0]
                    bpmap[idx + 1] = j
                    bpmap[j] = idx + 1
                    break

    return bpmap
    
def _dot_to_ct(name, seq, structure, energy):
    """Writes a pseudoknot-aware .ct (b2ct silently drops the crossing pairs).
    VARNA renders the crossing from this .ct
    """
    bpmap = _dot_to_bpmap(structure)
    n = len(seq)

    lines = [f"{n:5d} dG = {energy}\t{name}\n"]
    for i in range(1, n+1):
        pair = bpmap.get(i, 0)
        lines.append(f"{i:5d} {seq[i-1]} {i-1:7d} {i+1 if i < n else 0:4d} {pair:4d} {i:4d}\n")

    return "".join(lines)

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
    p.add_argument("-r", dest="fold", default="vienna", choices=["vienna", "unafold", "cplfold", "1", "0"], help="folding backend: (1=vienna, 0=unafold are legacy aliases)")
    p.add_argument("-p", dest="basepair_scores", default=None, metavar="BASEPAIR_SCORES.TXT", help="cplfold support matrix (i j count); required for -r cplfold")
    p.add_argument("-b", dest="begin", type=int, default=None, help="cplfold window start (required for -r cplfold)")
    p.add_argument("-e", dest="end", type=int, default=None, help="cplfold window end (required for -r cplfold)")
    p.add_argument("--alpha", dest="alpha", type=float, default=config.CPL_DEFAULTS["alpha"], help="cplfold bonus weight")
    p.add_argument("--beta", dest="beta", type=float, default=config.CPL_DEFAULTS["beta"], help="cplfold bonus weight")
    p.add_argument("--normalize", dest="normalize", choices=["raw", "log"], default=config.CPL_DEFAULTS["normalize"], help="cplfold bonus normalization")
    p.add_argument("--beam-size", dest="beam_size", type=int, default=config.CPL_DEFAULTS["beam_size"], help="cplfold beam size")
    p.add_argument("--energy-delta", dest="energy_delta", type=float, default=config.CPL_DEFAULTS["energy_delta"], help="cplfold energy delta")
    p.add_argument("--max-phase1", dest="max_phase1", type=int, default=config.CPL_DEFAULTS["max_phase1"], help="cplfold max phase 1")
    p.add_argument("--max-phase2", dest="max_phase2", type=int, default=config.CPL_DEFAULTS["max_phase2"], help="cplfold max phase 2")
    p.add_argument("--energy-model", dest="energy_model", choices=["DP09", "DP03", "CC06", "CC09", "RE"], default=config.CPL_DEFAULTS["energy_model"], help="cplfold energy model")


    return p

def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run(
        args.in_constraints,
        args.in_fasta,
        output_id=args.output_id,
        shuffling=(args.shuffling=="1"),
        fold=args.fold,
        basepair_scores=args.basepair_scores,
        begin=args.begin,
        end=args.end,
        alpha=args.alpha,
        beta=args.beta,
        normalize=args.normalize,
        beam_size=args.beam_size,
        energy_delta=args.energy_delta,
        max_phase1=args.max_phase1,
        max_phase2=args.max_phase2,
        energy_model=args.energy_model
    )

    return 0

if __name__ == "__main__":
    sys.exit(main())