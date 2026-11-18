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
from hyb2.tools.logsetup import is_quiet
import contextlib

log = logging.getLogger(__name__)

def run(in_constraints, in_fasta, *, output_id=None, shuffling=False, fold="vienna",
        vienna_bin=None, basepair_scores=None, begin=None, end=None, alpha=config.CPL_DEFAULTS["alpha"],
        beta=config.CPL_DEFAULTS["beta"], normalize=config.CPL_DEFAULTS["normalize"], beam_size=config.CPL_DEFAULTS["beam_size"],
        energy_delta=config.CPL_DEFAULTS["energy_delta"], max_phase1=config.CPL_DEFAULTS["max_phase1"], 
        max_phase2=config.CPL_DEFAULTS["max_phase2"], allow_pseudoknot=config.CPL_DEFAULTS["allow_pseudoknot"], energy_model=config.CPL_DEFAULTS["energy_model"],
        support=config.CPL_DEFAULTS["support"], support_hyb=None):
    
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
                      support=support, support_hyb=support_hyb, alpha=alpha, beta=beta, normalize=normalize, beam_size=beam_size,
                      energy_delta=energy_delta, max_phase1=max_phase1, max_phase2=max_phase2,
                      allow_pseudoknot=allow_pseudoknot, energy_model=energy_model, vienna_bin=vienna_bin)
        
        if output_id:
            for path in (ct_output, vienna_output):
                shutil.copy(path, f"{path}.{output_id}")
        
        return (vienna_output, ct_output)

    log.debug(f"Input constraints file: {in_constraints}")
    log.debug(f"Creating constraints file: {current_constraints}")
    log.debug(f"Creating shuffled constraints file: {shuffled_constraints}")

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
        elif fold == "mfold":
            _fold_mfold(in_fasta, current_constraints, ct_output, vienna_output=None)
        else:
            raise ValueError("not a valid folding algorithm")
        
        _, i, j, _ = current[-1].split()
        satisfied = False

        if not os.path.exists(ct_output):   # mfold writes no .ct when the constraints are infeasible
            current = accepted[:]
            continue

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

    # keep a readable record of which constraints the fold accepted
    accepted_log = f"{in_constraints}.accepted"
    with open(accepted_log, "w") as fout:
        fout.writelines(accepted)
    log.info(f"{len(accepted)} of {len(lines)} constraints accepted ({fold}); list: {accepted_log}")

    if fold == "vienna":
        _fold_vienna_constrained(in_fasta, current_constraints, ct_output, vienna_bin, vienna_output=vienna_output)
    elif fold == "unafold":
        _fold_unafold(in_fasta, ct_output, vienna_output=vienna_output)
    elif fold == "mfold":
        _fold_mfold(in_fasta, current_constraints, ct_output, vienna_output=vienna_output)
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
    proc = subprocess.run(
        [b + "RNAfold", "--noconv", "--noPS", f"--commands={constraints_file}"],
        input=fasta_text, capture_output=True, text=True
    )

    # exit 1 + "empty solution set" is the expected over constrained rejection
    # the greedy loop detects it via the absent pair. Any other failure is real
    if proc.returncode != 0 and "empty solution set" not in proc.stderr:
        raise RuntimeError(f"RNAfold failed (exit {proc.returncode}): {proc.stderr.strip()}")
    
    fold = proc.stdout

    vienna = fold.replace("&>", "-").replace("&", "")

    if vienna_output is not None:
        with open(vienna_output, "w") as fout:
            fout.write(vienna)

    ct = subprocess.run([b + "b2ct"], input=vienna,
                        capture_output=True, text=True, check=True).stdout
    
    ct = ct.replace("ENERGY =", "dG =")

    with open(ct_output, "w") as fout:
        fout.write(ct)

def _fold_cplfold(in_fasta, basepair_scores, begin, end, ct_output, vienna_output,
                  alpha, beta, normalize, beam_size, energy_delta, max_phase1,
                  max_phase2, allow_pseudoknot, energy_model, vienna_bin,
                  support="cofold", support_hyb=None):
    quiet = is_quiet()

    if config.CPLFOLD_DIR not in sys.path:
        sys.path.insert(0, config.CPLFOLD_DIR)

    from CPLfold import two_phase_pseudoknot_fold, compute_energy_hotknots, HotKnotsEnergy
    from CPLfold_parser import BeamCKYParserHyper

    with open(in_fasta) as f:
        lines = f.readlines()
        gene_name = lines[0][1:].rstrip()
        seq = lines[1].rstrip().replace("T", "U")
        n = len(seq)

        if n != end - begin + 1:
            raise ValueError("wrong format")

    matrix = support_matrix(support, n, begin, end, normalize,
                            basepair_scores=basepair_scores, support_hyb=support_hyb)

    with ui.spinner("folding (cplfold) ") if quiet else contextlib.nullcontext():
        if allow_pseudoknot:
            results = two_phase_pseudoknot_fold(
                seq,
                bonus_matrix=matrix,
                alpha=alpha,
                beta=beta,
                beam_size=beam_size,
                energy_delta=energy_delta,
                max_phase1=max_phase1,
                max_phase2=max_phase2,
                allow_pseudoknot=True,
                energy_model=energy_model,
                verbose=not quiet
            )
            if not results:
                raise ValueError("CPLfold returned no structures")
            structure, energy = results[0]["structure"], results[0].get("energy")
        else:
            # Without pseudoknots only the single best structure is needed. Calling the
            # parser directly skips the suboptimal-structure step (an outside pass written
            # in plain Python, ~99% of the runtime on long RNAs); the fold itself is the
            # same numba-compiled LinearFold beam search.
            parser = BeamCKYParserHyper(beam_size=beam_size, lv=True)
            if matrix is not None:
                parser.set_alpha(alpha)
                parser.set_bonus_matrix(matrix, n)
            structure = parser.parse(seq)[0]
            if not structure:
                raise ValueError("CPLfold returned no structures")
            energy = compute_energy_hotknots(seq, structure, HotKnotsEnergy(), energy_model)

    if energy is None:
        raise ValueError(f"CPLfold returned no energy for {gene_name}")

    vienna = f">{gene_name}\n{seq}\n{structure} ({energy})\n"

    if vienna_output is not None:
        with open(vienna_output, "w") as fout:
            fout.write(vienna)

    ct = _dot_to_ct(gene_name, seq, structure, energy)

    with open(ct_output, "w") as fout:
        fout.write(ct)


def support_matrix(support, n, begin, end, normalize, basepair_scores=None, support_hyb=None):
    """CPLfold bonus matrix (n x n, symmetric) from chimera data.

    cofold   : base pairs predicted by cofolding each chimera's two arms (basepair_scores
               file, "i j count", count = number of chimeras predicting that pair). Counts
               for i-j and j-i are added together.
    coverage : every position pair spanned by a chimera's two arms gets +1 per chimera.
    gaussian : IRIS/CPLfold PARIS model (CPLfold/extract_paris_scores.py): each chimera adds
               the outer product of two normal curves centred on its arms, peak-normalised,
               sd = mean arm length / 6.
    Chimeras are counted once each (read counts are not used), as in the cofold method.
    normalize: "log" -> log(x + 1), "raw" -> x.
    Returns None when there is no data to use.
    """
    import numpy as np

    raw = np.zeros((n, n))

    if support == "cofold":
        if basepair_scores is None:
            return None
        with open(basepair_scores) as f:
            for line in f:
                elements = line.split()
                if len(elements) < 3:
                    continue
                i, j, count = int(elements[0]), int(elements[1]), float(elements[2])
                if not (begin <= i <= end and begin <= j <= end):
                    continue
                raw[i - begin, j - begin] += count
                if i != j:
                    raw[j - begin, i - begin] += count

    elif support in ("coverage", "gaussian"):
        if support_hyb is None:
            raise ValueError(f"--support {support} needs the chimeras (.hyb) for this window")
        arms = []
        with open(support_hyb) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 14:
                    continue
                a = (int(c[6]) - begin, int(c[7]) - begin)
                b = (int(c[12]) - begin, int(c[13]) - begin)
                if min(a + b) < 0 or max(a + b) > n - 1:
                    continue
                arms.append((a, b))

        if support == "coverage":
            for (s1, e1), (s2, e2) in arms:
                raw[s1:e1 + 1, s2:e2 + 1] += 1
                raw[s2:e2 + 1, s1:e1 + 1] += 1
            raw[np.diag_indices(n)] /= 2      # an arm overlapping itself was counted twice
        elif arms:
            sd = np.mean([e - s for arm in arms for s, e in arm]) / 6
            sd = max(sd, 1e-3)
            half = int(np.ceil(6 * sd)) + 1   # beyond 6 sd the curve is < 1e-7 of its peak
            for (s1, e1), (s2, e2) in arms:
                mu1, mu2 = (s1 + e1) / 2, (s2 + e2) / 2
                x = np.arange(max(0, int(mu1) - half), min(n, int(mu1) + half + 2))
                y = np.arange(max(0, int(mu2) - half), min(n, int(mu2) + half + 2))
                gx = np.exp(-0.5 * ((x - mu1) / sd) ** 2)
                gy = np.exp(-0.5 * ((y - mu2) / sd) ** 2)
                block = np.outer(gx, gy)
                raw[np.ix_(x, y)] += block
                raw[np.ix_(y, x)] += block.T
            raw[raw < 1e-6] = 0
    else:
        raise ValueError(f"unknown support method: {support}")

    return np.log1p(raw) if normalize == "log" else raw


def _fold_mfold(in_fasta, constraints_file, ct_output, vienna_output):
    """Fold with mfold 3.6 using UNAFold-style constraints (F i j k / P i j k).

    Replaces hybrid-ss-min, whose bioconda build (OligoArrayAux) ignores F constraints.
    mfold is run in a scratch folder: it needs its energy tables in the working folder,
    deletes them afterwards, and writes many side files. When the constraints cannot
    all be satisfied mfold writes no structure; ct_output is then removed so the
    caller treats the constraint as rejected.
    """
    import glob, tempfile

    with open(in_fasta) as f:
        fasta_text = f.read()

    with tempfile.TemporaryDirectory(prefix="mfold_") as d:
        for dat in glob.glob(os.path.join(config.MFOLD_DAT_DIR, "*.dat")):
            shutil.copy(dat, d)
        Path(d, "seq.fasta").write_text(fasta_text)
        args = [config.MFOLD_BIN, "SEQ=seq.fasta", "MAX=1"]
        if constraints_file and os.path.getsize(constraints_file) > 0:
            with open(constraints_file) as fin:
                # mfold's aux parser expects space-separated fields
                aux = [" ".join(l.split()) + "\n" for l in fin if l.strip()]
            Path(d, "seq.aux").write_text("".join(aux))
            args.append("AUX=seq.aux")
        # mfold's console output is Latin-1 (degree signs), so keep it as bytes
        subprocess.run(args, cwd=d, capture_output=True)

        result = Path(d, "seq.ct")
        if not result.exists():
            if os.path.exists(ct_output):
                os.remove(ct_output)
            if vienna_output is not None:
                raise RuntimeError(f"mfold produced no structure for {in_fasta}; see constraints in {constraints_file}")
            return
        shutil.copy(result, ct_output)

    if vienna_output is not None:
        Path(vienna_output).write_text(ct2b_gk3(Path(ct_output).read_text(), hybrid_ss_min=True))


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

from hyb2.tools.config import str2bool

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
    p.add_argument("-r", dest="fold", default="vienna", choices=["vienna", "mfold", "unafold", "cplfold", "1", "0"], help="folding backend: (1=vienna, 0=unafold are legacy aliases)")
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
    p.add_argument("--allow-pseudoknot", dest="allow_pseudoknot", type=str2bool, default=config.CPL_DEFAULTS["allow_pseudoknot"], metavar="BOOL", help="cplfold pseudoknot search; slow on long RNAs (default false)")
    p.add_argument("--energy-model", dest="energy_model", choices=["DP09", "DP03", "CC06", "CC09", "RE"], default=config.CPL_DEFAULTS["energy_model"], help="cplfold energy model")
    p.add_argument("--support", dest="support", choices=["cofold", "coverage", "gaussian"], default=config.CPL_DEFAULTS["support"], help="cplfold support matrix: cofold (default), coverage, gaussian")
    p.add_argument("--support-hyb", dest="support_hyb", default=None, metavar="WINDOW.HYB", help="chimeras for --support coverage/gaussian")

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
        allow_pseudoknot=args.allow_pseudoknot,
        energy_model=args.energy_model,
        support=args.support,
        support_hyb=args.support_hyb
    )

    return 0

if __name__ == "__main__":
    sys.exit(main())