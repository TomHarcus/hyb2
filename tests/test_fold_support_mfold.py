"""Support matrices for CPLfold (cofold / coverage / gaussian) and the mfold backend."""
import shutil

import numpy as np
import pytest

from hyb2.folding.comrades_fold import support_matrix, _fold_mfold
from hyb2.tools import config


def _hyb(tmp_path, rows):
    p = tmp_path / "w.hyb"
    p.write_text("".join(
        f"{k}_1\t.\t.\tRNA\t.\t.\t{a}\t{b}\t.\tRNA\t.\t.\t{c}\t{d}\t.\n" for k, (a, b, c, d) in enumerate(rows, 1)))
    return str(p)


def test_cofold_adds_both_orientations(tmp_path):
    scores = tmp_path / "bp.txt"
    scores.write_text("3\t18\t5\n18\t3\t3\n")
    m = support_matrix("cofold", 20, 1, 20, "raw", basepair_scores=str(scores))
    assert m[2, 17] == m[17, 2] == 8


def test_coverage_counts_each_chimera_over_its_arms(tmp_path):
    hyb = _hyb(tmp_path, [(2, 4, 15, 17), (3, 4, 16, 17)])
    m = support_matrix("coverage", 20, 1, 20, "raw", support_hyb=hyb)
    assert m[1, 14] == 1 and m[2, 15] == 2 and m[15, 2] == 2 and m[0, 14] == 0
    assert m.sum() == 2 * (9 + 4)


def test_gaussian_peaks_at_arm_centres(tmp_path):
    hyb = _hyb(tmp_path, [(10, 16, 60, 66)])
    m = support_matrix("gaussian", 80, 1, 80, "raw", support_hyb=hyb)
    assert np.unravel_index(np.argmax(np.triu(m)), m.shape) == (12, 62)
    assert np.allclose(m, m.T) and m[12, 62] == pytest.approx(1.0)


@pytest.mark.skipif(shutil.which(config.MFOLD_BIN) is None, reason="mfold not installed")
def test_mfold_applies_force_and_reports_infeasible(tmp_path):
    seq = "GGGGAAAACCCCUUUUUUUUGGGGAAAACCCC"
    fa = tmp_path / "s.fasta"
    fa.write_text(f">s\n{seq}\n")
    aux, ct = tmp_path / "s.aux", tmp_path / "s.ct"

    aux.write_text("F\t1\t32\t4\n")          # tab-separated, as hyb2 writes constraints
    _fold_mfold(str(fa), str(aux), str(ct), None)
    partner = {int(l.split()[0]): int(l.split()[4]) for l in ct.read_text().splitlines()[1:]}
    assert partner[1] == 32

    aux.write_text("F\t1\t32\t4\nF\t1\t12\t4\n")   # both cannot hold
    _fold_mfold(str(fa), str(aux), str(ct), None)
    assert not ct.exists()
