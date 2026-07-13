"""Shared data structures for the HYB2 port.

NOTE: this is a deliberate exception to the "one module per bin/ script" rule
that governs stages/ -- Hybrid_long.pm is a data class used across several
stages (remove_duplicate_hybrids now, combine_hyb_merge later), not a pipeline
transform, so it lives here rather than in stages/.

`Hybrid` ports bin/Hybrid_long.pm. Only the methods actually exercised by
remove_duplicate_hybrids_hOH5_2.pl are ported (the .pm has 28; these ~7 are
the used subset):
    new / initialize_hyb / seq_ID / get_bit_names / match_bit_name /
    sum_e_values / line
When Tier 2 needs Hybrid_long_2.pm (used by combine_hyb_merge_touching.pl),
add a second class here.
"""
import re

class Hybrid:
    """A single .hyb record. Port of bin/Hybrid_long.pm."""

    def __init__(self) -> None:
        self._seq_ID = None
        self._bit1_nm = None
        self._bit2_nm = None
        self._bit1_eval = None
        self._bit2_eval = None
        self._line = None
        self._experiment = None

    def initialize_hyb(self, line: str, exp=None, *_ignored) -> bool:
        """Parse one .hyb line into this object; return None to skip (as the
        Perl `initialize_hyb ... or next` idiom does)."""
        line = line.rstrip("\n")
        f = line.split("\t")

        self._seq_ID = f[0] if len(f) > 0 else None
        self._bit1_nm = f[3] if len(f) > 3 else None
        self._bit1_eval = f[8] if len(f) > 8 else None
        self._bit2_nm = f[9] if len(f) > 9 else None
        self._bit2_eval = f[14] if len(f) > 14 else None

        self._line = line
        self._experiment = exp 

        return True


    def seq_ID(self) -> str:
        return self._seq_ID

    def get_bit_names(self) -> tuple[str, str]:
        return (self._bit1_nm, self._bit2_nm)

    def match_bit_name(self, nm) -> str | None:
        """Return the bit name matching regex nm (e.g: _mRNA), else None"""
        if re.search(nm, self._bit1_nm):
            return self._bit1_nm
        elif re.search(nm, self._bit2_nm):
            return self._bit2_nm
        return None

    def sum_e_values(self) -> float:
        return float(self._bit1_eval) + float(self._bit2_eval)

    def line(self) -> str:
        """The original record text, for output."""
        return self._line
