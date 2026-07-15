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

class Hybrid2:
    def __init__(self) -> None:
        self._seq_ID = None
        self._seq = None
        self._dG = None
        self._bit1_nm = None
        self._bit2_nm = None
        self._bit1_st_in_rd = None
        self._bit1_end_in_rd = None
        self._bit2_st_in_rd = None
        self._bit2_end_in_rd = None
        self._bit1_eval = None
        self._bit2_eval = None
        self._line = None
        self._experiment = None
        self._sorted_bit_nm = None
        self._bit1_st = None
        self._bit1_end = None
        self._bit2_st = None
        self._bit2_end = None
        self._count = None
        self._found_overlap = None
        self._twoway_overlap = None
        self._seq_ID_list = []

    def initialize_hyb(self, line: str, exp=None, *_ignored) -> bool:
        line = line.rstrip("\n")
        f = line.split("\t")

        self._seq_ID = f[0] if len(f) > 0 else None
        self._seq = f[1] if len(f) > 1 else None
        self._dG = f[2] if len(f) > 2 else None
        self._bit1_nm = f[3] if len(f) > 3 else None
        self._bit1_st_in_rd = f[4] if len(f) > 4 else None
        self._bit1_end_in_rd = f[5] if len(f) > 5 else None
        self._bit1_st = int(f[6]) if len(f) > 6 else None
        self._bit1_end = int(f[7]) if len(f) > 7 else None
        self._bit1_eval = f[8] if len(f) > 8 else None
        self._bit2_nm = f[9] if len(f) > 9 else None
        self._bit2_st_in_rd = f[10] if len(f) > 10 else None
        self._bit2_end_in_rd = f[11] if len(f) > 11 else None
        self._bit2_st = int(f[12]) if len(f) > 12 else None
        self._bit2_end = int(f[13]) if len(f) > 13 else None
        self._bit2_eval = f[14] if len(f) > 14 else None
        last_column = f[15] if len(f) > 15 else None

        self._line = line
        self._experiment = exp


        self._found_overlap = 1
        self._twoway_overlap = 0
        self._seq_ID_list = [self._seq_ID]

        self._sorted_bit_nm = (
            f"{self._bit1_nm}\t{self._bit2_nm}"
            if self._bit1_nm < self._bit2_nm
            else f"{self._bit2_nm}\t{self._bit1_nm}"
        )

 
        if last_column is not None and re.fullmatch(r"[0-9]+", last_column):
            self._count = int(last_column)
        elif last_column is not None and re.match(r"[a-zA-Z_0-9]+=.*", last_column):
            m = re.search(r"(^|;)count_total=([0-9]+)", last_column)
            if m:
                self._count = int(m.group(2))
            m = re.search(r"(^|;)two_way_merged=([0-9]+)", last_column)
            if m:
                self._twoway_overlap = int(m.group(2))
            m = re.search(r"(^|;)seq_IDs_in_cluster=([^;]+)", last_column)
            if m:
                self._seq_ID_list = m.group(2).split(";")
        else:
            self._count = 1

        return True
    
    def check_hit_names(self, nm1, nm2) -> bool:
        pass
    
    def count(self) -> int:
        pass

    def found_overlap(self, fnd=None) -> int:
        pass

    def get_sorted_bit_names(self) -> str:
        pass

    def merge_with(self, other) -> None:
        pass

    def print_hyb(self) -> None:
        pass

    def print_hyb_15_columns(self) -> None:
        pass

    def reverse_bit_order(self) -> bool:
        pass

    def seq_ID_list(self) -> list[str]:
        pass

    def twoway_overlap(self, ovlp=None) -> int:
        pass

    def touches(self, other) -> bool:
        pass
