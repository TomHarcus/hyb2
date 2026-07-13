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

from __future__ import annotations


class Hybrid:
    """A single .hyb record. Port of bin/Hybrid_long.pm."""

    def __init__(self) -> None:
        raise NotImplementedError

    def initialize_hyb(self, line: str, *args) -> "Hybrid | None":
        """Parse one .hyb line into this object; return None to skip (as the
        Perl `initialize_hyb ... or next` idiom does)."""
        raise NotImplementedError

    def seq_ID(self) -> str:
        raise NotImplementedError

    def get_bit_names(self) -> tuple[str, str]:
        raise NotImplementedError

    def match_bit_name(self, suffix: str) -> str | None:
        """Return a bit name matching `suffix` (e.g. '_mRNA', '_microRNA'), else None."""
        raise NotImplementedError

    def sum_e_values(self) -> float:
        raise NotImplementedError

    def line(self) -> str:
        """The original record text, for output."""
        raise NotImplementedError
