"""Python port of the HYB2 pipeline (in-progress migration from bin/*.pl, *.sh, *.awk).

Stages are ported incrementally; see fixtures/ for the golden-output
regression harness used to validate each ported stage against the legacy
Bash/Perl/AWK implementation in bin/.
"""

__version__ = "0.1.0"
