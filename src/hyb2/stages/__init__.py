"""Ported pipeline stages, one module per legacy bin/ script.

Each module here should correspond to a single script currently in bin/
(e.g. stages/fasta2tab.py ports bin/fasta2tab.awk) so a stage can be
swapped from legacy to Python independently of its neighbors. Add stages
here as they're ported and validated against fixtures/.
"""
