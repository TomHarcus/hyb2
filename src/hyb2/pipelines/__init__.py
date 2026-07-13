"""Pipeline orchestrators -- drivers that chain stages/ together.

Deliberate exception to the "one module per bin/ script" rule of stages/: a
driver wires stages and owns CLI flags / run modes, so it is not a transform.
Each module here ports one top-level bin/ orchestrator (sam_composition.py <-
hyb2_sam_composition.sh; later hyb2_fold, hyb2, ...).
"""
