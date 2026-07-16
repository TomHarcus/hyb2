# HYB2 Python migration — file-by-file order

Topological migration order for `bin/` → `src/hyb2/`. Every file's dependencies
appear before it, so when you reach a file everything it needs is already ported.
Within a tier the order is flexible **except** where `→` marks a hard dependency.

Parity oracle: `fixtures/sam_composition_run/` (regenerate with
`scripts/generate_sam_composition_baseline.sh`). Diff each ported stage against it.

**Testing approach (per Grzegorz): real-data golden diffs.** Each stage runs on
the captured real input and diffs against the golden output in
`fixtures/sam_composition_run/`; plus one end-to-end pipeline parity test. No
hand-written synthetic inputs, no re-running legacy scripts at test time. Tier 1
test files exist as `xfail` stubs (`tests/test_<stage>.py`) that flip to `xpass`
when a stage is implemented — then remove the marker. Layout: stages in
`src/hyb2/stages/`, the `Hybrid` class in `src/hyb2/models.py`, the orchestrator
in `src/hyb2/pipelines/sam_composition.py`.

Legend: `[ ]` todo · `[x]` done · `[~]` already Python (adopt/verify only) ·
`(!)` needs confirmation from Grzegorz.

---

## Excluded — do NOT migrate

**Stays as-is (call from Python, don't rewrite):**
- R scripts (10): `cdm_2genes.R`, `cdm_indiv_zoom.R`, `contact_density_map.R`,
  `contact_density_map_indiv.R`, `contact_density_map_zoom.R`, `DESeq_run.R`,
  `differential_coverage_map.R`, `hyb2_shiny.R`, `similarity_heatmap.R`,
  `viewpoint_graph.R`
- VARNA Java library (`VARNA/`, jar) — external tool
- Packaging/binaries: `faToTwoBit`, `hyb2_install*`, `*.yml`, `renv*`

**Already Python (adopt into package, no translation):**
- [~] `sam2blast_3` (SAM → blast; on the spine)
- [~] `DESeq_interaction_split_select.py`
- [x] `mtophits_blast` → `src/hyb2/stages/mtophits_blast.py` (verify vs oracle)
- [ ] `hyb2_composition_pies.py` — **missing locally, get from Grzegorz**

**Already ported by Tom (keep; not on the spine except mtophits):**
- [x] `fasta2tab.awk` · `solexa2fasta.awk` · `filter_homopolymers.awk` ·
  `txt2hyb.awk` · `blast_stats_2` · `hybrid_stats_2`

---

## TIER 1 — the spine (`hyb2_sam_composition.sh`) — DO NOW

- [x] 1. `create_reference_file.pl` — leaf
- [x] 2. `histogram.pl` — leaf (used twice; also needed in Tier 2)
- [x] 3. `collapse_blast_2.sh` — leaf (port, or keep shelling to Grzegorz's memory-safe version)
- [x] 4. `Hybrid_long.pm` → 5 — port only the ~7 methods `remove_duplicate` uses, as a small `Hybrid` class
- [x] 5. `remove_duplicate_hybrids_hOH5_2.pl` — needs #4
- [x] 6. `get_mtop_hybrids.pl` — leaf, self-contained; the core chimera caller (biggest single unit)
- [x] 7. `hyb2_sam_composition.sh` — orchestrator; needs 1–6 + the two already-Python.
       Fold `_antisense` (MODE=1) and `sense_antisense` variants in as flags.

**Milestone: validated Python chimera-calling core.**

## TIER 2 — folding + CoupleFold hook (`hyb2_fold`)

- [x] 8. `ct2bps_2.awk` — leaf
- [x] 9. `hyb2constraints.pl` — leaf
- [x] 10. `hyb2fasta_bits_allRNAs.awk` — leaf
- [ ] 11. `add_dG_hyb_2.pl` — leaf (free-energy annotation)
- [ ] 12. `make_nicer_vienna_hOH5.awk` — leaf (vienna formatting)
- [ ] 13. `Ct2B_GK_3.pl` — leaf
- [ ] 14. `make_VARNA_scores_2.sh` — leaf
- [ ] 15. `svg_mod_coord.sh` — leaf
- [ ] 16. `VARNAcmd.sh` — thin wrapper to the VARNA jar (keep/port wrapper, not VARNA)
- [x] 17. `Hybrid_long_2.pm` → 18 — v2 module, distinct from #4
- [x] 18. `combine_hyb_merge_touching.pl` — needs #17
- [ ] 19. `bp2hyb.sh` — needs #18
- [ ] 20. `bp_score.sh` — needs #14
- [ ] 21. `comradesFold2` — needs #13  ← **CoupleFold folding-backend switch point**
- [ ] 22. `comradesScore` — needs #20, #21
- [ ] 23. `comradesMakeConstraints_2` — needs #8, #9, #10, `histogram.pl`, `fasta2tab`, #19
       ← **CoupleFold support-matrix point (Phase 3)**
- [ ] 24. `plot_VARNA` — needs #15, #16
- [ ] 25. `hyb2_fold` — orchestrator; needs all above

## TIER 3 — remaining README commands + glue

*Mapping sub-branch:*
- [ ] 26. `make_comp_fasta.pl` — use Grzegorz's newer 98-line version
- [ ] 27. `make_hyb_db_2` — needs `fasta2tab` (done)
- [ ] 28. `bowtie2_fastq2sam`, `bowtie2_fastq.gz2sam`, `bowtie2_fasta2sam` — need #26, #27, `solexa2fasta` (done)

*Coverage / plotting wrappers (R stays):*
- [ ] 29. `plot_hybrids_3.awk`
- [ ] 30. `hyb2_coverage` — needs #29
- [ ] 31. `plot_cdm`, `plot_all_cdm` — thin R wrappers
- [ ] 32. `hyb2blast.awk`
- [ ] 33. `blast2gplot.pl`
- [ ] 34. `plot_viewpoint` — needs #32, #33
- [ ] 35. `make_hybrid_annotation_table.pl`
- [ ] 36. `plot_differential_map`, `plot_similarity_map` — thin R/py wrappers
- [ ] 37. `fasta_hyb2_formatting`, `sam_filt` — leaves

*Top orchestrators — last:*
- [ ] 38. `hyb2_app` — needs `hyb2_coverage`
- [ ] 39. `hyb2_compare` — needs `hyb2`, `hyb2_app`, `hyb2_coverage`, `hyb2_fold`, #35
- [ ] 40. `hyb2` — the mega-orchestrator, **dead last**

## PRUNE — confirm with Grzegorz, then delete (do NOT port)

- [ ] (!) `sam2hyb` + its exclusive deps: `getseqs_sam`, `collapse_hyb_2.sh`,
       `remove_duplicate_hits_blast.pl`, `rm_blast_dup`, `hyb_chim_types.awk`,
       `hyb_mapped_position_overlaps_3.awk`, `normalize_table_columns.awk`,
       `chim_types` (data file). Only reachable via `sam2hyb`. **If
       `hyb2_sam_composition.sh` replaces `sam2hyb`, this whole cluster is dead
       weight** — the bulk of the "redundant files" Grzegorz mentioned.
- [ ] (!) `combine_hyb_merge` (old) — superseded by `combine_hyb_merge_touching.pl` (#18)
- [ ] (!) `make_bubble_array.sh` — orphan, not referenced by anything

---

## Open decision that resizes this list

**Does `hyb2_sam_composition.sh` replace `sam2hyb`, or must `sam2hyb` also be
ported?** If it replaces it, the entire PRUNE cluster is deletion, not migration.
