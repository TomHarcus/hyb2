""" Port of hyb2

The main orchestrator that ties in the whole pipeline

"""

def print_help():
    print("No options specified!")
    print("Usage:")
    print("To run hyb2, make sure you've activated the conda environment:")
    print("\tconda activate hyb2")
    print("hyb2-py -i <input.fastq/sam -d <fasta_file> -o <output_ID> -a <gene_ID_1 -x <start_coord_1> -y <start_coord_2> -l <length_of_fragments>")
    print("Options:")
    print("\t --verbose (-V) add if you want verbose mode on")
    print("\t --config input args in one yaml file")
    print("\t--input (-i) input fastq/sam file/hyb")
    print("\t--reference (-d) fasta file used for mapping")
    print("\t--output-id (-o) output ID")
    print("\t--blast-threshold (-v) BLAST threshold (default=0.1)")
    print("\t--max-overlap (-m) maximum overlap (default=4)")
    print("\t--max-hits (-h) maximum hits per sequence (default=10)")
    print("\t--gene-1 (-a) gene ID of interest")
    print("\t--gene-2 (-b) second gene ID of interest (if different)")
    print("\t--heatmap-quantile (-q) upper limit for heatmap chimeric count (default=0.95)")
    print("\t--x-start (-x) start coordinate of 1st strand/gene for zoomed-in contact and folding")
    print("\t--y-start (-y) start coordinate of 2nd strand/gene for zoomed-in contact and folding")
    print("\t--length (-l) length of fragments for zoomed-in contact and folding")
    print("\t--varna-jar (-j) directory of VARNAcmd.jar (default set when installing hyb2)")
    print("\t--calc-energy (-e) calculate folding energy: 1 to calculate, 0 to skip and save on runtime (default=0)")
    print("\t--fold-backend (-r) folding algorithm: 'cplfold' for CPLfold, 'unafold' or '0' for UNAfold, 'vienna' or '1' for ViennaRNA (default='cplfold')")
    print("\t--interactive (-0) interactive mode for VARNA pop-up: 0 to disable, 1 to activate (default=0)")
    print("\t--reproducible deterministic run-to-run output (adds bowtie2 --reorder)")

    # folding parameters
    print("")
    print("Optional configs for cplfold")
    print("\t--alpha cplfold bonus weight (default=0.5)")
    print("\t--beta cplfold bonus weight (default=0.0)")
    print("\t--normalize cplfold bonus normalization either 'raw' or 'log' (default='log')")
    print("\t--beam-size cplfold beam size (default=100)")
    print("\t--energy-delta cplfold energy delta (default=5.0)")
    print("\t--max-phase1 cplfold max phase 1 (default=10)")
    print("\t--max-phase2 cplfold max phase 2 (default=5)")
    print("\t--energy-model cplfold energy model: 'DP09', 'DP03', 'CC06', 'CC09', 'RE' (default='DP09')")

    print("")

    print("To only plot contact density map after generating hyb output:")
    print("hyb2 -i <output.hyb> -a <gene_ID_1>")
    print("Or:")
    print("hyb2-coverage -i <output.hyb> -a <gene_ID_1> -w <x_coord_1> -x <x_coord_2> -y <y_coord_1> -z <y_coord_2>")
    print("\t --verbose (-V) add if you want verbose mode on")
    print("\t --config input args in one yaml file")
    print("\t--input (-i) output.hyb from running hyb2")
    print("\t--gene-1 (-a) gene ID of interest")
    print("\t--gene-2 (-b) second gene ID of interest (if different)")
    print("\t--limit (-q) upper limit for heatmap chimeric count (default=0.95)")
    print("\t--x1 (-w) start coordinate of x1 for zoomed-in contact")
    print("\t--x2 (-x) start coordinate of x2 for zoomed-in contact")
    print("\t--y1 (-y) start coordinate of y1 for zoomed-in contact (or 2nd gene coord_1 if -b)")
    print("\t--y2 (-z) start coordinate of y2 for zoomed-in contact (or 2nd gene coord_2 if -b)")

    print("")

    print("To only fold RNA structures:")
    print("hyb2-fold -i <output.hyb> -d <fasta_file> -a <gene_ID_1> -x <start_coord_1> -y <start_coord_2> -l <length_of_fragments>")
    print("\t --verbose (-V) add if you want verbose mode on")
    print("\t --config input args in one yaml file")
    print("\t--input (-i) output.hyb from running hyb2")
    print("\t--reference (-d) fasta file used for mapping")
    print("\t--gene-1 (-a) gene ID of interest")
    print("\t--gene-2 (-b) second gene ID of interest (if different)")
    print("\t--x-start (-x) start coordinate of 1st strand/gene for folding")
    print("\t--y-start (-y) start coordinate of 2nd strand/gene for folding")
    print("\t--length (-l) length of fragment")
    print("\t--varna-jar (-j) directory of VARNAcmd.jar (default set when installing hyb2)")
    print("\t--fold-backend (-r) folding algorithm: 'cplfold' for CPLfold, 'unafold' or '0' for UNAfold, 'vienna' or '1' for ViennaRNA (default='cplfold')")
    print("\t--interactive (-0) interactive mode for VARNA pop-up: 0 to disable, 1 to activate (default=0)")

    # folding parameters
    print("")
    print("Optional configs for cplfold")
    print("\t--alpha cplfold bonus weight (default=0.5)")
    print("\t--beta cplfold bonus weight (default=0.0)")
    print("\t--normalize cplfold bonus normalization either 'raw' or 'log' (default='log')")
    print("\t--beam-size cplfold beam size (default=100)")
    print("\t--energy-delta cplfold energy delta (default=5.0)")
    print("\t--max-phase1 cplfold max phase 1 (default=10)")
    print("\t--max-phase2 cplfold max phase 2 (default=5)")
    print("\t--energy-model cplfold energy model: 'DP09', 'DP03', 'CC06', 'CC09', 'RE' (default='DP09')")

    print("")

    print("To start GUI:")
    print("hyb2_app -i <output.hyb> -d <fasta_file> -a <gene_id_1>")
    print("\t-i output.hyb from running hyb2")
    print("\t-d fasta file used for mapping")
    print("\t-a gene ID of interest")
    print("\t-b second gene ID of interest (if different)")
    print("\t-j directory of VARNAcmd.jar (default set when installing hyb2)")

    print("")

    print("To compare between datasets after generating hyb2 files:")
    print("hyb2-compare -i <input.table> -o <output_id> -a <gene_id> -d <fasta_file>")
    print("Run hyb2-compare without options for more help details regarding hyb2-compare")

    print("")

    print("Any other queries, email me at laujianyou@live.com")

    return 0

from pathlib import Path

import sys

from hyb2.mapping.fasta_hyb2_formatting import fasta_hyb2_formatting
from hyb2.mapping.bowtie2_map import bowtie2_map
from hyb2.chimera import sam_composition
from hyb2.coverage.hyb2_coverage import hyb2_coverage
from hyb2.viewpoint.plot_viewpoint import plot_viewpoint
from hyb2.folding import hyb2_fold
from hyb2.tools.config import CPL_DEFAULTS

import logging

log = logging.getLogger(__name__)

from hyb2.tools import ui
steps = ui.Steps()

def run(in_file, db, out, *, hmax, blast_threshold, max_overlap,
        gene_1, gene_2, limit, x_coord, y_coord, length, varna, fold,
        interactive=False, reproducible=False, alpha=CPL_DEFAULTS["alpha"], beta=CPL_DEFAULTS["beta"],
        normalize=CPL_DEFAULTS["normalize"], beam_size=CPL_DEFAULTS["beam_size"],
        energy_delta=CPL_DEFAULTS["energy_delta"], max_phase1=CPL_DEFAULTS["max_phase1"],
        max_phase2=CPL_DEFAULTS["max_phase2"], energy_model=CPL_DEFAULTS["energy_model"]):

    
    steps.start("Processing input")

    """
    Step 1: 
    process SAM files to Hyb format
    """

    db_text = Path(db).read_text()

    if "|" in db_text:
        formatted = fasta_hyb2_formatting(db_text)
        db = db.replace(".fasta", ".hyb.fasta", 1)
        Path(db).write_text(formatted)

    # split to determine what file type
    ext = in_file.split(".")
    mappable = ext[-1] in ("fasta", "fastq") or ext[-2:] == ["fastq", "gz"]
    hyb_path = f"{out}.hyb"

    if ext[-1] == "hyb":
        out = in_file.replace(".hyb", "", 1)
        hyb_path = in_file
        log.info("Hyb format as input detected")

    elif not Path(hyb_path).is_file():
        if mappable:
            bowtie2_map(in_file, db, out, reproducible=reproducible)
            log.info("SAM file generated")
            sam = f"{out}.sam"

        else:
            sam = in_file

        sam_composition.run(sam, out=out, hmax=hmax, blast_threshold=blast_threshold, max_overlap=max_overlap)
        log.info("Hyb file generated")

    else:
        log.debug("Hyb file exists. Next step.")

    """
    Step 2:
    Plot contact density map of selected genes
    """

    if gene_1:
        steps.start("Contact density maps")
        hyb2_coverage(hyb_path, gene_1, gene_2=None, limit=limit, x1=None, x2=None, y1=None, y2=None)

    if gene_2:
        hyb2_coverage(hyb_path, gene_1, gene_2, limit=limit, x1=None, x2=None, y1=None, y2=None)
    else:
        log.debug("To plot contact density map of second gene, add to command -b <gene_ID_2>")


    if not gene_2 and not y_coord and length:

        x_end = x_coord + length - 1
        hyb2_coverage(hyb_path, gene_1, gene_2=None, limit=limit, x1=x_coord, x2=x_end, y1=x_coord, y2=x_end)

    elif not gene_2 and y_coord and length:

        x_end = x_coord + length - 1
        y_end = y_coord + length - 1
        hyb2_coverage(hyb_path, gene_1, gene_2=None, limit=limit, x1=x_coord, x2=x_end, y1=y_coord, y2=y_end)

    elif gene_2 and y_coord and length:

        x_end = x_coord + length - 1
        y_end = y_coord + length - 1
        hyb2_coverage(hyb_path, gene_1, gene_2, limit=limit, x1=x_coord, x2=x_end, y1=y_coord, y2=y_end)

    else:
        log.info("To plot zoomed-in contact density map, add to command: -x <x_coord> -y <y_coord> -l <length>")

    # Plot viewpoint graphs of selected genes

    if not gene_2:
        steps.start("Viewpoint graph")
        plot_viewpoint(hyb_path, db, gene_1=gene_1, gene_2=None)

    else:
        steps.start("Viewpoint graph")
        plot_viewpoint(hyb_path, db, gene_1=gene_1, gene_2=gene_2)

    """
    Step 3:
    generate RNA secondary structure of selected strands
    """

    if x_coord and length and not gene_2 and not y_coord:
        steps.start(f"Folding {gene_1} (this can take a while)")
        hyb2_fold.run(hyb_path, GENE_1=gene_1, GENE_2=None, FASTA_1=db,
              x_coord=x_coord, y_coord=None, length=length,
              VARNA=varna, interactive=interactive, FOLD=fold,
              alpha=alpha, beta=beta, normalize=normalize, beam_size=beam_size,
              energy_delta=energy_delta, max_phase1=max_phase1,
              max_phase2=max_phase2, energy_model=energy_model)

    elif x_coord and length and not gene_2 and y_coord:
        steps.start(f"Folding {gene_1} (this can take a while)")
        hyb2_fold.run(hyb_path, GENE_1=gene_1, GENE_2=None, FASTA_1=db,
              x_coord=x_coord, y_coord=y_coord, length=length,
              VARNA=varna, interactive=interactive, FOLD=fold,
              alpha=alpha, beta=beta, normalize=normalize, beam_size=beam_size,
              energy_delta=energy_delta, max_phase1=max_phase1,
              max_phase2=max_phase2, energy_model=energy_model)

    elif x_coord and length and gene_2 and y_coord:
        steps.start(f"Folding {gene_1} with {gene_2} (this can take a while)")
        hyb2_fold.run(hyb_path, GENE_1=gene_1, GENE_2=gene_2, FASTA_1=db,
              x_coord=x_coord, y_coord=y_coord, length=length,
              VARNA=varna, interactive=interactive, FOLD=fold,
              alpha=alpha, beta=beta, normalize=normalize, beam_size=beam_size,
              energy_delta=energy_delta, max_phase1=max_phase1,
              max_phase2=max_phase2, energy_model=energy_model)

    else:
        log.info("No options specified to generate secondary structure")

    log.info("Analysis completed")
    log.info("To compare between different datasets and plot differential coverage map, similarity heatmap, and differential structures, use:")
    log.info("hyb2_compare -i <input.table> -a <gene_ID> -d <fasta>")
    log.info("For more details, run: hyb2_compare")