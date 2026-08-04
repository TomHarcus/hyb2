""" Port of hyb2

The main orchestrator that ties in the whole pipeline

"""

def print_help():
    print("No options specified!")
    print("Usage:")
    print("To run hyb2, make sure you've activated the conda environment:")
    print("\tconda activate hyb2")
    print("hyb2 -i <input.fastq/sam -d <fasta_file> -o <output_ID> -a <gene_ID_1 -x <start_coord_1> -y <start_coord_2> -l <length_of_fragments>")
    print("Options:")
    print("\t-i input fastq/sam file")
    print("\t-d fasta file used for mapping")
    print("\t-o output ID")
    print("\t-v BLAST threshold (default=0.1)")
    print("\t-m maximum overlap (default=4)")
    print("\t-h maximum hits per sequence (default=10)")
    print("\t-a gene ID of interest")
    print("\t-b second gene ID of interest (if different)")
    print("\t-q upper limit for heatmap chimeric count (default=0.95)")
    print("\t-x start coordinate of 1st strand/gene for zoomed-in contact and folding")
    print("\t-y start coordinate of 2nd strand/gene for zoomed-in contact and folding")
    print("\t-l length of fragments for zoomed-in contact and folding")
    print("\t-j directory of VARNAcmd.jar (default set when installing hyb2)")
    print("\t-e calculate folding energy: 1 to calculate, 0 to skip and save on runtime (default=0)")
    print("\t-r folding algorithm: 'cplfold' for CPLfold, 'unafold' or '0' for UNAfold, 'vienna' or '1' for ViennaRNA (default='cplfold')")

    print("")

    print("To only plot contact density map after generating hyb output:")
    print("hyb2 -i <output.hyb> -a <gene_ID_1>")
    print("Or:")
    print("hyb2_coverage -i <output.hyb> -a <gene_ID_1> -w <x_coord_1> -x <x_coord_2> -y <y_coord_1> -z <y_coord_2>")
    print("\t-i output.hyb from running hyb2")
    print("\t-a gene ID of interest")
    print("\t-b second gene ID of interest (if different)")
    print("\t-q upper limit for heatmap chimeric count (default=0.95)")
    print("\t-w start coordinate of x1 for zoomed-in contact")
    print("\t-x start coordinate of x2 for zoomed-in contact")
    print("\t-y start coordinate of y1 for zoomed-in contact (or 2nd gene coord_1 if -b)")
    print("\t-z start coordinate of y2 for zoomed-in contact (or 2nd gene coord_2 if -b)")

    print("")

    print("To only fold RNA structures:")
    print("hyb2_fold -i <output.hyb> -d <fasta_file> -a <gene_ID_1> -x <start_coord_1> -y <start_coord_2> -l <length_of_fragments>")
    print("\t-i output.hyb from running hyb2")
    print("\t-d fasta file used for mapping")
    print("\t-a gene ID of interest")
    print("\t-b second gene ID of interest (if different)")
    print("\t-x start coordinate of 1st strand/gene for folding")
    print("\t-y start coordinate of 2nd strand/gene for folding")
    print("\t-j directory of VARNAcmd.jar (default set when installing hyb2)")
    print("\t-r folding algorithm: 'cplfold' for CPLfold, 'unafold' or '0' for UNAfold, 'vienna' or '1' for ViennaRNA (default='cplfold')")
    print("\t-0 interactive mode for VARNA pop-up: 0 to disable, 1 to activate (default=0)")

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
    print("hyb2_compare -i <input.table> -o <output_id> -a <gene_id> -d <fasta_file>")
    print("Run hyb2_compare without options for more help details regarding hyb2_compare")

    print("")

    print("Any other queries, email me at laujianyou@live.com")

    return 0

from pathlib import Path

import sys

from hyb2.stages.fasta_hyb2_formatting import fasta_hyb2_formatting
from hyb2.pipelines.bowtie2_map import bowtie2_map
from hyb2.pipelines import sam_composition
from hyb2.pipelines.hyb2_coverage import hyb2_coverage
from hyb2.pipelines.plot_viewpoint import plot_viewpoint
from hyb2.pipelines import hyb2_fold



def run(in_file, db, out, *, hmax, blast_threshold, max_overlap,
        gene_1, gene_2, limit, x_coord, y_coord, length, varna, fold):

    """
    Step 1: 
    process SAM files to Hyb format
    """

    if "|" in open(db).read():
        formatted = fasta_hyb2_formatting(open(db).read())
        db = db.replace(".fasta", ".hyb.fasta", 1)
        Path(db).write_text(formatted)

    ext = in_file.split(".")
    mappable = ext[-1] in ("fasta", "fastq") or ext[-2:] == ["fastq", "gz"]
    hyb_path = f"{out}.hyb"

    if ext[-1] == "hyb":
        out = in_file.replace(".hyb", "", 1)
        hyb_path = in_file
        print("Hyb format as input detected")

    elif not Path(hyb_path).is_file():
        if mappable:
            bowtie2_map(in_file, db, out)
            print("SAM file generated")
            sam = f"{out}.sam"

        else:
            sam = in_file

        sam_composition.run(sam, out=out, hmax=hmax, blast_threshold=blast_threshold, max_overlap=max_overlap)
        print("Hyb file generated")

    else:
        print("Hyb file exists. Next step.")

    """
    Step 2:
    Plot contact density map of selected genes
    """

    if gene_1:

        hyb2_coverage(hyb_path, gene_1, gene_2=None, limit=limit, x1=None, x2=None, y1=None, y2=None)

    if gene_2:

        hyb2_coverage(hyb_path, gene_1, gene_2, limit=limit, x1=None, x2=None, y1=None, y2=None)
    else:
        print("To plot contact density map of second gene, add to command -b <gene_ID_2>", file=sys.stderr)


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
        print("To plot zoomed-in contact density map, add to command: -x <x_coord> -y <y_coord> -l <length>", file=sys.stderr)

    # Plot viewpoint graphs of selected genes

    if not gene_2:
        plot_viewpoint(hyb_path, db, gene_1=gene_1, gene_2=None)

    else:
        plot_viewpoint(hyb_path, db, gene_1=gene_1, gene_2=gene_2)

    """
    Step 3:
    generate RNA secondary structure of selected strands
    """

    if x_coord and length and not gene_2 and not y_coord:

        hyb2_fold.run(hyb_path, GENE_1=gene_1, GENE_2=None, FASTA_1=db,
              x_coord=x_coord, y_coord=None, length=length,
              VARNA=varna, interactive=0, FOLD=fold)

    elif x_coord and length and not gene_2 and y_coord:

        hyb2_fold.run(hyb_path, GENE_1=gene_1, GENE_2=None, FASTA_1=db,
              x_coord=x_coord, y_coord=y_coord, length=length,
              VARNA=varna, interactive=0, FOLD=fold)

    elif x_coord and length and gene_2 and y_coord:

        hyb2_fold.run(hyb_path, GENE_1=gene_1, GENE_2=gene_2, FASTA_1=db,
              x_coord=x_coord, y_coord=y_coord, length=length,
              VARNA=varna, interactive=0, FOLD=fold)

    else:
        print("No options specified to generate secondary structure", file=sys.stderr)

    print("Analysis completed")
    print("To compare between different datasets and plot differential coverage map, similarity heatmap, and differential structures, use:")
    print("hyb2_compare -i <input.table> -a <gene_ID> -d <fasta>")
    print("For more details, run: hyb2_compare")