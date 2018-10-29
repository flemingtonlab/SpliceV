
import argparse
import os
import sys

def parse():
   
    parser = argparse.ArgumentParser(description='Plot a transcript')

    parser.add_argument(
        "-bsj",
        nargs='*',
        help='''Path to backsplice junction bed-formatted files (listed in the same 
        order as the input bam files)'''
    )

    parser.add_argument(
        "-stranded",
        help='''If strand-specific sequencing, indicate 'forward' if upstream reads are
        forward strand, otherwise indicate 'reverse' (True-Seq is 'reverse').'''
    )
    
    parser.add_argument(
        "-sj",
        nargs='*',
        help='''Path to canonical splice junction bed-formatted files (listed in the same
        order as the input bam files)'''
    )

    parser.add_argument(
        "-is",
        "--intron-scale",
        type=float,
        help="The factor by which intron white space should be reduced"
    )

    parser.add_argument(
        "-b",
        "--bam",
        nargs='*',
        required=True,
        help="Path to bam file(s)"
    )

    parser.add_argument(
        "-c",
        "--color",
        default="#C21807",
        type=str,
        help='''Exon color. Hex colors (i.e. \"\#4286f4\". For hex, an escape \"\\\" must precede the argument), RGB (i.e. 211,19,23)
        or names (i.e. \"red\")'''
    )

    parser.add_argument(
        "-t",
        "--transcript",
        help='Name of transcript to plot (must match "transcript_id" field of gtf file)'
    )

    parser.add_argument(
        "-g",
        "--gene",
        help='''Name of gene to plot (overrides "-t" flag). 
        Will plot the longest transcript derived from that gene'''
    )

    parser.add_argument(
        "-f",
        "--filter",
        default=0,
        type=int,
        help='Filter out sj and circles that have fewer than this number of counts.'
    )

    parser.add_argument(
        "-n",
        "--normalize",
        action='store_true',
        help='Normalize coverage between samples'
    )
    
    parser.add_argument(
        "-rc",
        "--reduce_canonical",
        type=float,
        help='Factor by which to reduce canonical curves'
    )

    parser.add_argument(
        "-rbs",
        "--reduce_backsplice",
        type=float,
        help='Factor by which to reduce backsplice curves'
    )

    parser.add_argument(
        "-ro",
        "--repress_open",
        action='store_true',
        help='Do not automatically open plot in browser'
    )

    parser.add_argument(
        "-en",
        "--exon_numbering",
        action='store_true',
        help='Label exons'
    )

    parser.add_argument(
        "-gtf",
        required=True,
        help="Path to gtf file"
    )

    parser.add_argument(
        '-rnabp',
        nargs='*',
        help="List of RNA binding proteins to plot."
    )

    parser.add_argument(
        "-fa",
        help="Path to fasta file"
    )

    parser.add_argument(
        '-format',
        default='PNG',
        choices=['SVG', 'PDF', 'PNG', 'JPG', 'TIFF'],
        help="Output image format"
    )

    parser.add_argument(
        "-alu",
        help="Path to Alu bed file"
    )

    args = parser.parse_args()

    # Check GTF, BAM, SJ, and BSJ paths.
    if not os.path.exists(args.gtf):
        sys.exit("GTF: {} was not found".format(args.gtf))

    for path in args.bam: 
        if not os.path.exists(path):
            sys.exit("BAM: {} was not found".format(path))

    if args.sj:
        for path in args.sj: 
            if not os.path.exists(path):
                sys.exit("Splice junction file {} was not found.\n If you want to obtain splice junction reads from the input bam files, don't specify a splice junction file.".format(path))
    
    if args.bsj:
        for path in args.bsj: 
            if not os.path.exists(path):
                sys.exit("Backsplice junction file {} was not found.\n If you want to obtain backsplice junction reads from the input bam files, don't specify a backsplice junction file".format(path))

    if not (args.gene or args.transcript): 
        sys.exit('Either a gene or a transcript must be specified. (ex. "-t ENST00000390665" or "-g EGFR")')    

    
    return args

