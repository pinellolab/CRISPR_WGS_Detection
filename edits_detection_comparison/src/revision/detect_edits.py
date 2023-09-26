"""
"""

from mutect2 import detect_edits_mutect2
from utils import VCALLINGTOOLS, check_args

import multiprocessing
import argparse


def parse_commandline():
    """Parse command line arguments

    :return: input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Detect rare genome editing events on input target sites, "
        "using Mutect2, Strelka2, and Varscan2",
        usage="\n\tpython %(prog)s"
    )
    parser.add_argument(
        "--tool",
        type=str,
        required=True,
        metavar="TOOL-NAME",
        help="Variant calling tool <mutect2, strelka, varscan>",
    )
    parser.add_argument(
        "--targets",
        type=str,
        required=True,
        metavar="TARGETS-FILE",
        help="Target sites file",
    )
    parser.add_argument(
        "--genome",
        type=str,
        required=True,
        metavar="GENOME",
        help="Genome FASTA",   
    )
    parser.add_argument(
        "--bam1",
        type=str,
        required=True,
        metavar="BAM-FILE",
        help="Normal BAM"
    )
    parser.add_argument(
        "--bam2",
        type=str,
        required=True,
        metavar="BAM-FILE",
        help="Tumor BAM"
    )
    # TODO: check if mandatory for strelka and varscan
    parser.add_argument(
        "--normal-sample",
        type=str,
        dest="normal_sample",
        nargs="?",
        default="",
        metavar="NORMAL-SAMPLE",
        help="Normal sample",
    )
    parser.add_argument(
        "--out",
        type=str,
        required=True,
        metavar="OUTDIR",
        help="Output directory",
    )
    parser.add_argument(
        "--upstream",
        action="store_true",
        default=False,
        help="Detect genome editing events on 10 kb region upstream the 100 bp "
        "target site",
    )
    parser.add_argument(
        "--downstream",
        action="store_true",
        default=False,
        help="Detect genome editing events on 10 kb region downstream the 100 bp "
        "target site",
    )
    parser.add_argument(
        "--casoffinder",
        action="store_true",
        default=False,
        help="Targets file retrieved from CasOffinder",
    )
    parser.add_argument(
        "--circleseq",
        action="store_true",
        default=False,
        help="Targets file retrieved from CIRCLE-seq",
    )
    parser.add_argument(
        "--guideseq",
        action="store_true",
        default=False,
        help="Targets file retrieved from GUIDE-seq",
    )
    parser.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        nargs="?",
        default=1,
        help="Threads (0 to autodetect)",
    )
    args = parser.parse_args()
    check_args(args)
    return args    

def main():
    args = parse_commandline()
    threads = multiprocessing.cpu_count() if args.threads == 0 else args.threads
    if args.tool == VCALLINGTOOLS[0]:  # mutect2
        detect_edits_mutect2(
            args.targets,
            args.genome,
            args.bam1,
            args.bam2,
            args.normal_sample,
            args.out,
            args.casoffinder,
            args.circleseq,
            args.guideseq,
            args.upstream,
            args.downstream,
            threads,
        )

if __name__ == "__main__":
    main()
