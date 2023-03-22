"""
"""
from command_runners import run_commands
from targets_file_parser import parse_targets_coordinates

from glob import glob
from tqdm import tqdm

import argparse
import sys
import os

MUTECT2 = "gatk Mutect2"  # Mutect2

def parse_commandline():
    """The function parses the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified input regions",
        usage="\n\tpython %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--bam1 <BAM> --bam2 <BAM> --normal <NORMAL-SAMPLE> --out <OUTDIR> "
        "--thread <THREADS>"
    )
    parser.add_argument(
        "--targets", type=str, metavar="TARGETS-FILE", help="Targets coordinates file"
    )
    parser.add_argument("--genome", type=str, metavar="GENOME", help="Reference genome")
    parser.add_argument("--bam1", type=str, metavar="BAM", help="BAM file")
    parser.add_argument("--bam2", type=str, metavar="BAM", help="BAM file")
    parser.add_argument(
        "--normal", type=str, metavar="NORMAL-SAMPLE", help="Normal sample"
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    parser.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        nargs="?",
        default=1,
        help="Number of threads to use during the run",
    )
    return parser.parse_args()  # parse command line


def run_mutect(genome, bam1, bam2, normal, targets, threads, out):
    commands = [
        "%s -R %s -I %s -I %s -normal %s -L %s -O %s > /dev/null" 
        % (MUTECT2, genome, bam1, bam2, normal, target, os.path.join(out, "%s.vcf" % (target)))
        for target in targets
    ]
    run_commands(commands, threads)  # run mutect in parallel

def main():
    args = parse_commandline()
    targets = parse_targets_coordinates(args.targets)
    run_mutect(args.genome, args.bam1, args.bam2, args.normal, targets, args.threads, args.out)
    assert len(targets) == len(glob(os.path.join(args.out, "*.vcf")))

if __name__ == "__main__":
    main()


