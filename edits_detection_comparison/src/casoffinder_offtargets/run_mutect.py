"""Run GATK Mutect2 on the set of off-target regions identified by CasOffinder. 
The regions are specified in the input file.
"""

from command_runners import run_commands
from targets_file_parser import parse_targets_coordinates

from glob import glob

import argparse
import os

MUTECT2 = "gatk Mutect2"  # Mutect2


def parse_commandline():
    """Parse the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified input regions",
        usage="\n\tpython %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--bam1 <BAM> --bam2 <BAM> --normal <NORMAL-SAMPLE> --out <OUTDIR> "
        "--thread <THREADS>",
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
        "--offtarget-upstream",
        action="store_true",
        default=False,
        dest="offtarget_upstream",
        help="Detect edits on the 20bp upstream the target site",
    )
    parser.add_argument(
        "--offtarget-downstream",
        action="store_true",
        default=False,
        dest="offtarget_downstream",
        help="Detect edits on the 20bp downstream the target site",
    )
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
    """Wrapper function to call edits on the input regions using GATK Mutect2

    :param genome: reference genome
    :type genome: str
    :param bam1: BAM
    :type bam1: str
    :param bam2: BAM
    :type bam2: str
    :param normal: normal sample
    :type normal: str
    :param targets: target regions
    :type targets: List[str]
    :param threads: threads
    :type threads: int
    :param out: out
    :type out: _type_
    """
    commands = [
        "%s -R %s -I %s -I %s -normal %s -L %s -O %s > /dev/null"
        % (
            MUTECT2,
            genome,
            bam1,
            bam2,
            normal,
            target,
            os.path.join(out, "%s.vcf" % (target)),
        )
        for target in targets
    ]
    run_commands(commands, threads)  # run mutect in parallel


def main():
    args = parse_commandline()
    targets = parse_targets_coordinates(
        args.targets, args.offtarget_upstream, args.offtarget_downstream
    )
    run_mutect(
        args.genome, args.bam1, args.bam2, args.normal, targets, args.threads, args.out
    )
    assert len(targets) == len(glob(os.path.join(args.out, "*.vcf")))


if __name__ == "__main__":
    main()
