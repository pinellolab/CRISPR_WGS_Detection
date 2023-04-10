from utils import (
    VCALLINGTOOLS,
    GENOME,
    OFFTARGETS,
    GUIDES,
    CELLTYPES,
    OUTDIR,
    BAMS,
    create_result_dirtree
)

import multiprocessing
import subprocess
import tempfile
import argparse
import sys
import os

MUTECTPY = "run_mutect.py"
STRELKAPY = "run_strelka.py"
STRELKARUNDIR = tempfile.mkdtemp()
VARSCANPY = "run_varscan.py"

def parse_commandline():
    """Parse the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    # parse input arguments using argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Run variant calling tools to detect CRISPR-Cas9 edits",
        usage="\n\tpython %(prog)s --tool <TOOL-NAME> --threads "
        "<THREADS>",
    )
    parser.add_argument(
        "--tool",
        type=str,
        metavar="TOOL-NAME",
        help="Variant calling tool. Available values: <mutect2, strelka, "
        "varscan, pindel>",
    )
    parser.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        nargs="?",
        default=1,
        help="Number of threads to use during the run (use 0 to autodetect "
        "and use all the available resources)",
    )
    args = parser.parse_args()
    # check arguments consistency
    if args.tool not in VCALLINGTOOLS:
        raise ValueError(
            "%s cannot run %s. Check the help for the available tools"
            % (__file__, args.tool)
        )
    if args.threads < 0:
        raise ValueError("Forbidden number of threads selected (%d)" % (args.threads))
    if args.threads == 0:  # autodetect
        args.threads = multiprocessing.cpu_count()
    return args

def run_commands(commands, tool):
    """Run the input commands on the specified input regions in parallel

    :param commands: commands to execute
    :type commands: List[str]
    :param tool: tool
    :type tool: str
    :raises OSError: raise on subprocess.call() failure
    """
    for cmd in commands:
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error ocuurred while running %s\n(%s)" % (tool, cmd))

def runinfo(args):
    """(PRIVATE)
    Print run information

    :param args: input arguments
    :type args: arparse.ArgumentParser
    """

    sys.stderr.write("-- RUN INFO --\n\n")
    sys.stderr.write("\tTOOL:\t%s\n" % (args.tool))
    sys.stderr.write("\tTHREADS:\t%d\n\n" % (args.threads))

def run_mutect(threads):
    """Run Mutect2 to call edits on the input regions

    :param threads: threads
    :type threads: int
    """
    commands = []
    for guide in GUIDES:
        outdir = os.path.join(OUTDIR, VCALLINGTOOLS[0])
        targets = os.path.join(OFFTARGETS, "casoffinder.%s.txt.out" % (guide.replace("Site4", "4").replace("Site3", "3")))
        for cell_type in CELLTYPES:
            cmd = str(
                "python %s --targets %s --genome %s --bam1 %s --bam2 %s "
                "--normal DNMT1Site3 --out %s --thread %d"
            )
            if cell_type == CELLTYPES[0]:  # GM12878
                bam1 = os.path.join(BAMS, "%s.cram" % (guide)) if guide == GUIDES[0] else os.path.join(BAMS, "%s.bam" % (guide))
                bam2 = os.path.join(BAMS, "DNMT1Site3.bam")
            else:  # K562
                bam1 = os.path.join(BAMS, "%s_%s.cram" % (cell_type, guide))
                bam2 = os.path.join(BAMS, "K562_DNMT1Site3.cram")
            odir = os.path.join(outdir, cell_type, guide)
            commands.append(
                cmd % (MUTECTPY, targets, GENOME, bam1, bam2, odir, threads)
            )
    run_commands(commands, VCALLINGTOOLS[0])

def run_strelka():
    """Run strelka to call edits on the input regions
    """
    commands = []
    for guide in GUIDES:
        outdir = os.path.join(OUTDIR, VCALLINGTOOLS[1])
        targets = os.path.join(OFFTARGETS, "casoffinder.%s.txt.out" % (guide.replace("Site4", "4").replace("Site3", "3")))
        for cell_type in CELLTYPES:
            cmd = (
                "python %s --targets %s --genome %s --normal-bam %s --tumor-bam "
                "%s --run-dir %s --out %s"
            )
            if cell_type == CELLTYPES[0]:  # GM12878
                tumor_bam = os.path.join(BAMS, "%s.cram" % (guide)) if guide == GUIDES[0] else os.path.join(BAMS, "%s.bam" % (guide))
                normal_bam = os.path.join(BAMS, "DNMT1Site3.bam")
            else:  # K562
                tumor_bam = os.path.join(BAMS, "%s_%s.cram" % (cell_type, guide))
                normal_bam = os.path.join(BAMS, "%s_DNMT1Site3.cram" % (cell_type))
            odir = os.path.join(outdir, cell_type, guide)
            commands.append(
                cmd % (STRELKAPY, targets, GENOME, normal_bam, tumor_bam, STRELKARUNDIR, odir)
            )
    run_commands(commands, VCALLINGTOOLS[1])

def run_varscan():
    """Run Varscan to call edits on the input regions
    """
    commands = []
    for guide in GUIDES:
        outdir = os.path.join(OUTDIR, VCALLINGTOOLS[3])
        targets = os.path.join(OFFTARGETS, "casoffinder.%s.txt.out" % (guide.replace("Site4", "4").replace("Site3", "3")))
        for cell_type in CELLTYPES:
            cmd = (
                "python %s --targets %s --genome %s --normal-bam %s --tumor-bam "
                "%s --out %s"
            )
            if cell_type == CELLTYPES[0]:  # GM12878
                tumor_bam = os.path.join(BAMS, "%s.cram" % (guide)) if guide == GUIDES[0] else os.path.join(BAMS, "%s.bam" % (guide))
                normal_bam = os.path.join(BAMS, "DNMT1Site3.bam")
            else:  # K562
                tumor_bam = os.path.join(BAMS, "%s_%s.cram" % (cell_type, guide))
                normal_bam = os.path.join(BAMS, "K562_DNMT1Site3.cram")
            odir = os.path.join(outdir, cell_type, guide)
            commands.append(
                cmd % (VARSCANPY, targets, GENOME, normal_bam, tumor_bam, odir)
            )
    run_commands(commands, VCALLINGTOOLS[3])


def main():
    args = parse_commandline()
    runinfo(args)  # print run info to stderr
    if args.tool == VCALLINGTOOLS[0]:  # mutect2
        create_result_dirtree(args.tool)
        run_mutect(args.threads)
    elif args.tool == VCALLINGTOOLS[1]:  # strelka
        create_result_dirtree(args.tool)
        run_strelka()
    elif args.tool == VCALLINGTOOLS[3]:  # varscan
        create_result_dirtree(args.tool)
        run_varscan()

if __name__ == "__main__":
    main()

