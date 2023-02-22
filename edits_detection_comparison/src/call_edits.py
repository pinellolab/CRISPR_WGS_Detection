"""This python script calls different variant calling tools to detect 
mutations induced by CRISPR-CAs9 editing.
"""

import multiprocessing
import subprocess
import argparse
import tempfile
import sys
import os

# tool wrappers
MUTECTPY = "run_mutect.py"
STRELKAPY = "run_strelka.py"
PINDELPY = "run_pindel.py"
VARSCANPY = "run_varscan.py"
STRELKARUNDIR = tempfile.mkdtemp()  # strelka run utils
# guides, cell types and validation experiment types
GUIDES = ["EMX1", "HEKSite4", "RNF2", "VEGFASite3"]
CELLTYPES = ["GM12878", "K562"]
EXPERIMENTTYPE = ["circleseq", "guideseq"]
# variant calling tools
VCALLINGTOOLS = ["mutect2", "strelka", "pindel", "varscan"]
# data directories
BASEDIR = "../data"
GUIDESEQ = os.path.join(BASEDIR, "offtarget_detection/guideseq")
CIRCLESEQ = os.path.join(BASEDIR, "offtarget_detection/circleseq/")
BAMS = os.path.join(BASEDIR, "wgs/bams/")
GENOME = os.path.join(BASEDIR, "genome/Homo_sapiens_assembly38.fasta")
OUTDIR = os.path.join(BASEDIR, "VCFs")


def parse_commandline():
    """The function parses the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    # parse input arguments using argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Run variant calling tools to detect CRISPR-Cas9 edits",
        usage="\n\tpython %(prog)s --tool <TOOL-NAME> --type <TYPE> --threads "
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
        "--type",
        type=str,
        metavar="TYPE",
        help="Target sites validation experiment type. Available values: "
        "<guideseq, circleseq>",
    )
    parser.add_argument(
        "--offregion",
        action="store_true",
        default=False,
        help="Shift the target regions 100bp upstream and downstream (not "
        "overlapping the original target site)",
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
    if args.type != "circleseq" and args.type != "guideseq":
        raise ValueError(
            "Forbidden target sites validation experiment type. Check the help "
            "for the available values"
        )
    return args


def _create_celltype_dirtree(root):
    """(PRIVATE)
    Build the directory tree for each cell type (GM12878 and K562)

    :param root: root directory
    :type root: str
    """
    assert isinstance(root, str)
    gm12878_dir = os.path.join(root, CELLTYPES[0])
    # check cell types directories
    if not os.path.isdir(gm12878_dir):
        os.mkdir(gm12878_dir)
    onregion_dir = os.path.join(gm12878_dir, "onregion")
    if not os.path.join(onregion_dir):
        os.mkdir(onregion_dir)
    offregion_dir = os.path.join(gm12878_dir, "offregion")
    if not os.path.join(offregion_dir):
        os.mkdir(offregion_dir)
    guides_dirs = [
        os.path.join(r, g) for r in [onregion_dir, offregion_dir] for g in GUIDES
    ]
    for gd in guides_dirs:
        if not os.path.isdir(gd):
            os.mkdir(gd)
    k562_dir = os.path.join(root, CELLTYPES[1])
    if not os.path.isdir(k562_dir):
        os.mkdir(k562_dir)
    onregion_dir = os.path.join(k562_dir, "onregion")
    if not os.path.join(onregion_dir):
        os.mkdir(onregion_dir)
    offregion_dir = os.path.join(k562_dir, "offregion")
    if not os.path.join(offregion_dir):
        os.mkdir(offregion_dir)
    guides_dirs = [
        os.path.join(r, g) for r in [onregion_dir, offregion_dir] for g in GUIDES
    ]


def create_result_dirtree(tool):
    """Build the directory tree for the tool

    :param tool: variant calling tool
    :type tool: str
    """
    assert tool in VCALLINGTOOLS
    # if there is no tool tree directory, build the tree
    tool_root_dir = os.path.join(OUTDIR, tool)
    if not os.path.isdir(tool_root_dir):
        os.mkdir(tool_root_dir)
    # check circleseq directories
    circleseq_dir = os.path.join(tool_root_dir, EXPERIMENTTYPE[0])
    if not os.path.isdir(circleseq_dir):
        os.mkdir(circleseq_dir)
    _create_celltype_dirtree(circleseq_dir)  # create cell types dir tree
    # check guideseq directories
    guideseq_dir = os.path.join(tool_root_dir, EXPERIMENTTYPE[1])
    if not os.path.isdir(guideseq_dir):
        os.mkdir(guideseq_dir)
    _create_celltype_dirtree(guideseq_dir)  # create cell types dir tree


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
            raise OSError("An error ocuurred while running %s" % (tool))


def run_mutect2(exp_type, offregion, threads):
    """Run Mutect2 to call edits on the input regions

    :param exp_type: validation experiment type <circleseq, guideseq>
    :type exp_type: str
    :param offregion: analyze offregions
    :type offregion: bool
    :param threads: threads
    :type threads: int
    """
    commands = []  # commands array
    for guide in GUIDES:
        outdir = os.path.join(OUTDIR, "mutect2")  # in mutect2 results directory
        if exp_type == "guideseq":
            targets = os.path.join(GUIDESEQ, "%s.guideseq" % (guide))
            outdir = os.path.join(outdir, "guideseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
            name_col = 6
        else:  # circleseq
            targets = os.path.join(
                CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide)
            )
            outdir = os.path.join(outdir, "circleseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
            name_col = 12
        for cell_type in CELLTYPES:
            cmd = str(
                "python3 %s --targets %s --genome %s --bam1 %s --bam2 %s "
                "--normal DNMT1Site3 --chrom-col %d --start-col %d --stop-col %d "
                "--name-col %d %s --out %s --threads %d"
            )
            if cell_type == "GM12878":
                if guide == "EMX1":  # CRAM file, others in BAM format
                    bam1 = os.path.join(BAMS, "%s.cram" % (guide))
                else:
                    bam1 = os.path.join(BAMS, "%s.bam" % (guide))
                bam2 = os.path.join(BAMS, "DNMT1Site3.bam")
            else:  # cell type K562
                bam1 = os.path.join(BAMS, "%s_%s.cram" % (cell_type, guide))
                bam2 = os.path.join(BAMS, "K562_DNMT1Site3.cram")
            odir = os.path.join(outdir, cell_type)
            offr = ""
            if offregion:
                odir = os.path.join(odir, "offregion", guide)
                offr = "--offregion"  # shift target sites
            else:
                odir = os.path.join(odir, "onregion", guide)
            commands.append(
                cmd
                % (
                    MUTECTPY,
                    targets,
                    GENOME,
                    bam1,
                    bam2,
                    chrom_col,
                    start_col,
                    stop_col,
                    name_col,
                    offr,
                    odir,
                    threads,
                )
            )
    # run variant calling
    run_commands(commands, threads, VCALLINGTOOLS[0])


def run_strelka(exp_type, offregion, threads):
    """Run Strelka to call edits on the input regions

    :param exp_type: validation experiment type <circleseq, guideseq>
    :type exp_type: str
    :param offregion: analyze offregions
    :type offregion: bool
    :param threads: threads
    :type threads: int
    """

    commands = []  # commands array
    if not os.path.isdir(STRELKARUNDIR):
        os.mkdir(STRELKARUNDIR)
    for guide in GUIDES:
        outdir = os.path.join(OUTDIR, "strelka")
        if exp_type == "guideseq":
            targets = os.path.join(GUIDESEQ, "%s.guideseq" % (guide))
            outdir = os.path.join(outdir, "guideseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
            name_col = 6
        else:  # circleseq
            targets = os.path.join(
                CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide)
            )
            outdir = os.path.join(outdir, "circleseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
            name_col = 12
        for cell_type in CELLTYPES:
            cmd = str(
                "python %s --targets %s --genome %s --normal-bam %s --tumor-bam "
                "%s --chrom-col %d --start-col %d --stop-col %d --name-col %d %s "
                "--run-dir %s --out %s"
            )
            if cell_type == "GM12878":
                if guide == "EMX1":  # CRAM file, others in BAM format
                    bam1 = os.path.join(BAMS, "%s.cram" % (guide))
                else:
                    bam1 = os.path.join(BAMS, "%s.bam" % (guide))
                bam2 = os.path.join(BAMS, "DNMT1Site3.bam")
            else:  # K562 cell type
                bam1 = os.path.join(BAMS, "%s_%s.cram" % (cell_type, guide))
                bam2 = os.path.join(BAMS, "K562_DNMT1Site3.cram")
            odir = os.path.join(outdir, cell_type)
            offr = ""
            if offregion:
                odir = os.path.join(odir, "offregion", guide)
                offr = "--offregion"
            else:
                odir = os.path.join(odir, "onregion", guide)
            commands.append(
                cmd
                % (
                    STRELKAPY,
                    targets,
                    GENOME,
                    bam1,
                    bam2,
                    chrom_col,
                    start_col,
                    stop_col,
                    name_col,
                    offr,
                    STRELKARUNDIR,
                    odir,
                )
            )
    # run variant calling
    run_commands(commands, threads, VCALLINGTOOLS[1])
    # remove strelka tmp wd
    code = subprocess.call("rm -rf %s" % (STRELKARUNDIR), shell=True)
    if code != 0:
        raise OSError("An error ocuurred while running %s" % (cmd))


def run_pindel(exp_type, offregion, threads):
    """Run Pindel to call edits on the input regions

    :param exp_type: validation experiment type <circleseq, guideseq>
    :type exp_type: str
    :param offregion: analyze offregions
    :type offregion: bool
    :param threads: threads
    :type threads: int
    """
    commands = []  # commands array
    for guide in GUIDES:
        outdir = os.path.join(OUTDIR, "pindel")
        if exp_type == "guideseq":
            targets = os.path.join(GUIDESEQ, "%s.guideseq" % (guide))
            outdir = os.path.join(outdir, "guideseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
        else:  # circleseq
            targets = os.path.join(
                CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide)
            )
            outdir = os.path.join(outdir, "circleseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
        for cell_type in CELLTYPES:
            cmd = str(
                "python %s --targets %s --genome %s --normal-bam %s --normal-sample "
                "%s --tumor-bam %s --tumor-sample %s --chrom-col %d --start-col %d "
                "--stop-col %d %s --out %s --threads %d"
            )
            if cell_type == "GM12878":
                if guide == "EMX1":  # CRAM file, others in BAM format
                    bam1 = os.path.join(BAMS, "%s.cram" % (guide))
                else:
                    bam1 = os.path.join(BAMS, "%s.bam" % (guide))
                bam2 = os.path.join(BAMS, "DNMT1Site3.bam")
            else:  # K562 cell type
                bam1 = os.path.join(BAMS, "%s_%s.cram" % (cell_type, guide))
                bam2 = os.path.join(BAMS, "K562_DNMT1Site3.cram")
            odir = os.path.join(outdir, cell_type)
            offr = ""
            if offregion:
                odir = os.path.join(odir, "offregion", guide)
                offr = "--offregion"
            else:
                odir = os.path.join(odir, "onregion", guide)
            commands.append(
                cmd
                % (
                    PINDELPY,
                    targets,
                    GENOME,
                    bam2,
                    "DNMT1Site3",
                    bam1,
                    guide,
                    chrom_col,
                    start_col,
                    stop_col,
                    offr,
                    odir,
                    threads,
                )
            )
    # run variant calling
    run_commands(commands, threads, VCALLINGTOOLS[2])


def run_varscan(exp_type, offregion, threads):
    """Run Varscan to call edits on the input regions

    :param exp_type: validation experiment type <circleseq, guideseq>
    :type exp_type: str
    :param offregion: analyze offregions
    :type offregion: bool
    :param threads: threads
    :type threads: int
    """
    commands = []  # commands array
    for guide in GUIDES:
        outdir = os.path.join(OUTDIR, "varscan")
        if exp_type == "guideseq":
            targets = os.path.join(GUIDESEQ, "%s.guideseq" % (guide))
            outdir = os.path.join(outdir, "guideseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
        else:  # circleseq
            targets = os.path.join(
                CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide)
            )
            outdir = os.path.join(outdir, "circleseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
        for cell_type in CELLTYPES:
            cmd = str(
                "python %s --targets %s --genome %s --normal-bam %s --tumor-bam %s "
                "--chrom-col %d --start-col %d --stop-col %d %s --out %s"
            )
            if cell_type == "GM12878":
                if guide == "EMX1":  # CRAM file, others in BAM format
                    bam1 = os.path.join(BAMS, "%s.cram" % (guide))
                else:
                    bam1 = os.path.join(BAMS, "%s.bam" % (guide))
                bam2 = os.path.join(BAMS, "DNMT1Site3.bam")
            else:  # K562 cell type
                bam1 = os.path.join(BAMS, "%s_%s.cram" % (cell_type, guide))
                bam2 = os.path.join(BAMS, "K562_DNMT1Site3.cram")
            odir = os.path.join(outdir, cell_type)
            offr = ""
            if offregion:
                odir = os.path.join(odir, "offregion", guide)
                offr = "--offregion"
            else:
                odir = os.path.join(odir, "onregion", guide)
            commands.append(
                cmd
                % (
                    VARSCANPY,
                    targets,
                    GENOME,
                    bam2,
                    bam1,
                    chrom_col,
                    start_col,
                    stop_col,
                    offr,
                    odir,
                )
            )
    # run variant calling
    run_commands(commands, threads, VCALLINGTOOLS[3])


def _runinfo(args):
    """(PRIVATE)
    Print run information

    :param args: input arguments
    :type args: arparse.ArgumentParser
    """

    sys.stderr.write("-- RUN INFO --\n\n")
    sys.stderr.write("\tTOOL:\t%s\n" % (args.tool))
    sys.stderr.write("\tTYPE:\t%s\n" % (args.type))
    sys.stderr.write("\tOFFREGION:\t%s\n" % (args.offregion))
    sys.stderr.write("\tTHREADS:\t%d\n\n" % (args.threads))


def main():
    args = parse_commandline()
    # check input arguments consistency
    if args.tool not in VCALLINGTOOLS:
        raise ValueError(
            "%s cannot run %s. Check the help for the available tools"
            % (__file__, args.tool)
        )
    if args.type not in EXPERIMENTTYPE:
        raise ValueError(
            "Forbidden target sites validation experiment type. Check the help "
            "for the available values"
        )
    if args.threads < 0:
        raise ValueError("Forbidden number of threads selected (%d)" % (args.threads))
    if args.threads == 0:  # autodetect
        args.threads = multiprocessing.cpu_count()
    _runinfo(args)  # print run info to stderr
    if args.tool == VCALLINGTOOLS[0]:  # mutect2
        create_result_dirtree(args.tool)  # create mutect2 result directory tree
        run_mutect2(args.type, args.offregion, args.threads)  # run mutect2
    elif args.tool == VCALLINGTOOLS[1]:  # strelka
        create_result_dirtree(args.tool)  # create strelka result directory tree
        run_strelka(args.type, args.offregion)  # run strelka
    elif args.tool == VCALLINGTOOLS[2]:  # pindel
        create_result_dirtree(args.tool)  # create pindel result directory tree
        run_pindel(args.type, args.offregion)  # run pindel
    elif args.tool == VCALLINGTOOLS[3]:  # varscan
        create_result_dirtree(args.tool)  # create varscan result directory tree
        run_varscan(args.type, args.offregion)  # run varscan
    else:
        raise ValueError(
            "%s cannot run %s. Check the help for the available tools"
            % (__file__, args.tool)
        )


if __name__ == "__main__":
    main()
