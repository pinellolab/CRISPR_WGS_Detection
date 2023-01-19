"""This python script calls different variant calling tools to detect 
mutations induced by CRISPR-CAs9 editing.

NB. Strelka runs only on Python2, therefore the script syntax is suited to 
be compatible with both python3 and python2.
"""

import subprocess
import argparse
import tempfile
import sys
import os


SCRIPTS = "/PHShome/mi825/Desktop/wgs_crisprcas9/src/"
MUTECTPY = os.path.join(SCRIPTS, "run_mutect.py")
STRELKAPY = os.path.join(SCRIPTS, "run_strelka.py")
PINDELPY = os.path.join(SCRIPTS, "run_pindel.py")
VARSCANPY = os.path.join(SCRIPTS, "run_varscan.py")
# GUIDES = ["EMX1", "HEKSite4", "RNF2", "VEGFASite3"]
GUIDES = ["VEGFASite3"]
VCALLINGTOOLS = ["mutect2", "strelka", "pindel", "varscan"]
# CELLTYPES = ["GM12878", "K562"]
CELLTYPES = ["K562"]
BASEDIR = "/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/"
GUIDESEQ = os.path.join(
    BASEDIR, "offtargetDetection/casoffinder/offby6/CRISPRessoWGS/guideseq_anno"
)
STRELKARUNDIR = tempfile.mkdtemp()
CIRCLESEQ = os.path.join(BASEDIR, "offtargetDetection/circleseq/")
BAMS = os.path.join(BASEDIR, "wgs/GM12878-Cas9/WGS1000/data/")
GENOME = "/data/pinello/COMMON_DATA/REFERENCE_GENOMES/Broad/hg38/Homo_sapiens_assembly38.fasta"
OUTDIR = "/PHShome/mi825/Desktop/wgs_crisprcas9/VCFs"


def parse_commandline():
    """The function parses the input command line arguments."""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Run variant calling tools to detect CRISPR-Cas9 edits",
        usage="\n\tpython3 %(prog)s --tool <TOOL-NAME> --type <TYPE>"
    )
    parser.add_argument(
        "--tool", 
        type=str, 
        metavar="TOOL-NAME", 
        help="Variant calling tool. Available values: <mutect2, strelka, "
             "varscan, pindel>"
    )
    parser.add_argument(
        "--type", 
        type=str, 
        metavar="TYPE",
        help="Target sites validation experiment type. Available values: "
             "<guideseq, circleseq>"
    )
    parser.add_argument(
        "--offregion",
        action="store_true",
        default=False,
        help="Shift the target regions 100bp upstream and downstream (not "
             "overlapping the original target site)"
    )
    args = parser.parse_args()
    # check arguments consistency
    if args.tool not in VCALLINGTOOLS:
        raise ValueError(
            "%s cannot run %s. Check the help for the available tools" % (
                __file__, args.tool
            )
        )
    if args.type != "circleseq" and args.type != "guideseq":
        raise ValueError(
            "Forbidden target sites validation experiment type. Check the help "
            "for the available values"
        )
    return args


def run_mutect2(exp_type, offregion):
    """The function build the commands to run Mutect2."""

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
            targets = os.path.join(CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide))
            outdir = os.path.join(outdir, "circleseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
            name_col = 12
        for cell_type in CELLTYPES:
            cmd = str(
                "python3 %s --targets %s --genome %s --bam1 %s --bam2 %s "
                "--normal DNMT1Site3 --chrom-col %d --start-col %d --stop-col %d "
                "--name-col %d %s --out %s"
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
                cmd % (
                    MUTECTPY, targets, GENOME, bam1, bam2, chrom_col, start_col, stop_col, name_col, offr, odir
                )
            )
    # run variant calling
    for cmd in commands:
        sys.stderr.write("\n\n%s\n\n" % (cmd))  # TODO: remove this line
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error ocuurred while running %s" % (cmd))


def run_strelka(exp_type, offregion):
    """The function builds the commands to run strelka."""
    
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
            targets = os.path.join(CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide))
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
                cmd % (
                    STRELKAPY, targets, GENOME, bam1, bam2, chrom_col, start_col, stop_col, name_col, offr, STRELKARUNDIR, odir
                )
            )
    # run variant calling
    for cmd in commands:
        sys.stderr.write("\n\n%s\n\n" % (cmd))  # TODO: remove this line
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error ocuurred while running %s" % (cmd))
    code = subprocess.call("rm -rf %s" % (STRELKARUNDIR), shell=True)
    if code != 0:
        raise OSError("An error ocuurred while running %s" % (cmd))


def run_pindel(exp_type, offregion):
    """The function builds the commands to run pindel."""

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
            targets = os.path.join(CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide))
            outdir = os.path.join(outdir, "circleseq")
            chrom_col = 0
            start_col = 1
            stop_col = 2
        for cell_type in CELLTYPES:
            cmd = str(
                "python %s --targets %s --genome %s --normal-bam %s --normal-sample "
                "%s --tumor-bam %s --tumor-sample %s --chrom-col %d --start-col %d "
                "--stop-col %d %s --out %s"
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
                cmd % (
                    PINDELPY, targets, GENOME, bam2, "DNMT1Site3", bam1, guide, chrom_col, start_col, stop_col, offr, odir
                )
            )
    # run variant calling
    for cmd in commands:
        sys.stderr.write("\n\n%s\n\n" % (cmd))  # TODO: remove this line
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error ocuurred while running %s" % (cmd))


def run_varscan(exp_type, offregion):
    """The function builds the command to run varscan."""

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
            targets = os.path.join(CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide))
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
                cmd % (
                    VARSCANPY, targets, GENOME, bam2, bam1, chrom_col, start_col, stop_col, offr, odir
                )
            )
    # run variant calling
    for cmd in commands:
        sys.stderr.write("\n\n%s\n\n" % (cmd))  # TODO: remove this line
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error ocuurred while running %s" % (cmd))


def __runinfo(args):
    """Print the current run info."""

    sys.stderr.write("-- RUN INFO --\n\n")
    sys.stderr.write("\tTOOL:\t%s\n" % (args.tool))
    sys.stderr.write("\tTYPE:\t%s\n" % (args.type))
    sys.stderr.write("\tOFFREGION:\t%s\n\n" % (args.offregion))


def main():
    args = parse_commandline()
    __runinfo(args)
    if args.tool == VCALLINGTOOLS[0]:  # mutect2
        run_mutect2(args.type, args.offregion)  # run mutect2
    elif args.tool == VCALLINGTOOLS[1]:  # strelka
        run_strelka(args.type, args.offregion)  # run strelka
    elif args.tool == VCALLINGTOOLS[2]:  # pindel
        run_pindel(args.type, args.offregion)  # run pindel
    elif args.tool == VCALLINGTOOLS[3]:  # varscan
        run_varscan(args.type, args.offregion)  # run varscan
    else:
        raise ValueError(
            "%s cannot run %s. Check the help for the available tools" % (
                __file__, args.tool 
            )
        )


if __name__ == "__main__":
    main()

