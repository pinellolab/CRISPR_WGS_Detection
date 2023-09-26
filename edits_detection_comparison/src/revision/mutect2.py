"""
"""

from command_runners import run_commands
from targets import read_targets_file
from utils import bgzip, tabix

from glob import glob

import os

MUTECT2 = "gatk Mutect2"  # Mutect2

def mutect2(genome, bam1, bam2, normal, targets, threads, outdir):
    """Detect rare genome editing events using GATK Mutect2

    :param genome: genome
    :type genome: str
    :param bam1: normal BAM
    :type bam1: str
    :param bam2: tumor BAM
    :type bam2: str
    :param normal: normal sample
    :type normal: str
    :param targets: target sites
    :type targets: List[str]
    :param threads: threads
    :type threads: int
    :param outdir: output directory
    :type outdir: str
    """
    commands = [
        "%s -R %s -I %s -I %s -normal %s -L %s -O %s > /dev/null" % (
            MUTECT2, genome, bam1, bam2, normal, target, os.path.join(outdir, "%s.vcf" % (target))
        )
        for target in targets
    ]
    run_commands(commands, threads)


def detect_edits_mutect2(
    targets_file,
    genome,
    bam1,
    bam2,
    normal,
    outdir,
    casoffinder, 
    circleseq, 
    guideseq, 
    upstream, 
    downstream,
    threads
):
    # read targets from input file
    targets = read_targets_file(
        targets_file, casoffinder, circleseq, guideseq, upstream, downstream
    )
    # detect editing events with mutect2
    mutect2(genome, bam1, bam2, normal, targets, threads, outdir)
    vcfs = glob(os.path.join(outdir, "*.vcf"))
    assert len(targets) == len(vcfs)  # one VCF per target site
    # compress and index VCFs 
    bgzip(vcfs, threads)
    tabix(glob(os.path.join(outdir, "*.vcf.gz")), threads)



