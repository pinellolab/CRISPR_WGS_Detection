"""Utilities functions and variables
"""

from command_runners import run_commands

import os

VCALLINGTOOLS = ["mutect2", "strelka", "varscan"]
ONTARGETWIDTH = 100  # 100 bp on-target region
OFFTARGETWIDTH = 10000  # 10 kb off-target regions

BGZIP = "bgzip"
TABIX = "tabix"

def bgzip(vcf_files, threads):
    """Compress VCF file

    :param vcf_files: VCF
    :type vcf_files: str
    :param threads: threads
    :type threads: int
    """
    commands = ["%s %s" % (BGZIP, vcf_file) for vcf_file in vcf_files]
    run_commands(commands, threads)

def tabix(vcf_files, threads):
    """Index VCF file

    :param vcf_files: VCF
    :type vcf_files: str
    :param threads: threads
    :type threads: int
    """
    commands = ["%s %s" % (TABIX, vcf_file) for vcf_file in vcf_files]
    run_commands(commands, threads)

def check_args(args):
    if args.tool not in VCALLINGTOOLS:
        raise ValueError("Wrong variant caller (%s)!" % (args.tool))
    if not os.path.isfile(args.targets):
        raise FileNotFoundError("Targets file not found!")
    if not os.path.isfile(args.genome):
        raise FileNotFoundError("Genome FASTA not found!") 
    if not os.path.isfile(args.bam1):
        raise FileNotFoundError("BAM file (%s) not found!" % (args.bam1))
    if not os.path.isfile(args.bam2):
        raise FileNotFoundError("BAM file (%s) not found!" % (args.bam2))
    if args.tool == VCALLINGTOOLS[0] and not bool(args.normal_sample):  # mutect2
        raise ValueError("Normal sample name not given!")
    if not os.path.isdir(args.out):  
        os.mkdir(args.out)  # create directory
    if args.threads < 0 or not isinstance(args.threads, int):
        raise ValueError("Forbidden threads number (%s)" % (args.threads))
    