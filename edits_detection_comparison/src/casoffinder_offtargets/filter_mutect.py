"""Filter edits called by Mutect2, using FilterMutectCalls functionality from 
GATK toolkit
"""

from glob import glob

import subprocess
import argparse
import os

FILTER = "gatk FilterMutectCalls"


def parse_commandline():
    """Parse the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to filter edits called by Mutect2",
        usage="\n\tpython %(prog)s --input-dir <INPUT-DIR> --output-dir "
        "<OUTPUT-DIR> --genome <GENOME>",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        metavar="INPUT-DIR",
        dest="input_dir",
        required=True,
        help="Folder storing the VCF files returned by Mutect2",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        metavar="OUTPUT-DIR",
        dest="output_dir",
        required=True,
        help="Folder where the filtered VCFs should be stored",
    )
    parser.add_argument(
        "--genome", type=str, metavar="GENOME", required=True, help="Reference genome"
    )
    return parser.parse_args()


def recover_vcfs(vcfdir):
    """Recover VCF files to filter

    :return: VCF files
    :rtype: List[str]
    """
    if not os.path.isdir(vcfdir):
        raise FileNotFoundError("Unable to locate %s" % (vcfdir))
    vcfs = glob(os.path.join(vcfdir, "*.vcf"))  # collect VCFs
    assert all([os.path.isfile(f) for f in vcfs])  # check VCF files existence
    return vcfs


def filter_mutect(vcf, genome, vcfout):
    """Filter edits called by Mutect2"""
    code = subprocess.call(
        "%s -V %s -R %s -O %s" % (FILTER, vcf, genome, vcfout), shell=True
    )
    if code != 0:
        raise OSError("Mutect2 filtering failed on %s" % (vcf))


def main():
    args = parse_commandline()  # parse commandline
    vcffiles = recover_vcfs(args.input_dir)  # recover VCFs
    if not os.path.exists(args.output_dir):  # if does not exist, create the out folder
        os.mkdir(args.output_dir)
    for vcf in vcffiles:
        vcfout = os.path.join(
            args.output_dir, "%s.filtered.vcf" % (os.path.basename(vcf))
        )
        filter_mutect(vcf, args.genome, vcfout)


if __name__ == "__main__":
    main()
