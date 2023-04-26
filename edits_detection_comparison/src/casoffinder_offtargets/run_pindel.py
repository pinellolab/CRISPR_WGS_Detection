"""Run Pindel on the set of off-target regions identified by CasOffinder. The 
regions are specified in the input file.
"""

from targets_file_parser import parse_targets_coordinates

from glob import glob

import subprocess
import tempfile
import argparse
import os

PINDEL = "pindel"
PINDEL2VCF = "pindel2vcf"

def parse_commandline():
    """Parse the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Script to call variants in the specified regions using pindel.",
        usage="\n\tpython3 %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--normal-bam <BAM> --normal-sample <NORMAL-SAMPLE> --tumor-bam <BAM> "
        " --tumor-sample <TUMOR-SAMPLE> --out <OUTDIR>",
    )
    parser.add_argument(
        "--targets", type=str, metavar="TARGETS-FILE", help="Targets coordinates file"
    )
    parser.add_argument(
        "--genome", type=str, metavar="GENOME", help="Reference genome"
    )
    parser.add_argument(
        "--normal-bam", type=str, metavar="BAM", dest="normal_bam", help="Normal BAM file"
    )
    parser.add_argument(
        "--normal-sample", 
        type=str, 
        metavar="NORMAL-SAMPLE", 
        dest="normal_sample", 
        help="Normal sample"
    )
    parser.add_argument(
        "--tumor-bam", type=str, metavar="BAM", dest="tumor_bam", help="Tumor BAM file"
    )
    parser.add_argument(
        "--tumor-sample", 
        type=str, 
        metavar="TUMOR-SAMPLE", 
        dest="tumor_sample", 
        help="Tumor sample"
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    return parser.parse_args()  # parse command line

def config_file(tumor_bam, tumor_sample, normal_bam, normal_sample):
    """Writes the config file required by Pindel to call edits

    :param tumor_bam: tumor BAM
    :type tumor_bam: str
    :param tumor_sample: tumor sample
    :type tumor_sample: str
    :param normal_bam: normal BAM
    :type normal_bam: str
    :param normal_sample: normal sample
    :type normal_sample: str
    :raises OSError: raise if config file writing fails
    :return: config file
    :rtype: str
    """
    configfile = tempfile.NamedTemporaryFile().name
    try:
        with open(configfile, mode="w") as outfile:
            outfile.write("%s\t250\t%s\n" % (tumor_bam, tumor_sample))
            outfile.write("%s\t250\t%s\n" % (normal_bam, normal_sample))
    except OSError:
        raise OSError("An error occurred while writing config file")
    assert os.path.isfile(configfile)
    return os.path.abspath(configfile)

def compute_vcf(genome, outfile):
    """Convert Pindel output file in VCF format

    :param genome: reference genome
    :type genome: str
    :param outfile: output file
    :type outfile: str
    :raises OSError: raise if VCF conversion fails
    :raises OSError: raise if VCF conversion fails
    """
    code = subprocess.call(
        "%s -p %s_D -r %s -R hg38 -d 20131217 -G" % (
            PINDEL2VCF, outfile, genome
        ),
        shell=True
    )
    if code != 0:
        raise OSError("An error occurred while computing %s_D VCF" % (outfile))
    code = subprocess.call(
        "%s -p %s_SI -r %s -R hg38 -d 20131217 -G" % (
            PINDEL2VCF, outfile, genome
        ),
        shell=True
    )
    if code != 0:
        raise OSError("An error occurred while computing %s_SI VCF" % (outfile))

def run_pindel(
    tumor_bam, tumor_sample, normal_bam, normal_sample, targets, outdir, genome
):
    """Wrapper function to call edits on the input regions using Pindel

    :param tumor_bam: tumor BAM
    :type tumor_bam: str
    :param tumor_sample: tumor sample
    :type tumor_sample: str
    :param normal_bam: normal BAM
    :type normal_bam: str
    :param normal_sample: normal sample
    :type normal_sample: str
    :param targets: target regions
    :type targets: List[str]
    :param outdir: output directory
    :type outdir: str
    :param genome: reference genome
    :type genome: str
    :raises OSError: raise if Pindel fails
    """
    configfile = config_file(tumor_bam, tumor_sample, normal_bam, normal_sample)
    for target in targets:
        # call edits
        outfile = os.path.join(outdir, target)
        code = subprocess.call(
            "%s -f %s -i %s -c %s -o %s" % (
                PINDEL, genome, configfile, target, outfile
            ),
            shell=True
        )
        if code != 0:
            raise OSError("An error occurred while calling edits on %s" % (target))
        # convert pindel output in VCF format
        compute_vcf(genome, outfile)

def main():
    args = parse_commandline()
    targets = parse_targets_coordinates(args.targets)
    run_pindel(args.tumor_bam, args.tumor_sample, args.normal_bam, args.normal_sample, targets, args.out, args.genome)

if __name__ == "__main__":
    main()
