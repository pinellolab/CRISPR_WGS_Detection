"""Run Varscan on the set of off-target regions identified by CasOffinder. 
The regions are specified in the input file.
"""

from targets_file_parser import parse_targets_coordinates

from glob import glob

import subprocess
import argparse
import tempfile
import os

SAMTOOLS = "samtools mpileup"
SAMTOOLSTMP = tempfile.mkdtemp()  # stores temporary MPILEUP files
VARSCAN = "varscan somatic"


def parse_commandline():
    """Parse the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified input regions",
        usage="\n\tpython %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--normal-bam <BAM> --tumor-bam <BAM> --out <OUTDIR>",
    )
    parser.add_argument(
        "--targets", type=str, metavar="TARGETS-FILE", help="Targets coordinate file"
    )
    parser.add_argument("--genome", type=str, metavar="GENOME", help="Reference genome")
    parser.add_argument(
        "--normal-bam",
        type=str,
        metavar="BAM",
        dest="normal_bam",
        help="Normal BAM file",
    )
    parser.add_argument(
        "--tumor-bam", type=str, metavar="BAM", dest="tumor_bam", help="Tumor BAM file"
    )
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
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    return parser.parse_args()  # parse command line


def mpileup(genome, region, bam, outfile):
    """Computes mpileup files using SAMtools for the input BAM. The resulting
    mpileup file covers the input region.

    :param genome: reference genome
    :type genome: str
    :param region: region
    :type region: str
    :param bam: BAM file
    :type bam: str
    :param outfile: output mpileup file
    :type outfile: str
    :raises OSError: raise if mpileup construction fails
    :return: mpileup file
    :rtype: str
    """
    print(genome, 0, region, bam, outfile)
    code = subprocess.call(
        "%s -f %s -d %d -r %s %s > %s" % (SAMTOOLS, genome, 0, region, bam, outfile),
        shell=True,
    )
    if code != 0:
        raise OSError("An error occurred while creating mpileup file for %s" % (region))
    return outfile


def run_varscan(genome, targets, normal_bam, tumor_bam, outdir):
    """Wrapper function to call edits on the input regions using VarScan

    :param genome: reference genome
    :type genome: str
    :param targets: target regions
    :type targets: List[str]
    :param normal_bam: normal BAM
    :type normal_bam: str
    :param tumor_bam: tumor BAM
    :type tumor_bam: str
    :param outdir: output directory
    :type outdir: str
    :raises OSError: raise if Varscan fails
    :raises OSError: raise if mpileup files deletion fails
    """
    for target in targets:
        # compute mpileup files for normal and tumor BAMs
        mpileup_normal = mpileup(
            genome,
            target,
            normal_bam,
            os.path.join(SAMTOOLSTMP, "%s_normal.mpileup" % (target)),
        )
        mpileup_tumor = mpileup(
            genome,
            target,
            tumor_bam,
            os.path.join(SAMTOOLSTMP, "%s_tumor.mpileup" % (target)),
        )
        # call edits with varscan
        outfile = os.path.join(outdir, target)
        code = subprocess.call(
            "%s %s %s %s --min-var-freq 0.0001 --output-vcf 1 --strand-filter 1"
            % (VARSCAN, mpileup_normal, mpileup_tumor, outfile),
            shell=True,
        )
        if code != 0:
            raise OSError("An error occurred while calling edits on %s" % (target))
        # remove mpileup files
        code = subprocess.call("rm %s %s" % (mpileup_normal, mpileup_tumor), shell=True)
        if code != 0:
            raise OSError(
                "Unable to delete %s and %s" % (mpileup_normal, mpileup_tumor)
            )


def main():
    args = parse_commandline()
    targets = parse_targets_coordinates(
        args.targets, args.offtarget_upstream, args.offtarget_downstream
    )
    run_varscan(args.genome, targets, args.normal_bam, args.tumor_bam, args.out)


if __name__ == "__main__":
    main()
