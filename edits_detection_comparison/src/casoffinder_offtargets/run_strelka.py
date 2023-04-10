"""Run Strelka on the set of off-target regions identified by CasOffinder. The 
regions are specified in the input file.
"""

from targets_file_parser import parse_targets_coordinates

from glob import glob

import subprocess
import argparse
import os

STRELKACONFIG = "configureStrelkaSomaticWorkflow.py"
STRELKARUN = "runWorkflow.py"

def parse_commandline():
    """Parse the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified regions using Strelka",
        usage="\n\tpython %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--normal-bam <BAM> --tumor-bam <BAM> --run-dir <RUN-DIR> --out <OUTDIR> ",
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
        "--tumor-bam", type=str, metavar="BAM", dest="tumor_bam", help="Tumor BAM file"
    )
    parser.add_argument(
        "--run-dir", type=str, metavar="RUN-DIR", dest="run_dir", help="Temporary running directory"
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    return parser.parse_args()  # parse command line

def config(normal_bam, tumor_bam, genome, region, rundir):
    """Create Strelka's config file

    :param normal_bam: normal BAM
    :type normal_bam: str
    :param tumor_bam: tumor BAM
    :type tumor_bam: str
    :param genome: reference genome
    :type genome: str
    :param region: region
    :type region: str
    :param rundir: temporary run directory
    :type rundir: str
    :raises OSError: raise if config file creation fails
    """
    cmd = (
        "%s --normalBam %s --tumorBam %s --referenceFasta %s --region %s --exome --runDir %s" 
        % (STRELKACONFIG, normal_bam, tumor_bam, genome, region, rundir)
    )
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError(
            "An error occurred while building config script for region %s" % (region)
        )
    
def recover_results(rundir, region, outdir):
    """Recover the VCFs created by Strelka

    :param rundir: temporary run directory
    :type rundir: str
    :param region: region
    :type region: str
    :param outdir: output directory
    :type outdir: str
    :raises OSError: raise if VCF unzipping fails
    :raises OSError: raise if unzipped VCF recovery fails
    :raises OSError: raise if VCF renaming fails
    :raises OSError: raise if temporary directory removal fails
    """
    # recover VCF.GZ files and unzip VCFs
    for vcfgz in glob(os.path.join(rundir, "results/variants/*.vcf.gz")):
        code = subprocess.call("gunzip %s" % (vcfgz), shell=True)
        if code != 0:
            raise OSError("Unable to unzip %s" % (vcfgz))
    for vcf in glob(os.path.join(rundir, "results/variants/*.vcf")):
        code = subprocess.call("mv %s %s" % (vcf, outdir), shell=True)
        if code != 0:
            raise OSError("Unable to recover %s" % (vcf))
    for vcf in glob(os.path.join(outdir, "somatic.*.vcf")):
        subprocess.call(
            "mv %s %s" % (
                vcf,
                os.path.join(outdir, "_".join([region, os.path.basename(vcf)]))
            ),
            shell=True
        )
        if code != 0:
            raise OSError("An error occurred while renaming %s" % (vcf))
    # remove temporary data
    for d in os.listdir(rundir):
        code = subprocess.call("rm -rf %s" % (os.path.join(rundir, d)), shell=True)
        if code != 0:
            raise OSError("An error occurred while removing %s" % (d))


def run_strelka(normal_bam, tumor_bam, genome, targets, rundir, outdir):
    """Wrapper function to call edits on the input regions using Strelka

    :param normal_bam: normal BAM
    :type normal_bam: str
    :param tumor_bam: tumor BAM
    :type tumor_bam: str
    :param genome: reference genome
    :type genome: str
    :param targets: target regions
    :type targets: List[str]
    :param rundir: temporary run directory
    :type rundir: str
    :param outdir: output directory
    :type outdir: str
    :raises OSError: raise if Strelka fails
    """
    for target in targets:
        config(normal_bam, tumor_bam, genome, target, rundir)  # run config script
        # call edits
        cmd = "python %s -m local --quiet" % (os.path.join(rundir, STRELKARUN))
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error occurred while calling edits on %s" % (target))
        recover_results(rundir, target, outdir)
    
def main():
    args = parse_commandline()
    targets = parse_targets_coordinates(args.targets)
    run_strelka(args.normal_bam, args.tumor_bam, args.genome, targets, args.run_dir, args.out)
    assert len(targets) == len(glob(os.path.join(args.out, "*.vcf")))

if __name__ == "__main__":
    main()
 