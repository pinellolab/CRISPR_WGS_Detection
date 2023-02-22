"""Run Strelka on the set of regions specified in the input file.
"""

from tqdm import tqdm
from glob import glob

import subprocess
import argparse
import sys
import os


STRELKACONFIG = "configureStrelkaSomaticWorkflow.py"  # strelka config file
STRELKARUN = "runWorkflow.py"  # strelka workflow run
PADSIZE = 10000  # region padding size (10kb)


def parse_commandline():
    """The function parses the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified regions using Strelka",
        usage="\n\tpython3 %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--normal-bam <BAM> --tumor-bam <BAM> --chrom-col <CHR-COL-NUM> "
        "--start-col <START-COL-NUM> --stop-col <STOP-COL-NUM> --name-col "
        "<NAME-COL-NUM> --run-dir <RUN-DIR> --out <OUTDIR>",
    )
    parser.add_argument(
        "--targets", type=str, metavar="TARGETS-FILE", help="Target sites file"
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
        "--run-dir",
        type=str,
        metavar="RUN-DIR",
        dest="run_dir",
        help="Temporary run directory",
    )
    parser.add_argument(
        "--chrom-col",
        type=int,
        metavar="CHR-COL-NUM",
        help="Column containing chromosome",
        dest="chrom_col",
    )
    parser.add_argument(
        "--start-col",
        type=int,
        metavar="START-COL-NUM",
        help="Column containing start coordinate",
        dest="start_col",
    )
    parser.add_argument(
        "--stop-col",
        type=int,
        metavar="STOP-COL-NUM",
        help="Column containing stop coordinate",
        dest="stop_col",
    )
    parser.add_argument(
        "--name-col",
        type=int,
        metavar="NAME-COL-NUM",
        help="Column containing target site name",
        dest="name_col",
    )
    parser.add_argument(
        "--offregion",
        action="store_true",
        default=False,
        help="Shift the target regions 100bp upstream and downstream (not "
        "overlapping the original site)",
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    return parser.parse_args()

def parse_targets(targetfile):
    """Parse the guide targets file

    :param targetfile: guide targets file
    :type targetfile: str
    :raises OSError: raise on read() failure
    :return: targets file lines
    :rtype: List[str]
    """
    try:
        handle = open(targetfile, mode="r")
        lines = [line.strip().split() for i, line in enumerate(handle) if i > 0]
    except OSError:
        raise OSError("An error occurred while reading %s" % (targetfile))
    finally:
        handle.close()
    return lines


def compute_regions(lines, chrom_col, start_col, stop_col, offregion):
    """Pad the target sites by 10kb upstream and downstream

    :param lines: targets file lines
    :type lines: List[str]
    :param chrom_col: chromosome data columns
    :type chrom_col: int
    :param start_col: start data columns
    :type start_col: int
    :param stop_col: stop data columns
    :type stop_col: int
    :param offregion: pad off regions
    :type offregion: bool
    :return: padded regions
    :rtype: List[str]
    """
    if offregion:  # upstream and downstream target sites shift
        regions = [
            "%s:%s-%s"
            % (
                line[chrom_col],
                (int(line[start_col]) - 100 - PADSIZE),
                (int(line[start_col]) - 100),
            )
            for line in lines
        ] + [
            "%s:%s-%s"
            % (
                line[chrom_col],
                (int(line[stop_col]) + 100),
                (int(line[stop_col]) + 100 + PADSIZE),
            )
            for line in lines
        ]
        assert len(regions) == (len(lines) * 2)
    else:
        regions = [
            "%s:%s-%s"
            % (
                line[chrom_col],
                (int(line[start_col]) - PADSIZE),
                (int(line[stop_col]) + PADSIZE),
            )
            for line in lines
        ]
        assert len(lines) == len(regions)
    return regions


def get_names(regions):
    """Recover target sites names

    :param regions: padded target sites
    :type regions: List[str]
    :return: target sites names 
    :rtype: List[str]
    """
    names_list = [region.replace(":", "_").replace("-", "_") for region in regions]
    assert len(names_list) == len(regions)
    return names_list


def strelka(normal_bam, tumor_bam, genome, region, name, rundir, outdir):
    """Run Strelka on the input target site

    :param normal_bam: BAM file
    :type normal_bam: str
    :param tumor_bam: BAM file
    :type tumor_bam: str
    :param genome: genome
    :type genome: str
    :param region: target site
    :type region: str
    :param name: region name
    :type name: str
    :param rundir: run directory
    :type rundir: str
    :param outdir: output directory
    :type outdir: str
    :raises OSError: raise on strelka config failure
    :raises OSError: raise on strelka analysis failure
    :raises OSError: raise on gunzip failure
    :raises OSError: raise on mv failure
    :raises OSError: raise on mv failure
    """
    # run config script
    cmd = (
        "%s --normalBam %s --tumorBam %s --referenceFasta %s --region %s --exome --runDir %s"
        % (STRELKACONFIG, normal_bam, tumor_bam, genome, region, rundir)
    )
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error occurred while running %s" % (cmd))
    # call variants
    cmd = "python %s -m local --quiet" % (os.path.join(rundir, STRELKARUN))
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error occurred while running %s" % (cmd))
    # recover results
    vcfsgz = glob(os.path.join(rundir, "results/variants/*.vcf.gz"))
    for vcfgz in vcfsgz:
        code = subprocess.call("gunzip %s" % (vcfgz), shell=True)
        if code != 0:
            raise OSError("An error occurred while running %s" % (cmd))
    vcfs = glob(os.path.join(rundir, "results/variants/*.vcf"))
    for vcf in vcfs:
        cmd = "mv %s %s" % (vcf, outdir)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error occurred while running %s" % (cmd))
    vcfs = glob(os.path.join(outdir, "somatic.*.vcf"))
    for vcf in vcfs:
        cmd = "mv %s %s" % (
            vcf,
            os.path.join(outdir, "_".join([name, os.path.basename(vcf)])),
        )
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error occurred while running %s" % (cmd))
    # remove temporary results
    tmpdata = os.listdir(rundir)
    for d in tmpdata:
        code = subprocess.call("rm -rf %s" % (os.path.join(rundir, d)), shell=True)
        if code != 0:
            raise OSError("An error occurred while running %s" % (cmd))


def main():
    args = parse_commandline()
    # parse input targets file
    targets = parse_targets(args.targets)
    # recover regions
    regions = compute_regions(
        targets, args.chrom_col, args.start_col, args.stop_col, args.offregion
    )
    # run strelka
    names = get_names(regions)
    for i, region in tqdm(enumerate(regions)):
        strelka(
            args.normal_bam,
            args.tumor_bam,
            args.genome,
            region,
            names[i],
            args.run_dir,
            args.out,
        )


if __name__ == "__main__":
    main()
