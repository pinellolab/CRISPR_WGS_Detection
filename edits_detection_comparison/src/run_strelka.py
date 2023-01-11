"""Run Strelka on the set of regions specified in the input file.
"""

from glob import glob

import subprocess
import argparse
import sys
import os


STRELKACONFIG = "configureStrelkaSomaticWorkflow.py"
STRELKARUN = "runWorkflow.py"
PADSIZE = 10000


def parse_commandline():
    """The function parses the command line arguments.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified regions using Strelka",
        usage="\n\tpython3 %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
              "--normal-bam <BAM> --tumor-bam <BAM> --chrom-col <CHR-COL-NUM> "
              "--start-col <START-COL-NUM> --stop-col <STOP-COL-NUM> --name-col "
              "<NAME-COL-NUM> --run-dir <RUN-DIR> --out <OUTDIR>"
    )
    parser.add_argument(
        "--targets", type=str, metavar="TARGETS-FILE", help="Target sites file"
    )
    parser.add_argument("--genome", type=str, metavar="GENOME", help="Reference genome")
    parser.add_argument(
        "--normal-bam", type=str, metavar="BAM", dest="normal_bam", help="Normal BAM file"
    )
    parser.add_argument(
        "--tumor-bam", type=str, metavar="BAM", dest="tumor_bam", help="Tumor BAM file"
    )
    parser.add_argument(
        "--run-dir", type=str, metavar="RUN-DIR", dest="run_dir", help="Temporary run directory"
    )
    parser.add_argument(
        "--chrom-col", type=int, metavar="CHR-COL-NUM", help="Column containing chromosome", dest="chrom_col"
    )
    parser.add_argument(
        "--start-col", type=int, metavar="START-COL-NUM", help="Column containing start coordinate", dest="start_col"
    )
    parser.add_argument(
        "--stop-col", type=int, metavar="STOP-COL-NUM", help="Column containing stop coordinate", dest="stop_col"
    )
    parser.add_argument(
        "--name-col", type=int, metavar="NAME-COL-NUM", help="Column containing target site name", dest="name_col"
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    args = parser.parse_args()
    return args


def parse_targets(targetfile):
    """The function parses the target sites file."""
    
    try:
        handle = open(targetfile, mode="r")
        lines = [line.strip().split() for i, line in enumerate(handle) if i > 0]
    except OSError:
        raise OSError("An error occurred while reading %s" % (targetfile))
    finally:
        handle.close()
    return lines


def compute_regions(lines, chrom_col, start_col, stop_col):
    """The function computes the padded target regions"""

    regions = [
        "%s:%s-%s" % (
            line[chrom_col], (int(line[start_col]) - PADSIZE), (int(line[start_col]) + PADSIZE)
        )
        for line in lines
    ]
    assert len(lines) == len(regions)
    return regions


def get_names(lines, name_col):
    """The function returns a list containing the targets names."""
     
    names_list = [line[name_col] for line in lines]
    assert len(names_list) == len(lines)
    return names_list


def strelka(normal_bam, tumor_bam, genome, region, name, rundir, outdir):
    """The function runs Strelka on the input region."""

    # run config script
    cmd = "%s --normalBam %s --tumorBam %s --referenceFasta %s --region %s --exome --runDir %s" % (
        STRELKACONFIG, normal_bam, tumor_bam, genome, region, rundir
    )
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error occurred while running \"%s\"" % (cmd))
    # call variants
    cmd = "python %s -m local --quiet" % (os.path.join(rundir, STRELKARUN))
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error occurred while running \"%s\"" % (cmd))
    # recover results
    vcfsgz = glob(os.path.join(rundir, "results/variants/*.vcf.gz"))
    for vcfgz in vcfsgz:
        code = subprocess.call("gunzip %s" % (vcfgz), shell=True)
        if code != 0:
            raise OSError("An error occurred while running \"%s\"" % (cmd))
    vcfs = glob(os.path.join(rundir, "results/variants/*.vcf"))
    for vcf in vcfs:
        cmd = "mv %s %s" % (vcf, outdir)
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error occurred while running \"%s\"" % (cmd))
    vcfs = glob(os.path.join(outdir, "somatic.*.vcf"))
    for vcf in vcfs:
        cmd = "mv %s %s" % (vcf, os.path.join(outdir, "_".join([name, os.path.basename(vcf)])))
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise OSError("An error occurred while running \"%s\"" % (cmd))
    # remove temporary results
    tmpdata = os.listdir(rundir)
    for d in tmpdata:
        code = subprocess.call("rm -rf %s" % (os.path.join(rundir, d)), shell=True)
        if code != 0:
            raise OSError("An error occurred while running \"%s\"" % (cmd))


def main():
    args = parse_commandline()
    # parse input targets file
    targets = parse_targets(args.targets)
    # recover regions
    regions = compute_regions(targets, args.chrom_col, args.start_col, args.stop_col)
    # run strelka
    names = get_names(targets, args.name_col)
    for i, region in enumerate(regions):
        strelka(
            args.normal_bam, args.tumor_bam, args.genome, region, names[i], args.run_dir, args.out
        )


if __name__ == "__main__":
    main()

