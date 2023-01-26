"""Run pindel on the set of regions specified in the input file.
"""

import subprocess
import argparse
import tempfile
import sys
import os

PINDEL = "pindel"
PINDEL2VCF = "pindel2vcf"
PADSIZE = 10000


def parse_commandline():
    """The function parses the command line arguments."""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Script to call variants in the specified regions using pindel.",
        usage="\n\tpython3 %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--normal-bam <BAM> --normal-sample <NORMAL-SAMPLE> --tumor-bam <BAM> "
        " --tumor-sample <TUMOR-SAMPLE> --chrom-col <CHROM-COL-NUM> --start-col "
        "<START-COL-NUM> --stop-col <STOP-COL-NUM> --out <OUTDIR>",
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
        "--normal-sample",
        type=str,
        metavar="NORMAL-SAMPLE",
        dest="normal_sample",
        help="Normal sample",
    )
    parser.add_argument(
        "--tumor-bam", type=str, metavar="BAM", dest="tumor_bam", help="Tumor BAM file"
    )
    parser.add_argument(
        "--tumor-sample",
        type=str,
        metavar="TUMOR-SAMPLE",
        dest="tumor_sample",
        help="Tumor sample",
    )
    parser.add_argument(
        "--chrom-col",
        type=int,
        metavar="CHROM-COL-NUM",
        help="Column containing chromosome",
        dest="chrom_col",
    )
    parser.add_argument(
        "--start-col",
        type=int,
        metavar="START-COL-NUM",
        help="Column cotaining the start coordinates",
        dest="start_col",
    )
    parser.add_argument(
        "--stop-col",
        type=int,
        metavar="STOP-COL-NUM",
        help="Column containing the stop coordinates",
        dest="stop_col",
    )
    parser.add_argument(
        "--offregion",
        action="store_true",
        default=False,
        help="Shift the target regions 100bp upstream and downstream (not overlapping "
        "the original site)",
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    args = parser.parse_args()
    return args


def parse_targets_coordinates(targets_file, chrom_col, start_col, stop_col, offregion):
    """The function parses the target sites file."""

    handle = open(targets_file, mode="r")
    lines = [
        line.strip().split() for i, line in enumerate(handle) if i > 0
    ]  # skip header
    if offregion:
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
    """The function returns a list with the target names."""

    names_list = [region.replace(":", "_").replace("-", "_") for region in regions]
    assert len(names_list) == len(regions)
    return names_list


def config_file(tumor_bam, tumor_sample, normal_bam, normal_sample):
    """The function writes the config file required by pindel."""

    configfile = tempfile.NamedTemporaryFile().name
    handle = open(configfile, mode="w")
    handle.write("%s\t250\t%s\n" % (tumor_bam, tumor_sample))
    handle.write("%s\t250\t%s\n" % (normal_bam, normal_sample))
    handle.close()
    assert os.path.isfile(configfile)
    return os.path.abspath(configfile)


def pindel(configfile, genome, region, outdir, name):
    """The function runs pindel on the input regions."""

    # call variants
    outfile = os.path.join(outdir, name)
    cmd = "%s -f %s -i %s -c %s -o %s" % (PINDEL, genome, configfile, region, outfile)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError('An error occurred while running "%s"' % (cmd))
    # convert pindel output in VCF format
    cmd = "%s -p %s_D -r %s -R hg38 -d 20131217 -G" % (PINDEL2VCF, outfile, genome)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError('An error occurred while running "%s"' % (cmd))
    cmd = "%s -p %s_SI -r %s -R hg38 -d 20131217 -G" % (PINDEL2VCF, outfile, genome)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError('An error occurred while running "%s"' % (cmd))


def main():
    # parse commandline
    args = parse_commandline()
    # parse the input target coordinates
    regions = parse_targets_coordinates(
        args.targets, args.chrom_col, args.start_col, args.stop_col, args.offregion
    )
    # get target names
    names = get_names(regions)
    # create the config file to run pindel on the current data
    configfile = config_file(
        args.tumor_bam, args.tumor_sample, args.normal_bam, args.normal_sample
    )
    # run pindel
    for i, region in enumerate(regions):
        pindel(configfile, args.genome, region, args.out, names[i])
    cmd = "rm %s" % (configfile)  # delete the config file
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error ocuurred while running %s" % (cmd))


if __name__ == "__main__":
    main()
