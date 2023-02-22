"""Run pindel on the set of regions specified in the input file.
"""

from tqdm import tqdm
import multiprocessing
import subprocess
import argparse
import tempfile
import os

PINDEL = "pindel"  # Pindel
PINDEL2VCF = "pindel2vcf"  # convert pindel's output to VCF files
PADSIZE = 10000  # pad regions by 10Kb


def parse_commandline():
    """The function parses the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Script to call variants in the specified regions using pindel.",
        usage="\n\tpython3 %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--normal-bam <BAM> --normal-sample <NORMAL-SAMPLE> --tumor-bam <BAM> "
        " --tumor-sample <TUMOR-SAMPLE> --chrom-col <CHROM-COL-NUM> --start-col "
        "<START-COL-NUM> --stop-col <STOP-COL-NUM> --out <OUTDIR> --threads <THREADS>",
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
    parser.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        nargs="?",
        default=1,
        help="Number of threads to use during the run",
    )
    return parser.parse_args()


def parse_targets_coordinates(targets_file, chrom_col, start_col, stop_col, offregion):
    """Parse the guide target sites file

    :param targets_file: guide target sites file
    :type targets_file: str
    :param chrom_col: chromosome data columns
    :type chrom_col: int
    :param start_col: start data columns
    :type start_col: int
    :param stop_col: stop data columns
    :type stop_col: int
    :param offregion: pad off regions
    :type offregion: bool
    :return: padded guide target sites
    :rtype: List[str]
    """
    handle = open(targets_file, mode="r")  # read the target sites file
    lines = [
        line.strip().split() for i, line in enumerate(handle) if i > 0
    ]  # skip header
    if offregion:  # pad offregions
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
        assert len(regions) == (len(lines) * 2)  # 2 times the original input (up + downstream)
    else:  # pad onregions
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
    """Recover guide target sites names  

    :param regions: padded genomic regions
    :type regions: List[str]
    :return: guide target sites names
    :rtype: List[str]
    """
    names_list = [region.replace(":", "_").replace("-", "_") for region in regions]
    assert len(names_list) == len(regions)
    return names_list


def config_file(tumor_bam, tumor_sample, normal_bam, normal_sample):
    """Write the configfile required by Pindel

    :param tumor_bam: BAM file
    :type tumor_bam: str
    :param tumor_sample: tumor sample name
    :type tumor_sample: str
    :param normal_bam: BAM file
    :type normal_bam: str
    :param normal_sample: normal sample name
    :type normal_sample: str
    :return: configfile
    :rtype: str
    """
    configfile = tempfile.NamedTemporaryFile().name
    handle = open(configfile, mode="w")
    handle.write("%s\t250\t%s\n" % (tumor_bam, tumor_sample))
    handle.write("%s\t250\t%s\n" % (normal_bam, normal_sample))
    handle.close()
    assert os.path.isfile(configfile)
    return os.path.abspath(configfile)

def run_command(command):
    """Run the command

    :param command: command
    :type command: str
    :return: command signal
    :rtype: int
    """
    return subprocess.call(command, shell=True)

def run_commands(commands, threads):
    """Run the input commands

    :param commands: input commands
    :type commands: List[str]
    :param threads: threads
    :type threads: int
    :raises OSError: raise on Mutect2 failure
    """
    pool = multiprocessing.Pool(processes=threads)  # create `threads` threads
    result = pool.map_async(run_command, commands)  # run commands in parallel
    with tqdm(total=len(commands)) as progress:
        while not result.ready():
            remaining = result._number_left
            progress.update(len(commands) - remaining)
    codes = result.get()
    for code, cmd in zip(codes, commands):
        if code != 0:
            raise OSError("An error occurred while running %s" % (cmd))

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
    commands = [
        "%s -f %s -i %s -c %s -o %s" % (
            PINDEL, args.genome, configfile, region, os.path.join(args.out, names[i])
        )
        for i, region in enumerate(regions)
    ]
    run_commands(commands, args.threads)
    # convert pindel files in VCFs (deletions)
    commands = [
        "%s -p %s_D -r %s -R hg38 -d 20131217 -G" % (
            PINDEL2VCF, os.path.join(args.out, name), args.genome
        )
        for name in names
    ]
    run_commands(commands, args.threads)
    # convert pindel files in VCFs (insertions)
    commands = [
        "%s -p %s_SI -r %s -R hg38 -d 20131217 -G" % (
            PINDEL2VCF, os.path.join(args.out, name), args.genome
        )
        for name in names
    ]
    run_commands(commands, args.threads)
    cmd = "rm %s" % (configfile)  # delete the config file
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error ocuurred while running %s" % (cmd))


if __name__ == "__main__":
    main()
