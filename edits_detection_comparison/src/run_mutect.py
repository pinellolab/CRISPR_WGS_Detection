"""Run GATK Mutect2 on the set of regions specified in the input file.
"""

from glob import glob
from tqdm import tqdm

import pandas as pd

import multiprocessing
import subprocess
import argparse
import os

MUTECT2 = "gatk Mutect2"  # Mutect2
FILTER = "gatk FilterMutectCalls"  # Filter Mutect variants
PADSIZE = 10000  # pad regions by 10Kb upstream and downstream
COORDSCOL = "coords_mutect_region"  # coordinates column name


def parse_commandline():
    """The function parses the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified input regions",
        usage="\n\tpython %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
        "--bam1 <BAM> --bam2 <BAM> --normal <NORMAL-SAMPLE> --chrom-col "
        "<CHR-COL-NUM> --start-col <START-COL-NUM> --stop-col <STOP-COL-NUM> "
        "--name-col <NAME-COL-NUM> --offregion --out <OUTDIR> --threads <THREADS>",
    )
    parser.add_argument(
        "--targets", type=str, metavar="TARGETS-FILE", help="Targets coordinates file"
    )
    parser.add_argument("--genome", type=str, metavar="GENOME", help="Reference genome")
    parser.add_argument("--bam1", type=str, metavar="BAM", help="BAM file")
    parser.add_argument("--bam2", type=str, metavar="BAM", help="BAM file")
    parser.add_argument(
        "--normal", type=str, metavar="NORMAL-SAMPLE", help="Normal sample"
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
    parser.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        nargs="?",
        default=1,
        help="Number of threads to use during the run",
    )
    args = parser.parse_args()  # parse command line
    _check_args(args)  # check arguments consistency
    return args


def _check_args(args):
    """(PRIVATE)
    Check input arguments consistency

    :param args: input arguments
    :type args: argparse.Namespace
    """
    if not isinstance(args, argparse.Namespace):
        raise TypeError(
            "Expected %s, got %s"
            % (argparse.ArgumentParser.__name__, type(args).__name__)
        )
    if not isinstance(args.targets, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(args.targets).__name__)
        )
    if not os.path.isfile(args.targets):
        raise FileNotFoundError("Unable to locate %s" % (args.targets))
    if not isinstance(args.genome, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(args.genome).__name__)
        )
    if not os.path.isfile(args.genome):
        raise FileNotFoundError("Unable to locate %s" % (args.genome))
    if not isinstance(args.bam1, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(args.bam1).__name__)
        )
    if not os.path.isfile(args.bam1):
        raise FileNotFoundError("Unable to locate %s" % (args.bam1))
    if not isinstance(args.bam2, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(args.bam2).__name__)
        )
    if not os.path.isfile(args.bam2):
        raise FileNotFoundError("Unable to locate %s" % (args.bam2))
    if not isinstance(args.chrom_col, int):
        raise TypeError(
            "Expected %s, got %s" % (int.__name__, type(args.chrom_col).__name__)
        )
    if args.chrom_col < 0:
        raise ValueError("Wrong chromosome column number (%d)" % (args.chrom_col))
    if not isinstance(args.start_col, int):
        raise TypeError(
            "Expected %s, got %s" % (int.__name__, type(args.start_col).__name__)
        )
    if args.start_col < 0:
        raise ValueError("Wrong start column number (%d)" % (args.start_col))
    if not isinstance(args.stop_col, int):
        raise TypeError(
            "Expected %s, got %s" % (int.__name__, type(args.stop_col).__name__)
        )
    if args.stop_col < 0:
        raise ValueError("Wrong stop column number (%d)" % (args.stop_col))
    if not isinstance(args.name_col, int):
        raise TypeError(
            "Expected %s, got %s" % (int.__name__, type(args.name_col).__name__)
        )
    if args.name_col < 0:
        raise ValueError("Wrong name column number (%d)" % (args.name_col))


def compute_coordinates(chrom, start, stop):
    """Computes genomic coordinates string

    :param chrom: chromosome
    :type chrom: str
    :param start: start
    :type start: int
    :param stop: stop
    :type stop: int
    :raises TypeError: raise on chrom type mismatch
    :raises TypeError: raise on start type mismatch
    :raises TypeError: raise on stop type mismatch
    :raises ValueError: raise if start > stop
    :return: genomic coordinate
    :rtype: str
    """
    if not isinstance(chrom, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(chrom).__name__))
    if not isinstance(start, int):
        raise TypeError("Expected %s, got %s" % (int.__name__, type(start).__name__))
    if not isinstance(stop, int):
        raise TypeError("Expected %s, got %s" % (int.__name__, type(stop).__name__))
    if stop < start:
        raise ValueError(
            "Wrong start and stop coordinates (stop < start -- %d < %d)" % (stop, start)
        )
    return "%s:%d-%d" % (chrom, start, stop)


def parse_targets_coordinates(targets_file, chrom_col, start_col, stop_col, offregion):
    """Reda the input guide targets regions and pad them by 10kb upstream and
    downstream

    :param targets_file: guide target regions file
    :type targets_file: str
    :param chrom_col: chromosome data column
    :type chrom_col: int
    :param start_col: start data column
    :type start_col: int
    :param stop_col: stop data column
    :type stop_col: int
    :param offregion: compute off regions
    :type offregion: bool
    :raises TypeError: raise on targets_file type mismatch
    :raises FileNotFoundError: raise if targets_file cannot be found
    :raises IndexError: chromosome column index out of dataframe columns
    :raises IndexError: start column index out of dataframe columns
    :raises IndexError: stop column index out of dataframe columns
    :return: processed guide targets dataset
    :rtype: pd.DataFrame
    """
    if not isinstance(targets_file, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(targets_file).__name__)
        )
    if not os.path.isfile(targets_file):
        raise FileNotFoundError("Unable to locate %s" % (targets_file))
    targets = pd.read_csv(targets_file, sep="\t")  # load the target TSV file
    # get chromosome, start and stop column names in the dataframe
    columns = targets.columns.tolist()
    try:
        chrom_col = columns[chrom_col]
    except IndexError:
        raise IndexError(
            "The chromosome column number %d is not in %s" % (chrom_col, targets_file)
        )
    try:
        start_col = columns[start_col]
    except IndexError:
        raise IndexError(
            "The start coordinate column number %d is not in %s"
            % (start_col, targets_file)
        )
    try:
        stop_col = columns[stop_col]
    except IndexError:
        raise IndexError(
            "The stop coordinate column number %d is not in %s"
            % (stop_col, targets_file)
        )
    if offregion:  # off target regions
        targets1 = targets.copy()  # upstream shift
        targets1[stop_col] = targets1.apply(lambda x: x[start_col] - 100, axis=1)
        targets1[start_col] = targets1.apply(
            lambda x: x[start_col] - 100 - PADSIZE, axis=1
        )
        targets2 = targets.copy()  # downstream shift
        targets2[start_col] = targets2.apply(lambda x: x[stop_col] + 100, axis=1)
        targets2[stop_col] = targets2.apply(
            lambda x: x[stop_col] + 100 + PADSIZE, axis=1
        )
        targets = pd.concat([targets1, targets2])
    else:  # on target regions
        # pad start and stop coordinates
        targets[start_col] = targets.apply(lambda x: x[start_col] - PADSIZE, axis=1)
        targets[stop_col] = targets.apply(lambda x: x[stop_col] + PADSIZE, axis=1)
    # compute coordinates strings
    assert COORDSCOL not in columns
    targets[COORDSCOL] = targets.apply(
        lambda x: compute_coordinates(x[chrom_col], x[start_col], x[stop_col]), axis=1
    )
    return targets


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
    # parse the input targets file
    targets = parse_targets_coordinates(
        args.targets, args.chrom_col, args.start_col, args.stop_col, args.offregion
    )
    assert COORDSCOL in targets.columns.tolist()
    # run Mutect2 on each padded guide region
    coords = targets[COORDSCOL].tolist()
    names = targets[args.name_col].tolist()
    commands = [  # force stdout to /dev/null
        "%s -R %s -I %s -I %s -normal %s -L %s -O %s > /dev/null"
        % (
            MUTECT2,
            args.genome,
            args.bam1,
            args.bam2,
            args.normal,
            coord,
            os.path.join(args.out, "%s.%s.vcf" % (names[i], coord)),
        )
        for i, coord in enumerate(coords)
    ]
    run_commands(commands, args.threads)
    # filter variants called
    vcfs = glob(os.path.join(args.out, "*.vcf"))
    commands = [
        "%s -V %s -R %s -O %s" % (FILTER, vcf, args.genome, "%s.filtered.vcf" % (vcf))
        for vcf in vcfs
    ]
    run_commands(commands, args.threads)


if __name__ == "__main__":
    main()
