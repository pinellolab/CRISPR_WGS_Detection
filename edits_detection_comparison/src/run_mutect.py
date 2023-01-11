"""Run GATK Mutect2 on the set of regions specified in the input file.
"""

from typing import Optional
from glob import glob
from tqdm import tqdm

import pandas as pd
import numpy as np

import subprocess
import argparse
import sys
import os

MUTECT2 = "gatk Mutect2"  # Mutect2
FILTER = "gatk FilterMutectCalls"  # Filter Mutect variants
PADSIZE = 10000  # pad regions by 10Kb upstream and downstream
COORDSCOL = "coords_mutect_region"  # coordinates column name

def parse_commandline() -> argparse.ArgumentParser:
    """The function parses the command line arguments.

    ...
    
    Parameters
    ----------
    None

    Returns
    -------
    argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified input regions",
        usage="\n\tpython3 %(prog)s --targets <TARGETS-FILE> --genome <GENOME> --bam1 <BAM> --bam2 <BAM> --normal <NORMAL-SAMPLE> --chrom-col <CHR-COL-NUM> --start-col <START-COL-NUM> --stop-col <STOP-COL-NUM> --name-col <NAME-COL-NUM> --offregion --out <OUTDIR>",
    )
    parser.add_argument("--targets", type=str, metavar="TARGETS-FILE", help="Targets coordinates file")
    parser.add_argument("--genome", type=str, metavar="GENOME", help="Reference genome")
    parser.add_argument("--bam1", type=str, metavar="BAM", help="BAM file")
    parser.add_argument("--bam2", type=str, metavar="BAM", help="BAM file")
    parser.add_argument("--normal", type=str, metavar="NORMAL-SAMPLE", help="Normal sample")
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
    parser.add_argument(
        "--offregion", 
        action="store_true", 
        default=False, 
        help="Shift the target regions 100bp upstream and downstream (not "
             "overlapping the original site)"
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    args = parser.parse_args()
    # check arguments consistency
    __check_args(args)
    return args


def __check_args(args: argparse.Namespace) -> None:
    """ (PRIVATE)
    
    The function checks input arguments consistency.
    
    ...
    
    Parameters
    ----------
    args
        Input arguments

    Returns
    -------
    None
    """

    if not isinstance(args, argparse.Namespace):
        raise TypeError(f"Expected {argparse.ArgumentParser.__name__}, got {type(args).__name__}")
    if not isinstance(args.targets, str):
        raise TypeError(f"Expected {str.__name__}, got {type(args.targets).__name__}")
    if not os.path.isfile(args.targets):
        raise FileNotFoundError(f"Unable to locate {args.targets}")
    if not isinstance(args.genome, str):
        raise TypeError(f"Expected {str.__name__}, got {type(args.genome).__name__}")
    if not os.path.isfile(args.genome):
        raise FileNotFoundError(f"Unable to locate {args.genome}")
    if not isinstance(args.bam1, str):
        raise TypeError(f"Expected {str.__name__}, got {type(args.bam1).__name__}")
    if not os.path.isfile(args.bam1):
        raise FileNotFoundError(f"Unable to locate {args.bam1}")
    if not isinstance(args.bam2, str):
        raise TypeError(f"Expected {str.__name__}, got {type(args.bam2).__name__}")
    if not os.path.isfile(args.bam2):
        raise FileNotFoundError(f"Unable to locate {args.bam2}")
    if not isinstance(args.chrom_col, int):
        raise TypeError(f"Expected {int.__name__}, got {type(args.chrom_col).__name__}")
    if args.chrom_col < 0:
        raise ValueError(f"Wrong chromosome column number ({args.chrom_col})")
    if not isinstance(args.start_col, int):
        raise TypeError(f"Expected {int.__name__}, got {type(args.start_col).__name__}")
    if args.start_col < 0:
        raise ValueError(f"Wrong start column number ({args.start_col})")
    if not isinstance(args.stop_col, int):
        raise TypeError(f"Expected {int.__name__}, got {type(args.stop_col).__name__}")
    if args.stop_col < 0:
        raise ValueError(f"Wrong stop column number ({args.stop_col})")
    if not isinstance(args.name_col, int):
        raise TypeError(f"Expected {int.__name__}, got {type(args.name_col).__name__}")
    if args.name_col < 0:
        raise ValueError(f"Wrong name column number ({args.name_col})")


def compute_coordinates(chrom: str, start: int, stop: int) -> str:
    """The function computes the coordinate strings (format: <chrom:start-stop>
    
    ...

    Parameters
    ----------
    chrom
        Chromosome
    start
        Start position
    stop
        Stop position

    Returns
    -------
    str
    """

    if not isinstance(chrom, str):
        raise TypeError(f"Expected {str.__name__}, got {type(chrom).__name__}")
    if not isinstance(start, int):
        raise TypeError(f"Expected {int.__name__}, got {type(start).__name__}")
    if not isinstance(stop, int):
        raise TypeError(f"Expected {int.__name__}, got {type(stop).__name__}")
    if stop < start:
        raise ValueError(
            f"Wrong start and stop coordinates (stop < start -- {stop} < {start})"
        )
    coord = f"{chrom}:{start}-{stop}"
    return coord


def parse_targets_coordinates(targets_file: str, chrom_col: int, start_col: int, stop_col: int, offregion: bool) -> pd.DataFrame:
    """The function reads the input guide sequences file. The function pad the 
    guide coordinates by 10 Kbs upstream and downstream.

    ...

    Parameters
    ----------
    targets_file
        Input guide sequences file
    chrom_col
        Chromosome column number
    start_col
        Start coordinate columns number
    stop_col
        Stop coordinate column number
    offregion
        Shift the start and stop coordinates by 100 bp upstream and downstream 
        with respect to the original target site

    Returns
    -------
    pd.DataFrame
    """

    if not isinstance(targets_file, str):
        raise TypeError(f"Expected {str.__name__}, got {type(targets_file).__name__}")
    if not os.path.isfile(targets_file):
        raise FileNotFoundError(f"Unable to locate {targets_file}")
    targets = pd.read_csv(targets_file, sep="\t")  # load the target TSV file
    # get chromosome, start and stop column names in the dataframe
    columns = targets.columns.tolist()
    try:
        chrom_col = columns[chrom_col]
    except IndexError:
        raise IndexError(f"The chromosome column number {chrom_col} is not in {targets_file}")
    try:
        start_col = columns[start_col]
    except IndexError:
        raise IndexError(f"The start coordinate column number {start_col} is not in {targets_file}")
    try:
        stop_col = columns[stop_col]
    except IndexError:
        raise IndexError(f"The stop coordinate column number {stop_col} is not in {targets_file}")
    if offregion:  # off target regions
        targets1 = targets.copy()  # upstream shift
        targets1[stop_col] = targets1.apply(lambda x : x[start_col] - 100, axis=1)
        targets1[start_col] = targets1.apply(lambda x : x[start_col] - 100 - PADSIZE, axis=1)
        targets2 = targets.copy()  # downstream shift
        targets2[start_col] = targets2.apply(lambda x : x[stop_col] + 100, axis=1)
        targets2[stop_col] = targets2.apply(lambda x : x[stop_col] + 100 + PADSIZE, axis=1)
        targets = pd.concat([targets1, targets2])
    else:  # on target regions
        # pad start and stop coordinates
        targets[start_col] = targets.apply(lambda x : x[start_col] - PADSIZE, axis=1)
        targets[stop_col] = targets.apply(lambda x : x[stop_col] + PADSIZE, axis=1)
    # compute coordinates strings
    assert COORDSCOL not in columns
    targets[COORDSCOL] = targets.apply(lambda x : compute_coordinates(x[chrom_col], x[start_col], x[stop_col]), axis=1)
    return targets


def mutect2(
    genome: str, bam1: str, bam2: str, normal_sample: str, region: str, name: str, outdir: str
) -> None:
    """The function runs Mutect2 on the previously computed regions.

    ...

    Parameters
    ----------
    genome
        Genome
    bam1
        First BAM file
    bam2
        Second BAM file
    normal_sample
        Normal sample name
    region
        Genomic region
    name
        Output VCF filename
    outdir
        Output directory

    Returns
    -------
    None
    """

    outfile = os.path.join(outdir, f"{name}.{region}.vcf")
    cmd = f"{MUTECT2} -R {genome} -I {bam1} -I {bam2} -normal {normal_sample} -L {region} -O {outfile}"
    code = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if code != 0:
        raise subprocess.SubprocessError(f"An error occurred while running \"{cmd}\"")


def filter_mutect(vcf: str, genome: str, vcfout: str) -> None:
    """The function filters the mutations called by Mutect2 during the previous 
    step.

    ...

    Parameters
    ----------
    vcf
        VCF file
    genome
        Genome
    vcfout
        Ouput VCF filename

    Returns
    -------
    None
    """

    cmd = f"{FILTER} -V {vcf} -R {genome} -O {vcfout}"
    code = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if code != 0:
        raise subprocess.SubprocessError(f"An error occurred while running \"{cmd}\"")


def main():
    # parse commandline
    args = parse_commandline()
    # parse the input targets file
    targets = parse_targets_coordinates(args.targets, args.chrom_col, args.start_col, args.stop_col, args.offregion)
    assert COORDSCOL in targets.columns.tolist()
    # run Mutect2 on each padded guide region
    tqdm.pandas()
    targets.progress_apply(
        lambda x : mutect2(args.genome, args.bam1, args.bam2, args.normal, x[COORDSCOL], x[args.name_col], args.out),
        axis=1
    )
    # filter variants called 
    vcfs = glob(os.path.join(args.out, "*.vcf"))
    for vcf in tqdm(vcfs):
        vcfout = f"{vcf}.filtered.vcf"
        filter_mutect(vcf, args.genome, vcfout)

if __name__ == "__main__":
    main()
