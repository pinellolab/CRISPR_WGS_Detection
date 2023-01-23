"""Run Varscan on the set of regions specified in the input file.
"""

import subprocess
import tempfile
import argparse
import os


SAMTOOLS = "samtools mpileup"
SAMTOOLSTMP = tempfile.mkdtemp()  # stores temporary MPILEUP files
VARSCAN = "varscan somatic"
PADSIZE = 10000


def parse_commandline():
    """The function parses the command line arguments."""
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to call variants in the specified input regions",
        usage="\n\tpython3 %(prog)s --targets <TARGETS-FILE> --genome <GENOME> "
              "--normal-bam <BAM> --tumor-bam <BAM> --chrom-col <CHROM-COL-NUM> "
              "--start-col <START-COL-NUM> --stop-col <STOP-COL-NUM> --offregion "
              "--out <OUTDIR>"
    )
    parser.add_argument(
        "--targets", type=str, metavar="TARGETS-FILE", help="Targets coordinates file"
    )
    parser.add_argument("--genome", type=str, metavar="GENOME", help="Reference genome")
    parser.add_argument(
        "--normal-bam", type=str, metavar="BAM", help="Normal BAM file", dest="normal_bam"
    )
    parser.add_argument(
        "--tumor-bam", type=str, metavar="BAM", help="Tumor BAM file", dest="tumor_bam"
    )
    parser.add_argument(
        "--chrom-col", type=int, metavar="CHROM-COL-NUM", help="Column containing chromosome", dest="chrom_col"
    )
    parser.add_argument(
        "--start-col", type=int, metavar="START-COL-NUM", help="Column containing start coordinates", dest="start_col"
    )
    parser.add_argument(
        "--stop-col", type=int, metavar="STOP-COL-NUM", help="Column containing stop coordinates", dest="stop_col"
    )
    parser.add_argument(
        "--offregion",
        action="store_true", 
        default=False, 
        help="Shift the target regions 100bp upstream and dowstream (not "
             "overlapping the original site)"
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    args = parser.parse_args()
    return args


def parse_targets(targets):
    """The function parses the target sites file."""

    try:
        handle = open(targets, mode="r")
        lines = [line.strip().split() for i, line in enumerate(handle) if i > 0]  # skip header 
    except OSError:
        raise OSError("An error occurred while reading %s" % (targets))
    finally:
        handle.close()
    return lines


def get_names(regions):
    """The function assign a name to each input region."""
    
    names_list = [region.replace(":", "_").replace("-", "_") for region in regions]
    assert len(names_list) == len(regions)
    return names_list


def compute_regions(lines, chrom_col, start_col, stop_col, offregion):
    """The function computes the padded genomic regions."""

    if offregion:
        regions = [
            "%s:%s-%s" % (
                line[chrom_col], (int(line[start_col]) - 100 - PADSIZE), (int(line[start_col]) - 100)
            )
            for line in lines
        ] + [
            "%s:%s-%s" % (
                line[chrom_col], (int(line[stop_col]) + 100), (int(line[stop_col]) + 100 + PADSIZE)
            )
            for line in lines
        ]
        assert len(regions) == (len(lines) * 2)
    else:
        regions = [
            "%s:%s-%s" % (
                line[chrom_col], (int(line[start_col]) - PADSIZE), (int(line[stop_col]) + PADSIZE)
            )
            for line in lines
        ]
        assert len(lines) == len(regions)
    return regions


def mpileup(bam, genome, region, outfile):
    """The function computes the MPILEUP file using SAMtools for the input BAM 
    file. The MPILEUP file covers the input region.
    """

    cmd = "%s -f %s -d %d -r %s %s > %s" % (SAMTOOLS, genome, 0, region, bam, outfile)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error occurred while running \"%s\"" % (cmd))
    return outfile


def varscan(genome, normal_bam, tumor_bam, region, outdir, name):
    """The function call variants on the input region using varscan."""

    # compute the mpileup files for normal and tumor BAMs
    mpileup_normal = mpileup(
        normal_bam, genome, region, os.path.join(SAMTOOLSTMP, "%s_normal.mpileup" % (name))
    )
    mpileup_tumor = mpileup(
        tumor_bam, genome, region, os.path.join(SAMTOOLSTMP, "%s_tumor.mpileup" % (name))
    )
    # call varscan
    outfile = os.path.join(outdir, name)
    cmd = "%s %s %s %s --min-var-freq 0.0001" % (
        VARSCAN, mpileup_normal, mpileup_tumor, outfile
    )
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error occurred while running %s" % (cmd))
    # delete the mpileup files
    cmd = "rm %s %s" % (mpileup_normal, mpileup_tumor)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("An error occurred while running %s" % (cmd))
    


def main():
    args = parse_commandline()
    # parse the input targets
    targets = parse_targets(args.targets)
    # recover the target genomic regions
    regions = compute_regions(targets, args.chrom_col, args.start_col, args.stop_col, args.offregion)
    # run varscan
    names = get_names(regions)
    assert os.path.exists(SAMTOOLSTMP)
    for i, region in enumerate(regions):
        varscan(args.genome, args.normal_bam, args.tumor_bam, region, args.out, names[i])
    # delete the samtools temporary directory storing the mpileup files
    cmd = "rm -rf %s" % (SAMTOOLSTMP)
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError("AN error occurred while running %s" % (cmd))
    assert not os.path.exists(SAMTOOLSTMP)


if __name__ == "__main__":
    main()

