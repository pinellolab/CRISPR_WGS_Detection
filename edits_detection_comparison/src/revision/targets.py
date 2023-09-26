"""Functions to read and process target sites coordinates
"""

from utils import ONTARGETWIDTH, OFFTARGETWIDTH

import os

def read_targets_casoffinder(targets_file):
    """Read CasOffinder target file

    :param targets_file: targets
    :type targets_file: str
    :raises FileNotFoundError: raise if targets_file not found
    :raises OSError: raise if file read fails
    :return: targets
    :rtype: List[str]
    """
    if not os.path.isfile(targets_file):
        raise FileNotFoundError("%s not found!" % (targets_file))
    try:
        with open(targets_file, mode="r") as infile:
            infile.readline()  # skip header
            lines = [line.strip().split() for line in infile]
        assert lines
    except OSError as e:
        raise OSError("Reading %s failed!" % (targets_file))
    return lines

def process_coordinates_casoffinder(target, upstream, downstream):
    """Define the on-target and off-target (upstream and downstream) site 
    coordinates

    :param target: target
    :type target: List[str]
    :param upstream: upstream region
    :type upstream: bool
    :param downstream: downstream region
    :type downstream: bool
    :raises ValueError: raise if both conditions are True
    :return: padded target region
    :rtype: str
    """
    if sum([upstream, downstream]) not in [0, 1]:
        raise ValueError("Multiple editing events calling locations!")
    chrom, pos = target[3:5]  # chrom and target start
    strand = target[5]  # direction
    center = int(pos) + 17 if strand == "+" else int(pos) + 6  # target site center
    start = center - (ONTARGETWIDTH / 2)  # 100 bp region around the target site
    stop = center + (ONTARGETWIDTH / 2)
    if upstream:
        stop = start
        start = stop - OFFTARGETWIDTH
    elif downstream:
        start = stop
        stop = start + OFFTARGETWIDTH
    return "%s:%s-%s" % (chrom, start, stop)
    
    

def read_targets_file(targets_file, casoffinder, circleseq, guideseq, upstream, downstream):
    if sum([casoffinder, circleseq, guideseq]) != 1:
        raise ValueError("Multiple input targets data sources!")
    if casoffinder:
        return [
            process_coordinates_casoffinder(target, upstream, downstream)
            for target in read_targets_casoffinder(targets_file)
        ]

