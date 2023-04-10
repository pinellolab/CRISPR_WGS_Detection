"""Target files parsers. The parsers recover the offtarget regions identified
by CasOffinder and add 10bp upstream and downstream to each region. The padded 
regions are stored in the format `chr`:`start`-`stop`.
"""

def read_targets_file(targets_file):
    try:
        with open(targets_file, mode="r") as infile:
            lines = [line.strip().split() for line in infile]
    except OSError as e:
        raise OSError("Errors occurred while parsing the targets file")
    return lines

def get_coordinate(target):
    chrom, pos = target[1:3] # chrom and start on column 2 and 3 respectively
    strand = target[4]  # strand on fourth column
    center = int(pos) + 17 if strand == "+" else int(pos) + 6
    start = center - 10
    stop = center + 10
    return "%s:%s-%s" % (chrom, start, stop)

def compute_coordinates(targets):
    return [get_coordinate(target) for target in targets]

def parse_targets_coordinates(targets_file):
    lines = read_targets_file(targets_file)  # read targets file
    return compute_coordinates(lines)  # compute target regions
