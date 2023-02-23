"""The script constructs the reports storing the edits called by the different 
tools considered in the study. The reports are stored as TSV files.
"""

from tqdm import tqdm

import pandas as pd

import argparse
import os

# data directories
BASEDIR = "../data"
GUIDESEQ = os.path.join(BASEDIR, "offtarget_detection/guideseq")
CIRCLESEQ = os.path.join(BASEDIR, "offtarget_detection/circleseq/")
EDITS = os.path.join(BASEDIR, "VCFs")
REPORTS = os.path.join(BASEDIR, "edits_reports")
# variant calling tools
VCALLINGTOOLS = ["mutect2", "strelka", "pindel", "varscan"]
# guides, experiments type and cell types
GUIDES = ["EMX1", "HEKSite4", "RNF2", "VEGFASite3"]
EXPERIMENTS = ["circleseq", "guideseq"]
CELLTYPES = ["GM12878", "K562"]
# genomic regions padding size
PADSIZE = 10000
# variants type
INSERTION = "insertion"
DELETION = "deletion"
SNV = "snv"
# report columns
REPORT_COLS = [
    "SITE",
    "CHROM",
    "VAR-POS",
    "REF",
    "ALT",
    "FILTER",
    "DP-NORMAL",
    "DP-TUMOR",
    "AD-NORMAL-REF",
    "AD-NORMAL-ALT",
    "AD-TUMOR-REF",
    "AD-TUMOR-ALT",
    "TYPE",
    "TARGET-START",
    "TARGET-STOP",
    "STRAND",
    "MISMATCHES",
]
CIRCLESEQ_COLS = [
    "SITE",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER",
    "DP-NORMAL",
    "DP-TUMOR",
    "AD-NORMAL-REF",
    "AD-NORMAL-ALT",
    "AD-TUMOR-REF",
    "AD-TUMOR-ALT",
    "TYPE",
    "Start",
    "End",
    "Strand",
    "Distance",
]
GUIDESEQ_COLS = [
    "SITE",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER",
    "DP-NORMAL",
    "DP-TUMOR",
    "AD-NORMAL-REF",
    "AD-NORMAL-ALT",
    "AD-TUMOR-REF",
    "AD-TUMOR-ALT",
    "TYPE",
    "start",
    "end",
    "Strand",
    "Mismatch Total",
]


def parse_commandline():
    """The function parses the command line arguments provided as input

    :return: parsed input arguments
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Script to construct the reports storing the called edits",
        usage="\n\tpython3 %(prog)s --tool <TOOL-NAME> --type <TYPE> --out "
        "<OUTDIR> --offregion",
    )
    parser.add_argument(
        "--tool",
        type=str,
        metavar="TOOL-NAME",
        help="Variant calling tool. Available values: <mutect2, strelka, "
        "varscan, pindel>",
    )
    parser.add_argument(
        "--type",
        type=str,
        metavar="TYPE",
        help="Target sites validation experiment type. Available values: "
        "<guideseq, circleseq>",
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    parser.add_argument(
        "--offregion",
        action="store_true",
        default=False,
        help="Shift the target regions 100bp upstream and downstream (not "
        "overlapping the original target site)",
    )
    return parser.parse_args()


def read_targets(exp_type, guide):
    """Parse the guide target sites file

    :param exp_type: experiment type
    :type exp_type: str
    :param guide: guide
    :type guide: str
    :raises TypeError: raise on exp_type type mismatch
    :raises ValueError: raise on forbidden exp_type value
    :raises TypeError: raise on guide type mismatch
    :raises ValueError: raise on forbidden guide value
    :return: guide target sites
    :rtype: pd.DataFrame
    """

    if not isinstance(exp_type, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(exp_type).__name__))
    if exp_type not in EXPERIMENTS:
        raise ValueError("Forbidden experiment type (%s)" % (exp_type))
    if not isinstance(guide, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(guide).__name__))
    if not guide in GUIDES:
        raise ValueError("Forbidden guide (%s)" % (guide))
    if exp_type == EXPERIMENTS[0]:  # circleseq
        targets = pd.read_csv(
            os.path.join(CIRCLESEQ, "%s.circleseq.hg19.hg38.targetname" % (guide)),
            sep="\t",
        )
        cols = targets.columns.tolist()  # rename columns for later join
        targets.columns = cols[:-1] + ["SITE"]
    else:  # guideseq
        targets = pd.read_csv(os.path.join(GUIDESEQ, "%s.guideseq" % (guide)), sep="\t")
        cols = targets.columns.tolist()  # rename columns for later join
        targets.columns = cols[:6] + ["SITE"] + cols[7:]
    assert not targets.empty
    return targets


def read_vcf(vcf):
    """Parse the input VCF

    :param vcf: VCF file
    :type vcf: str
    :raises TypeError: raise on vcf type mismatch
    :raises FileNotFoundError: raise if vcf cannot be found
    :raises OSError: raise on read() failure
    :return: VCF lines
    :rtype: List[List]
    """
    if not isinstance(vcf, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(vcf).__name__))
    if not os.path.isfile(vcf):
        raise FileNotFoundError("Unable to locate %s" % (vcf))
    try:
        handle = open(vcf, mode="r")  # parse the VCF
        variants = [line.strip().split() for line in handle if not line.startswith("#")]
    except OSError:
        raise OSError("An error occurred while reading %s" % (vcf))
    finally:
        handle.close()
    return variants


def _recover_depth(edit, depth, normal, reference, pindel, varscan):
    """(PRIVATE)
    Recover read depth supporting the reference and alternative alleles

    :param edit: edit
    :type edit: List
    :param depth: edit read depth
    :type depth: str
    :param normal: normal sample, defaults to True
    :type normal: Optional[bool], optional
    :param reference: reference, defaults to True
    :type reference: Optional[bool], optional
    :param pindel: from pindel, defaults to False
    :type pindel: Optional[bool], optional
    :param varscan: from varscan, defaults to False
    :type varscan: Optional[bool], optional
    :raises TypeError: raise on depth type mismatch
    :raises ValueError: raise on forbidden depth value
    :return: edit read depth
    :rtype: int
    """
    if not isinstance(depth, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(depth).__name__))
    if depth not in ["AD", "DP"]:
        raise ValueError("Forbidden read depth read (%s)" % (depth))
    assert len(edit) == 11  # each edit contains exactly 11 fields
    depthidx = edit[8].split(":").index(depth)
    if depth == "DP":  # read depth
        if normal:  # normal condition
            if pindel:  # pindel's VCFs miss DP field
                depthidx = edit[8].split(":").index("AD")
                return sum(map(int, edit[9].split(":")[depthidx].split(",")))
            return int(edit[9].split(":")[depthidx])
        if pindel:  # pindel's VCFs miss DP field
            depthidx = edit[8].split(":").index("AD")
            return sum(map(int, edit[10].split(":")[depthidx].split(",")))
        return int(edit[10].split(":")[depthidx])  # tumor condition
    else:  # allele depth
        if normal:  # normal condition
            if reference:  # reads supporting ref allele
                if varscan:  # varscan has explicit allele read count
                    depthidx = edit[8].split(":").index("RD")
                    return int(edit[9].split(":")[depthidx])
                return int(edit[9].split(":")[depthidx].split(",")[0])
            # reads supporting alt allele
            if varscan:  # varscan has explicit allele read count
                depthidx = edit[8].split(":").index("AD")
                return int(edit[9].split(":")[depthidx])
            return int(edit[9].split(":")[depthidx].split(",")[1])
        else:  # tumor condition
            if reference:  # reads supporting ref allele
                if varscan:  # varscan has explicit allele read count
                    depthidx = edit[8].split(":").index("RD")
                    return int(edit[10].split(":")[depthidx])
                return int(edit[10].split(":")[depthidx].split(",")[0])
            # reads supporting alt allele
            if varscan:  # varscan has explicit allele read count
                depthidx = edit[8].split(":").index("AD")
                return int(edit[10].split(":")[depthidx])
            return int(edit[10].split(":")[depthidx].split(",")[1])


def edits_dataframe(edits, target_sites, strelka=False, pindel=False, varscan=False):
    """Construct dataframe from input variants

    :param edits: variants list
    :type edits: List
    :param target_sites: guide traget sites
    :type target_sites: List
    :param strelka: variants called by strelka, defaults to False
    :type strelka: Optional[bool], optional
    :param pindel: variants called by pindel, defaults to False
    :type pindel: Optional[bool], optional
    :param varscan: variants called by varscan, defaults to False
    :type varscan: Optional[bool], optional
    :return: edits dataset stored as dataframe
    :rtype: pd.DataFrame
    """
    assert len(edits) == len(target_sites)
    # edits dictionary
    edf = {
        "SITE": [],
        "CHROM": [],
        "POS": [],
        "ID": [],
        "REF": [],
        "ALT": [],
        "FILTER": [],
        "DP-NORMAL": [],
        "DP-TUMOR": [],
        "AD-NORMAL-REF": [],
        "AD-NORMAL-ALT": [],
        "AD-TUMOR-REF": [],
        "AD-TUMOR-ALT": [],
    }
    for i, tse in enumerate(edits):
        if bool(tse):  # skip target sites without called edits
            for e in tse:
                edf["SITE"].append(target_sites[i])
                edf["CHROM"].append(e[0])
                edf["POS"].append(e[1])
                edf["ID"].append(e[2])
                edf["REF"].append(e[3])
                edf["ALT"].append(e[4])
                edf["FILTER"].append(e[6])
                edf["DP-NORMAL"].append(
                    _recover_depth(e, "DP", normal=True, pindel=pindel)
                )
                edf["DP-TUMOR"].append(
                    _recover_depth(e, "DP", normal=False, pindel=pindel)
                )
                if strelka:  # strelka does not report AD
                    edf["AD-NORMAL-REF"].append(0)
                    edf["AD-NORMAL-ALT"].append(0)
                    edf["AD-TUMOR-REF"].append(0)
                    edf["AD-TUMOR-ALT"].append(0)
                else:
                    edf["AD-NORMAL-REF"].append(
                        _recover_depth(
                            e, "AD", normal=True, reference=True, varscan=varscan
                        )
                    )
                    edf["AD-NORMAL-ALT"].append(
                        _recover_depth(
                            e, "AD", normal=True, reference=False, varscan=varscan
                        )
                    )
                    edf["AD-TUMOR-REF"].append(
                        _recover_depth(
                            e, "AD", normal=False, reference=True, varscan=varscan
                        )
                    )
                    edf["AD-TUMOR-ALT"].append(
                        _recover_depth(
                            e, "AD", normal=False, reference=False, varscan=varscan
                        )
                    )
    edf = pd.DataFrame(edf)  # build dataframe from dict
    if edf.empty and any([bool(tse) for tse in edits]):
        raise ValueError(f"DataFrame is empty, but some edits have been called")
    return edf


def _etype(ref, alt):
    """(PRIVATE)
    Assign variant type

    :param ref: reference allele
    :type ref: str
    :param alt: alternative allele
    :type alt: str
    :return: edits type
    :rtype: str
    """
    if len(ref) < len(alt):
        return INSERTION
    elif len(alt) < len(ref):
        return DELETION
    return SNV


def assign_etype(ref, alt):
    """Assign type to inpu edit <insertion, deletion, snv>

    :param ref: reference allele
    :type ref: str
    :param alt: alternative allele
    :type alt: str
    :raises TypeError: raise on ref type mismatch
    :raises TypeError: raise on alt type mismatch
    :return: edit type
    :rtype: str
    """
    if not isinstance(ref, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(ref).__name__))
    if not isinstance(alt, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(alt).__name__))
    if "," in alt:  # polyploid alternative allele
        etypes = [_etype(ref, aa) for aa in alt.split(",")]
        return "-".join(list(set(etypes)))
    else:  # regular alternative allele
        return _etype(ref, alt)


def compute_distance(epos, target_pos):
    """Compute the distance between the edit site and their target site

    :param epos: edit position
    :type epos: int
    :param target_pos: target site position
    :type target_pos: int
    :raises TypeError: raise on epos type mismatch
    :raises ValueError: raise on epos forbidden value
    :raises TypeError: raise on target_pos type mismatch
    :raises ValueError: raise on target_pos forbidden value
    :return: distance
    :rtype: int
    """
    if not isinstance(epos, int):
        raise TypeError("Expected %s, got %s" % (int.__name__, type(epos).__name__))
    if epos < 0:
        raise ValueError("Forbidden edit position (%d)" % (epos))
    if not isinstance(target_pos, int):
        raise TypeError(
            "Expected %s, got %s" % (int.__name__, type(target_pos).__name__)
        )
    if target_pos < 0:
        raise ValueError("Forbidden edit position (%d)" % (target_pos))
    return epos - target_pos


def edit_location(epos, target_start, target_stop, ref, alt, etype):
    """Assess if the variant occurred within the target site

    :param epos: _description_
    :type epos: int
    :param target_start: _description_
    :type target_start: int
    :param target_stop: _description_
    :type target_stop: int
    :param ref: _description_
    :type ref: str
    :param alt: _description_
    :type alt: str
    :param etype: _description_
    :type etype: str
    :return: _description_
    :rtype: str
    """
    if epos >= target_start and epos <= target_stop:
        return "TP"
    else:
        if "," in alt:  # polyploid alternative allele
            if INSERTION in etype or DELETION in etype:
                for aa in alt.split(","):
                    if len(aa) < len(ref):  # deletion:
                        padpos = list(range(epos, epos + len(ref)))
                        if any(
                            [p <= target_stop and p >= target_start for p in padpos]
                        ):
                            return "TP"
                    elif len(ref) < len(aa):  # insertion
                        padpos = list(range(epos, epos + len(aa)))
                        if any(
                            [p <= target_stop and p >= target_start for p in padpos]
                        ):
                            return "TP"
        else:
            if DELETION in etype:
                padpos = list(range(epos, epos + len(ref)))
                if any([p <= target_stop and p >= target_start for p in padpos]):
                    return "TP"
            elif INSERTION in etype:
                padpos = list(range(epos, epos + len(alt)))
                if any([p <= target_stop and p >= target_start for p in padpos]):
                    return "TP"
    return "FP"


def process_edits_dataframe(edits, targets, exp_type) -> pd.DataFrame:
    """Refine edits report. Assign to each edit the type, the distance form
    start and stop coordinates, and a flag stating if the edit occurs inside or
    outside the target site

    :param edits: edits dataset
    :type edits: pd.DataFrame
    :param targets: targets dataset
    :type targets: pd.DataFrame
    :param exp_type: experiment type
    :type exp_type: str
    :raises TypeError: raise on edits type mismatch
    :raises TypeError: raise on targets type mismatch
    :return: refined edits dataset
    :rtype: pd.DataFrame
    """
    if not isinstance(edits, pd.DataFrame):
        raise TypeError(
            "Expected %s, got %s" % (pd.DataFrame.__name__, type(edits).__name__)
        )
    if not isinstance(targets, pd.DataFrame):
        raise TypeError(
            "Expected %s, got %s" % (pd.DataFrame.__name__, type(targets).__name__)
        )
    if edits.empty:  # empty edits dataframe
        edits = pd.DataFrame(
            columns=REPORT_COLS + ["START-DISTANCE", "STOP-DISTANCE", "FLAG"]
        )
        return edits
    # assign edit type (insertion, deletion, or snv)
    edits["TYPE"] = edits.apply(lambda x: assign_etype(x[4], x[5]), axis=1)
    # join targets and edits datasets
    edits = edits.merge(targets, on="SITE")
    # keep only the columns of interest
    if exp_type == EXPERIMENTS[0]:  # circleseq
        edits = edits[CIRCLESEQ_COLS]
    else:  # guideseq
        edits = edits[GUIDESEQ_COLS]
    edits.columns = REPORT_COLS  # rename columns
    # compute the distance between the called edits and their target site
    edits["START-DISTANCE"] = edits.apply(
        lambda x: compute_distance(int(x[2]), int(x[13])), axis=1
    )
    edits["STOP-DISTANCE"] = edits.apply(
        lambda x: compute_distance(int(x[2]), int(x[14])), axis=1
    )
    # assess if the eidts occurred inside the expected target site
    edits["FLAG"] = edits.apply(
        lambda x: edit_location(int(x[2]), int(x[13]), int(x[14]), x[3], x[4], x[12]),
        axis=1,
    )
    return edits


def _recover_variants_mutect2(targets, names_col, edits_dir, offregion, upstream):
    """(PRIVATE)
    Recover variants from VCF files returned by Mutect

    :param targets: guide target sites dataset
    :type targets: pd.DataFrame
    :param names_col: names column index
    :type names_col: int
    :param edits_dir: edits data directory
    :type edits_dir: str
    :param offregion: off regions
    :type offregion: bool
    :param upstream: off regions upstream
    :type upstream: bool
    :return: variants
    :rtype: List
    """
    file_suffix = "vcf.filtered.vcf"
    if offregion:  # offergion
        if upstream:  # upstream offregions
            fnames = [
                "%s.%s:%d-%d.%s"
                % (
                    x[names_col],
                    x[0],
                    int(x[1]) - 100 - PADSIZE,
                    int(x[1]) - 100,
                    file_suffix,
                )
                for _, x in targets.iterrows()
            ]
        else:  # downstream offregions
            fnames = [
                "%s.%s:%d-%d.%s"
                % (
                    x[names_col],
                    x[0],
                    int(x[2]) + 100,
                    int(x[2]) + 100 + PADSIZE,
                    file_suffix,
                )
                for _, x in targets.iterrows()
            ]
    else:  # onregions
        assert not upstream
        fnames = [
            "%s.%s:%d-%d.%s"
            % (
                x[names_col],
                x[0],
                int(x[1]) - PADSIZE,
                int(x[2]) + PADSIZE,
                file_suffix,
            )
            for _, x in targets.iterrows()
        ]
    return [read_vcf(os.path.join(edits_dir, fname)) for fname in fnames]


def report_mutect2(exp_type, guide, cell_type, outdir, offregion):
    """Construct the edits report from Mutect2 results

    :param exp_type: experiment type
    :type exp_type: str
    :param guide: guide
    :type guide: str
    :param cell_type: cell type
    :type cell_type: str
    :param outdir: output directory
    :type outdir: str
    :param offregion: off region report
    :type offregion: bool
    :raises TypeError: raise on exp_type type mismatch
    :raises ValueError: raise on forbidden exp_type value
    :raises TypeError: raise on guide type mismatch
    :raises ValueError: raise on forbidden guide value
    :raises TypeError: raise on cell_type type mismatch
    :raises ValueError: raise on forbidden cell_type value
    """
    if not isinstance(exp_type, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(exp_type).__name__))
    if exp_type not in EXPERIMENTS:
        raise ValueError("Forbidden experiment type (%s)" % (exp_type))
    if not isinstance(guide, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(guide).__name__))
    if not guide in GUIDES:
        raise ValueError("Forbidden guide (%s)" % (guide))
    if not isinstance(cell_type, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(cell_type).__name__)
        )
    if cell_type not in CELLTYPES:
        raise ValueError("Forbidden cell type (%s)" % (cell_type))
    # read target sites
    targets = read_targets(exp_type, guide)
    # parse VCFs
    region = "offregion" if offregion else "onregion"
    edits_dir = os.path.join(
        EDITS, VCALLINGTOOLS[0], exp_type, cell_type, region, guide
    )  # edits data
    if exp_type == EXPERIMENTS[0]:  # circleseq
        if offregion:
            # parse edits upstream (offregion)
            edits_upstream = _recover_variants_mutect2(
                targets, -1, edits_dir, offregion, True
            )
            # parse edits downstream (offregion)
            edits_downstream = _recover_variants_mutect2(
                targets, -1, edits_dir, offregion, False
            )
        else:
            edits = _recover_variants_mutect2(targets, -1, edits_dir, offregion, False)
    else:  # guideseq
        if offregion:
            # parse edits upstream (offregion)
            edits_upstream = _recover_variants_mutect2(
                targets, 6, edits_dir, offregion, True
            )
            # parse edits downstream (offregion)
            edits_downstream = _recover_variants_mutect2(
                targets, 6, edits_dir, offregion, False
            )
        else:
            edits = _recover_variants_mutect2(targets, 6, edits_dir, offregion, False)
    # build the edits dataframe
    if offregion:
        edits = pd.concat(
            [
                edits_dataframe(edits_upstream, targets.SITE.tolist()),
                edits_dataframe(edits_downstream, targets.SITE.tolist()),
            ]
        )
    else:
        edits = edits_dataframe(edits, targets.SITE.tolist())
    edits = process_edits_dataframe(edits, targets, exp_type)
    # write the report
    outfile = os.path.join(
        outdir, f"{VCALLINGTOOLS[0]}_{exp_type}_{cell_type}_{guide}_{region}.tsv"
    )
    edits.to_csv(outfile, sep="\t", index=False)


def recover_variants(targets, edits_dir, file_suffix, offregion, upstream):
    """Recover variants from VCF files returned by Strelka, Pindel and Varscan

    :param targets: guide target sites dataset
    :type targets: pd.DataFrame
    :param edits_dir: edits data directory
    :type edits_dir: str
    :param file_suffix: VCF file suffix
    :type file_suffix: str
    :param offregion: off regions
    :type offregion: bool
    :param upstream: off regions upstream
    :type upstream: bool
    :return: variants
    :rtype: List
    """
    if offregion:  # offregions
        if upstream:  # upstream offregions
            fnames = [
                "%s_%d_%d%s"
                % (x[0], int(x[1]) - 100 - PADSIZE, int(x[1]) - 100, file_suffix)
                for _, x in targets.iterrows()
            ]
        else:  # dowstream offregios
            fnames = [
                "%s_%d_%d%s"
                % (x[0], int(x[2]) + 100, int(x[2]) + 100 + PADSIZE, file_suffix)
                for _, x in targets.iterrows()
            ]
    else:  # onregions
        assert not upstream
        fnames = [
            "%s_%d_%d%s" % (x[0], int(x[1]) - PADSIZE, int(x[2]) + PADSIZE, file_suffix)
            for _, x in targets.iterrows()
        ]
    return [read_vcf(os.path.join(edits_dir, fname)) for fname in fnames]


def report_strelka(exp_type, guide, cell_type, outdir, offregion):
    """Construct the edits report from Strelka results

    :param exp_type: experiment type
    :type exp_type: str
    :param guide: guide
    :type guide: str
    :param cell_type: cell type
    :type cell_type: str
    :param outdir: output directory
    :type outdir: str
    :param offregion: off region report
    :type offregion: bool
    :raises TypeError: raise on exp_type type mismatch
    :raises ValueError: raise on forbidden exp_type value
    :raises TypeError: raise on guide type mismatch
    :raises ValueError: raise on forbidden guide value
    :raises TypeError: raise on cell_type type mismatch
    :raises ValueError: raise on forbidden cell_type value
    """
    if not isinstance(exp_type, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(exp_type).__name__))
    if exp_type not in EXPERIMENTS:
        raise ValueError("Forbidden experiment type (%s)" % (exp_type))
    if not isinstance(guide, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(guide).__name__))
    if not guide in GUIDES:
        raise ValueError("Forbidden guide (%s)" % (guide))
    if not isinstance(cell_type, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(cell_type).__name__)
        )
    if cell_type not in CELLTYPES:
        raise ValueError("Forbidden cell type (%s)" % (cell_type))
    # read the target sites
    targets = read_targets(exp_type, guide)
    # parse VCFs
    region = "offregion" if offregion else "onregion"
    edits_dir = os.path.join(
        EDITS, VCALLINGTOOLS[1], exp_type, cell_type, region, guide
    )
    if offregion:
        # recover upstream variants
        edits_upstream_snvs = recover_variants(
            targets, edits_dir, "_somatic.snvs.vcf", offregion, True
        )
        edits_upstream_indels = recover_variants(
            targets, edits_dir, "_somatic.indels.vcf", offregion, True
        )
        # recover downstream variants
        edits_downstream_snvs = recover_variants(
            targets, edits_dir, "_somatic.snvs.vcf", offregion, False
        )
        edits_downstream_indels = recover_variants(
            targets, edits_dir, "_somatic.indels.vcf", offregion, False
        )
    else:  # onregions
        edits_snvs = recover_variants(
            targets, edits_dir, "_somatic.snvs.vcf", offregion, False
        )
        edits_indels = recover_variants(
            targets, edits_dir, "_somatic.indels.vcf", offregion, False
        )
    # build the edits dataframe
    if offregion:
        edits = pd.concat(
            [
                pd.concat(
                    [
                        edits_dataframe(
                            edits_upstream_snvs, targets.SITE.tolist(), strelka=True
                        ),
                        edits_dataframe(
                            edits_upstream_indels, targets.SITE.tolist(), strelka=True
                        ),
                    ]
                ),
                pd.concat(
                    [
                        edits_dataframe(
                            edits_downstream_snvs, targets.SITE.tolist(), strelka=True
                        ),
                        edits_dataframe(
                            edits_downstream_indels, targets.SITE.tolist(), strelka=True
                        ),
                    ]
                ),
            ]
        )
    else:
        edits = pd.concat(
            [
                edits_dataframe(edits_snvs, targets.SITE.tolist(), strelka=True),
                edits_dataframe(edits_indels, targets.SITE.tolist(), strelka=True),
            ]
        )
    # build the edits dataset
    edits = process_edits_dataframe(edits, targets, exp_type)
    # write the report
    outfile = os.path.join(
        outdir, f"{VCALLINGTOOLS[1]}_{exp_type}_{cell_type}_{guide}_{region}.tsv"
    )
    edits.to_csv(outfile, sep="\t", index=False)


def report_pindel(exp_type, guide, cell_type, outdir, offregion):
    """Construct the edits report from Pindel results

    :param exp_type: experiment type
    :type exp_type: str
    :param guide: guide
    :type guide: str
    :param cell_type: cell type
    :type cell_type: str
    :param outdir: output directory
    :type outdir: str
    :param offregion: off region report
    :type offregion: bool
    :raises TypeError: raise on exp_type type mismatch
    :raises ValueError: raise on forbidden exp_type value
    :raises TypeError: raise on guide type mismatch
    :raises ValueError: raise on forbidden guide value
    :raises TypeError: raise on cell_type type mismatch
    :raises ValueError: raise on forbidden cell_type value
    """
    if not isinstance(exp_type, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(exp_type).__name__))
    if exp_type not in EXPERIMENTS:
        raise ValueError("Forbidden experiment type (%s)" % (exp_type))
    if not isinstance(guide, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(guide).__name__))
    if not guide in GUIDES:
        raise ValueError("Forbidden guide (%s)" % (guide))
    if not isinstance(cell_type, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(cell_type).__name__)
        )
    if cell_type not in CELLTYPES:
        raise ValueError("Forbidden cell type (%s)" % (cell_type))
    # read target sites
    targets = read_targets(exp_type, guide)
    # parse VCFs
    region = "offregion" if offregion else "onregion"
    edits_dir = os.path.join(
        EDITS, VCALLINGTOOLS[2], exp_type, cell_type, region, guide
    )
    if offregion:
        # recover upstream indels
        edits_upstream_deletions = recover_variants(
            targets, edits_dir, "_D.vcf", offregion, True
        )
        edits_upstream_insertions = recover_variants(
            targets, edits_dir, "_SI.vcf", offregion, True
        )
        # recover downstream indels
        edits_downstream_deletions = recover_variants(
            targets, edits_dir, "_D.vcf", offregion, False
        )
        edits_downstream_insertions = recover_variants(
            targets, edits_dir, "_SI.vcf", offregion, False
        )
    else:  # offregions
        edits_deletions = recover_variants(
            targets, edits_dir, "_D.vcf", offregion, False
        )
        edits_insertions = recover_variants(
            targets, edits_dir, "_SI.vcf", offregion, False
        )
    if offregion:
        edits = pd.concat(
            [
                pd.concat(
                    [
                        edits_dataframe(
                            edits_upstream_deletions, targets.SITE.tolist(), pindel=True
                        ),
                        edits_dataframe(
                            edits_upstream_insertions,
                            targets.SITE.tolist(),
                            pindel=True,
                        ),
                    ]
                ),
                pd.concat(
                    [
                        edits_dataframe(
                            edits_downstream_deletions,
                            targets.SITE.tolist(),
                            pindel=True,
                        ),
                        edits_dataframe(
                            edits_downstream_insertions,
                            targets.SITE.tolist(),
                            pindel=True,
                        ),
                    ]
                ),
            ]
        )
    else:
        edits = pd.concat(
            [
                edits_dataframe(edits_deletions, targets.SITE.tolist(), pindel=True),
                edits_dataframe(edits_insertions, targets.SITE.tolist(), pindel=True),
            ]
        )
    # build the edits dataframe
    edits = process_edits_dataframe(edits, targets, exp_type)
    outfile = os.path.join(
        outdir, f"{VCALLINGTOOLS[2]}_{exp_type}_{cell_type}_{guide}_{region}.tsv"
    )
    edits.to_csv(outfile, sep="\t", index=False)


def report_varscan(exp_type, guide, cell_type, outdir, offregion):
    """Construct the edits report from Varscan results

    :param exp_type: experiment type
    :type exp_type: str
    :param guide: guide
    :type guide: str
    :param cell_type: cell type
    :type cell_type: str
    :param outdir: output directory
    :type outdir: str
    :param offregion: off region report
    :type offregion: bool
    :raises TypeError: raise on exp_type type mismatch
    :raises ValueError: raise on forbidden exp_type value
    :raises TypeError: raise on guide type mismatch
    :raises ValueError: raise on forbidden guide value
    :raises TypeError: raise on cell_type type mismatch
    :raises ValueError: raise on forbidden cell_type value
    """
    if not isinstance(exp_type, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(exp_type).__name__))
    if exp_type not in EXPERIMENTS:
        raise ValueError("Forbidden experiment type (%s)" % (exp_type))
    if not isinstance(guide, str):
        raise TypeError("Expected %s, got %s" % (str.__name__, type(guide).__name__))
    if not guide in GUIDES:
        raise ValueError("Forbidden guide (%s)" % (guide))
    if not isinstance(cell_type, str):
        raise TypeError(
            "Expected %s, got %s" % (str.__name__, type(cell_type).__name__)
        )
    if cell_type not in CELLTYPES:
        raise ValueError("Forbidden cell type (%s)" % (cell_type))
    # read the target sites
    targets = read_targets(exp_type, guide)
    # parse VCFs
    region = "offregion" if offregion else "onregion"
    edits_dir = os.path.join(
        EDITS, VCALLINGTOOLS[3], exp_type, cell_type, region, guide
    )
    if offregion:  # offergions
        # upstream variants
        edits_upstream_snvs = recover_variants(
            targets, edits_dir, ".snp.vcf", offregion, True
        )
        edits_upstream_indels = recover_variants(
            targets, edits_dir, ".indel.vcf", offregion, True
        )
        # downstream variants
        edits_downstream_snvs = recover_variants(
            targets, edits_dir, ".snp.vcf", offregion, False
        )
        edits_downstream_indels = recover_variants(
            targets, edits_dir, ".indel.vcf", offregion, False
        )
    else:  # onregions
        edits_snvs = recover_variants(targets, edits_dir, ".snp.vcf", offregion, False)
        edits_indels = recover_variants(
            targets, edits_dir, ".indel.vcf", offregion, False
        )
    # build the edits dataframe
    if offregion:
        edits = pd.concat(
            [
                pd.concat(
                    [
                        edits_dataframe(
                            edits_upstream_snvs, targets.SITE.tolist(), varscan=True
                        ),
                        edits_dataframe(
                            edits_upstream_indels, targets.SITE.tolist(), varscan=True
                        ),
                    ]
                ),
                pd.concat(
                    [
                        edits_dataframe(
                            edits_downstream_snvs, targets.SITE.tolist(), varscan=True
                        ),
                        edits_dataframe(
                            edits_downstream_indels, targets.SITE.tolist(), varscan=True
                        ),
                    ]
                ),
            ]
        )
    else:
        edits = pd.concat(
            [
                edits_dataframe(edits_snvs, targets.SITE.tolist(), varscan=True),
                edits_dataframe(edits_indels, targets.SITE.tolist(), varscan=True),
            ]
        )
    # build the edits dataset
    edits = process_edits_dataframe(edits, targets, exp_type)
    # write the report
    outfile = os.path.join(
        outdir, f"{VCALLINGTOOLS[3]}_{exp_type}_{cell_type}_{guide}_{region}.tsv"
    )
    edits.to_csv(outfile, sep="\t", index=False)


def main():
    args = parse_commandline()
    if args.tool == VCALLINGTOOLS[0]:  # mutect2
        for guide in tqdm(GUIDES):
            for cell_type in CELLTYPES:
                report_mutect2(args.type, guide, cell_type, args.out, args.offregion)
    elif args.tool == VCALLINGTOOLS[1]:  # strelka
        for guide in tqdm(GUIDES):
            for cell_type in CELLTYPES:
                report_strelka(args.type, guide, cell_type, args.out, args.offregion)
    elif args.tool == VCALLINGTOOLS[2]:  # pindel
        for guide in tqdm(GUIDES):
            for cell_type in CELLTYPES:
                report_pindel(args.type, guide, cell_type, args.out, args.offregion)
    elif args.tool == VCALLINGTOOLS[3]:  # varscan
        for guide in tqdm(GUIDES):
            for cell_type in CELLTYPES:
                report_varscan(args.type, guide, cell_type, args.out, args.offregion)
    else:
        raise ValueError(
            f"Unable to find results for {args.tool}. Check the help for the available tools"
        )


if __name__ == "__main__":
    main()
