"""The script constructs the reports storing the edits called by the different 
tools considered in the study. The reports are stored as TSV files.
"""

from typing import List, Optional
from tqdm import tqdm

import pandas as pd

import argparse
import sys
import os


# data directories
BASEDIR = "/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS"
GUIDESEQ = os.path.join(
    BASEDIR, "offtargetDetection/casoffinder/offby6/CRISPRessoWGS/guideseq_anno"
)
CIRCLESEQ = os.path.join(BASEDIR, "offtargetDetection/circleseq/")
EDITS = os.path.join(
    BASEDIR, 
    "wgs/GM12878-Cas9/WGS1000/detectWithOtherTools/manuel_experiments/VCFs"
)
# TODO: remove --out 
# REPORTS = os.path.join(
#     BASEDIR, 
#     "wgs/GM12878-Cas9/WGS1000/detectWithOtherTools/manuel_experiments/reports"
# )
# variant calling tools
VCALLINGTOOLS = ["mutect2", "strelka", "pindel", "varscan"]
# guides
GUIDES = ["EMX1", "HEKSite4", "RNF2", "VEGFASite3"]
# experiments type
EXPERIMENTS = ["circleseq", "guideseq"]
# cell types
CELLTYPES = ["GM12878", "K562"]
# genomic regions padding size
PADSIZE = 10000
# variants type
INSERTION = "insertion"
DELETION = "deletion"
SNV = "snv"


def parse_commandline() -> argparse.ArgumentParser:
    """The function parses the command line input arguments.

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
        description="Script to construct the reports storing the called edits",
        usage="\n\tpython3 %(prog)s --tool <TOOL-NAME> --type <TYPE> --out "
              "<OUTDIR> --offregion"
    )
    parser.add_argument(
        "--tool", 
        type=str, 
        metavar="TOOL-NAME", 
        help="Variant calling tool. Available values: <mutect2, strelka, "
             "varscan, pindel>"
    )
    parser.add_argument(
        "--type", 
        type=str, 
        metavar="TYPE",
        help="Target sites validation experiment type. Available values: "
             "<guideseq, circleseq>"
    )
    parser.add_argument("--out", type=str, metavar="OUTDIR", help="Output directory")
    parser.add_argument(
        "--offregion",
        action="store_true",
        default=False,
        help="Shift the target regions 100bp upstream and downstream (not "
             "overlapping the original target site)"
    )
    args = parser.parse_args()
    return args


def read_targets(exp_type: str, guide: str) -> pd.DataFrame:
    """ (PRIVATE) 

    The function parses the target sites file.

    ...

    Parameters
    ----------
    exp_type
        Experiment type
    
    """

    if not isinstance(exp_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(exp_type).__name__}")
    if exp_type not in EXPERIMENTS:
        raise ValueError(f"Forbidden experiment type ({exp_type})")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not guide in GUIDES:
        raise ValueError(f"Forbidden guide ({guide})")
    if exp_type == EXPERIMENTS[0]:  # circleseq
        targets = pd.read_csv(
            os.path.join(CIRCLESEQ, f"{guide}.circleseq.hg19.hg38.targetname"),
            sep="\t"
        )
        cols = targets.columns.tolist()  # rename columns for later join
        targets.columns = cols[:-1] + ["SITE"]
    else:  # guideseq
        targets = pd.read_csv(os.path.join(GUIDESEQ, f"{guide}.guideseq"), sep="\t")  
        cols = targets.columns.tolist()  # rename columns for later join
        targets.columns = cols[:6] + ["SITE"] + cols[7:]
    assert not targets.empty
    return targets


def read_vcf(vcf: str) -> List[List]:
    """ (PRIVATE)

    The function parses the input VCF.

    ...

    Parameters
    ----------
    vcf
        Input VCF

    Returns
    -------
    List[List]
    """
    
    if not isinstance(vcf, str):
        raise TypeError(f"Expected {str.__name__}, got {type(vcf).__name__}")
    if not os.path.isfile(vcf):
        raise FileNotFoundError(f"Unable to locate {vcf}")
    try:
        handle = open(vcf, mode="r")
        variants = [
            line.strip().split() for line in handle if not line.startswith("#")
        ]
    except OSError:
        raise OSError(f"An error occurred while reading {vcf}")
    finally:
        handle.close()
    return variants


def edits_dataframe(edits: List, target_sites: List) -> pd.DataFrame:
    """ (PRIVATE) 

    The function constructs a dataframe from the parsed edits.

    ...

    Parameters
    ----------
    edits
        List of parsed edits
    target_sites
        List of target sites

    Returns
    -------
    pd.DataFrame
    """

    assert len(edits) == len(target_sites)
    # edits dictionary
    edf = {
        "SITE": [], "CHROM": [], "POS": [], "ID": [], "REF": [], "ALT": [], "FILTER": []
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
    edf = pd.DataFrame(edf)  # build dataframe from dict
    if edf.empty and any([bool(tse) for tse in edits]):
        raise ValueError(f"DataFrame is empty, but some edits have been called")
    return edf


def __etype(ref: str, alt: str) -> str:
    """ (PRIVATE)
    
    The function decide wheter the input edit is a deletion, insertion or snv.
    
    ...

    Parameters
    ----------
    ref
        Reference allele
    alt
        Alternative allele

    Returns
    -------
    str
    """

    if len(ref) < len(alt):
        return INSERTION
    elif len(alt) < len(ref):
        return DELETION
    return SNV


def assign_etype(ref: str, alt: str) -> str:
    """The function assigns a type to the input edit. The available types are 
    insertion, deletion, or snv.
    
    ...

    Parameters
    ----------
    ref
        Reference allele
    alt
        Alternative allele
    
    Returns
    -------
    str
    """

    if not isinstance(ref, str):
        raise TypeError(f"Expected {str.__name__}, got {type(ref).__name__}")
    if not isinstance(alt, str):
        raise TypeError(f"Expected {str.__name__}, got {type(alt).__name__}")
    if "," in alt:  # polyploid alternative allele
        etypes = []
        for aa in alt.split(","):
            etypes.append(__etype(ref, aa))
        return "-".join(list(set(etypes)))
    else:  # regular alternative allele
        return __etype(ref, alt)


def compute_distance(epos: int, target_pos: int) -> int:
    """The function computes the distance between the called edit edit site and 
    their supposed target site.
    
    ...

    Parameters
    ----------
    epos
        Edit position
    target_pos
        Target site coordinate

    Returns
    -------
    int
    """

    if not isinstance(epos, int):
        raise TypeError(f"Expected {int.__name__}, got {type(epos).__name__}")
    if epos < 0:
        raise ValueError(f"Forbidden edit position ({epos})")
    if not isinstance(target_pos, int):
        raise TypeError(f"Expected {int.__name__}, got {type(target_pos).__name__}")
    if target_pos < 0:
        raise ValueError(f"Forbidden edit position ({target_pos})")
    return epos - target_pos


def edit_location(
    epos: int, target_start: int, target_stop: int, ref: str, alt: str, etype: str
) -> str:
    """The function assesses if the input edit occurred inside or outside the 
    expected target site. If the variant is a deletion or insertion, the edit is 
    flagged as TP if it overlaps the target site positions.
    
    ...

    Parameters
    ----------
    epos
        Edit position
    target_start
        Target site start
    target_stop
        Target site stop
    ref
        Reference allele
    alt
        Alternative allele
    etype
        Edit type

    Returns
    -------
    str
    """

    if epos >= target_start and epos <= target_stop:
        return "TP"
    else:
        if "," in alt:  # polyploid alternative allele
            if INSERTION in etype or DELETION in etype:
                for aa in alt.split(","):
                    if len(aa) < len(ref):  # deletion:
                        padpos = list(range(epos, epos + len(ref)))
                        if any([p <= target_stop and p >= target_start for p in padpos]):
                            return "TP"
                    elif len(ref) < len(aa):  # insertion
                        padpos = list(range(epos, epos + len(aa)))
                        if any([p <= target_stop and p >= target_start for p in padpos]):
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


def report_mutect2(
    exp_type: str, guide: str, cell_type: str, outdir: str, offregion: bool
) -> None:
    """ The function construct the edits report using the variants called by 
    Mutect2.

    ...

    Parameters
    ----------
    exp_type
        Experiment type
    guide 
        Input guide
    cell_type
        Cell type
    outdir
        Output directory
    offregion   

    Returns
    -------
    pd.DataFrame
    """

    if not isinstance(exp_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(exp_type).__name__}")
    if exp_type not in EXPERIMENTS:
        raise ValueError(f"Forbidden experiment type ({exp_type})")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not guide in GUIDES:
        raise ValueError(f"Forbidden guide ({guide})")
    if not isinstance(cell_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(cell_type).__name__}")
    if cell_type not in CELLTYPES:
        raise ValueError(f"Forbidden cell type ({cell_type})")
    # read target sites
    targets = read_targets(exp_type, guide)
    # parse VCFs 
    region = "offregion" if offregion else "onregion"
    edits_dir = os.path.join(EDITS, VCALLINGTOOLS[0], exp_type, cell_type, region, guide)
    if exp_type == EXPERIMENTS[0]:  # circleseq
        if offregion:
            # parse edits upstream (offregion)
            edits_upstream = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir,
                        f"{x[-1]}.{x[0]}:{int(x[1]) - 100 - PADSIZE}-{int(x[1]) - 100}.vcf.filtered.vcf"
                    )
                ),
                axis=1
            )
            # parse edits downstream (offregion)
            edits_downstream = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir,
                        f"{x[-1]}.{x[0]}:{int(x[2]) + 100}-{int(x[2]) + 100 + PADSIZE}.vcf.filtered.vcf"
                    )
                ),
                axis=1
            )
        else:
            edits = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir, 
                        f"{x[-1]}.{x[0]}:{int(x[1]) - PADSIZE}-{int(x[2]) + PADSIZE}.vcf.filtered.vcf"
                    )
                ),
                axis=1
            )
    else:  # guideseq
        if offregion:
            # parse edits upstream (offregion)
            edits_upstream = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir,
                        f"{x[6]}.{x[0]}:{int(x[1]) - 100 - PADSIZE}-{int(x[1]) - 100}.vcf.filtered.vcf"
                    )
                ),
                axis=1
            )
            # parse edits downstream (offregion)
            edits_downstream = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir,
                        f"{x[6]}.{x[0]}:{int(x[2]) + 100}-{int(x[2]) + 100 + PADSIZE}.vcf.filtered.vcf"
                    )
                ),
                axis=1
            )
        else:
            edits = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir, 
                        f"{x[6]}.{x[0]}:{int(x[1]) - PADSIZE}-{int(x[2]) + PADSIZE}.vcf.filtered.vcf"
                    )
                ),
                axis=1
            )
    # build the edits dataframe 
    if offregion:
        edits = pd.concat(
                [
                    edits_dataframe(edits_upstream, targets.SITE.tolist()),
                    edits_dataframe(edits_downstream, targets.SITE.tolist())
                ]
            )
    else:
        edits = edits_dataframe(edits, targets.SITE.tolist())
    # assign edit type (insertion, deletion, or snv)
    edits["TYPE"] = edits.apply(lambda x : assign_etype(x[4], x[5]), axis=1)
    # join targets and edits datasets
    edits = edits.merge(targets, on="SITE")
    # keep only the columns of interest
    if exp_type == EXPERIMENTS[0]:  # circleseq
        keep = [
            "SITE", "CHROM", "POS", "REF", "ALT", "FILTER", "TYPE", "Start", "End", "Strand", "Distance"
        ]
    else:  # guideseq
        keep = [
            "SITE", "CHROM", "POS", "REF", "ALT", "FILTER", "TYPE", "start", "end", "Strand", "Mismatch Total"
        ]
    edits = edits[keep]
    edits.columns = [
        "SITE", "CHROM", "VAR-POS", "REF", "ALT", "FILTER", "TYPE", "TARGET-START", "TARGET-STOP", "STRAND", "MISMATCHES"
    ]  # rename columns
    # compute the distance between the called edits and their target site
    edits["START-DISTANCE"] = edits.apply(
        lambda x : compute_distance(int(x[2]), int(x[7])), axis=1
    )
    edits["STOP-DISTANCE"] = edits.apply(
        lambda x : compute_distance(int(x[2]), int(x[8])), axis=1
    )
    # assess if the edits occurred inside the expected target sites
    edits["FLAG"] = edits.apply(
        lambda x : edit_location(int(x[2]), int(x[7]), int(x[8]), x[3], x[4], x[6]),
        axis=1
    )
    # write the report
    outfile = os.path.join(
        outdir, f"{VCALLINGTOOLS[0]}_{exp_type}_{cell_type}_{guide}_{region}.tsv"
    )
    edits.to_csv(outfile, sep="\t", index=False)


def report_strelka(
    exp_type: str, guide: str, cell_type: str, outdir: str, offregion: bool
) -> None:
    """ The function construct the edits report using the variants called by 
    Mutect2.

    ...

    Parameters
    ----------
    exp_type
        Experiment type
    guide 
        Input guide
    cell_type
        Cell type
    outdir
        Output directory
    offregion   

    Returns
    -------
    pd.DataFrame
    """
    
    if not isinstance(exp_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(exp_type).__name__}")
    if exp_type not in EXPERIMENTS:
        raise ValueError(f"Forbidden experiment type ({exp_type})")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not guide in GUIDES:
        raise ValueError(f"Forbidden guide ({guide})")
    if not isinstance(cell_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(cell_type).__name__}")
    if cell_type not in CELLTYPES:
        raise ValueError(f"Forbidden cell type ({cell_type})")
    # read the target sites
    targets = read_targets(exp_type, guide)
    # parse VCFs
    region = "offregion" if offregion else "onregion"
    edits_dir = os.path.join(EDITS, VCALLINGTOOLS[1], exp_type, cell_type, region, guide)
    if offregion:
        edits_upstream_snvs = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir,
                        f"{x[0]}_{int(x[1]) - 100 - PADSIZE}_{int(x[1]) - 100}_somatic.snvs.vcf"
                    )
                ),
                axis=1
            )
        edits_upstream_indels = targets.apply(
            lambda x : read_vcf(
                os.path.join(
                    edits_dir,
                    f"{x[0]}_{int(x[1]) - 100 - PADSIZE}_{int(x[1]) - 100}_somatic.indels.vcf"
                )
            ),
            axis=1
        )
        edits_downstream_snvs = targets.apply(
            lambda x : read_vcf(
                os.path.join(
                    edits_dir,
                    f"{x[0]}_{int(x[2]) + 100}_{int(x[2]) + 100 + PADSIZE}_somatic.snvs.vcf"
                )
            ),
            axis=1
        )
        edits_downstream_indels = targets.apply(
            lambda x : read_vcf(
                os.path.join(
                    edits_dir,
                    f"{x[0]}_{int(x[2]) + 100}_{int(x[2]) + 100 + PADSIZE}_somatic.indels.vcf"
                )
            ),
            axis=1
        )
    else:
        edits_snvs = targets.apply(
                lambda x : read_vcf(
                    os.path.join(
                        edits_dir,
                        f"{x[0]}_{int(x[1]) - PADSIZE}_{int(x[2]) + PADSIZE}_somatic.snvs.vcf"
                    )
                ),
                axis=1 
            )
        edits_indels = targets.apply(
            lambda x : read_vcf(
                os.path.join(
                    edits_dir,
                    f"{x[0]}_{int(x[1]) - PADSIZE}_{int(x[2]) + PADSIZE}_somatic.indels.vcf"
                )
            ),
            axis=1
        )
    # build the edits dataframe
    if offregion:
        edits = pd.concat(
            pd.concat(
                [
                    edits_dataframe(edits_upstream_snvs, targets.SITE.tolist()),
                    edits_dataframe(edits_upstream_indels, targets.SITE.tolist())
                ]
            ),
            pd.concat(
                [
                    edits_dataframe(edits_downstream_snvs, targets.SITE.tolist()),
                    edits_dataframe(edits_downstream_indels, targets.SITE.tolist())
                ]
            )
        )
    else:
        edits = pd.concat(
            [
                edits_dataframe(edits_snvs, targets.SITE.tolist()),
                edits_dataframe(edits_indels, targets.SITE.tolist())
            ]
        )
    if not edits.empty:  # following operations only if edits have been called
        # assign edit type (insertion, deletion, or snv)
        edits["TYPE"] = edits.apply(lambda x : assign_etype(x[4], x[5]), axis=1)
        # join targets and edits dataset
        edits = edits.merge(targets, on="SITE")
        # keep only columns of interest
        if exp_type == EXPERIMENTS[0]:  # circleseq
            keep = [
                "SITE", "CHROM", "POS", "REF", "ALT", "FILTER", "TYPE", "Start", "End", "Strand", "Distance"
            ]
        else:  # guideseq
            keep = [
                "SITE", "CHROM", "POS", "REF", "ALT", "FILTER", "TYPE", "start", "end", "Strand", "Mismatch Total"
            ]
        edits = edits[keep]
        edits.columns = [
            "SITE", "CHROM", "VAR-POS", "REF", "ALT", "FILTER", "TYPE", "TARGET-START", "TARGET-STOP", "STRAND", "MISMATCHES"
        ]  # rename columns
        # compute the distance between the called edits and their target site
        edits["START-DISTANCE"] = edits.apply(
            lambda x : compute_distance(int(x[2]), int(x[7])), axis=1
        )
        edits["STOP-DISTANCE"] = edits.apply(
            lambda x : compute_distance(int(x[2]), int(x[8])), axis=1
        )
        # assess if the edits occurred inside the expected target sites
        edits["FLAG"] = edits.apply(
            lambda x : edit_location(int(x[2]), int(x[7]), int(x[8]), x[3], x[4], x[6]),
            axis=1
        )
    else:
        columns = [
            "SITE", "CHROM", "VAR-POS", "REF", "ALT", "FILTER", "TYPE", "TARGET-START", "TARGET-STOP", "STRAND", "MISMATCHES"
        ]
        edits = pd.DataFrame(columns=columns)
    # write the report
    outfile = os.path.join(
        outdir, f"{VCALLINGTOOLS[1]}_{exp_type}_{cell_type}_{guide}_{region}.tsv"
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
        pass
    elif args.tool == VCALLINGTOOLS[3]:  # varscan
        pass
    else:
        raise ValueError(
            f"Unable to find results for {args.tool}. Check the help for the available tools"
        )


if __name__ == "__main__":
    main()
        

