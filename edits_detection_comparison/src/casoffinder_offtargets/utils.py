import os

BASEDIR = "/path/to/root/folder/"
VCALLINGTOOLS = ["mutect2", "strelka", "pindel", "varscan"]
GUIDES = ["EMX1", "HEKSite4", "RNF2", "VEGFASite3"]
CELLTYPES = ["GM12878", "K562"]
OFFTARGETS = os.path.join(BASEDIR, "offtargetDetection/casoffinder")
GENOME = "/path/to/genome/fasta/"
BAMS = os.path.join(BASEDIR, "/path/to/bam/")
PINDEL_BAMS = os.path.join(BASEDIR, "/path/to/bam/")
OUTDIR = "/path/to/out/folder/"

def _create_celltype_dirtree(root):
    """(PRIVATE)
    Build the directory tree for each cell type (GM12878 and K562)

    :param root: root directory
    :type root: str
    """
    assert isinstance(root, str)
    for guide in GUIDES:
        guide_dir = os.path.join(root, guide)
        if not os.path.isdir(guide_dir):
            os.mkdir(guide_dir)


def create_result_dirtree(tool):
    """Build the directory tree for the tool

    :param tool: variant calling tool
    :type tool: str
    """
    assert tool in VCALLINGTOOLS
    # if there is no tool tree directory, build the tree
    tool_root_dir = os.path.join(OUTDIR, tool)
    if not os.path.isdir(tool_root_dir):
        os.mkdir(tool_root_dir)
    # check GM12878 directories
    gm12878_dir = os.path.join(tool_root_dir, CELLTYPES[0])
    if not os.path.isdir(gm12878_dir):
        os.mkdir(gm12878_dir)
    _create_celltype_dirtree(gm12878_dir)  # create cell types dir tree
    # check K562 directories
    k562_dir = os.path.join(tool_root_dir, CELLTYPES[1])
    if not os.path.isdir(k562_dir):
        os.mkdir(k562_dir)
    _create_celltype_dirtree(k562_dir)  # create cell types dir tree
