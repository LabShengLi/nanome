#!/usr/bin/env python3

"""
Define names and global variables

"""

import os

# import nanocompare.legacy.performance_plots as pp
from nanocompare.global_config import data_base_dir

# importlib.reload(pp)

# define the small error of 0 and 1, for fully-meth and unmeth eval
epslong = 1e-5

# can be 1.0, or 0.9 for the level of fully-meth definition
fully_meth_level = 1.0

# chr 1-22 X and Y
humanChrSet = [f'chr{k}' for k in range(1, 23)] + ['chrX', 'chrY']
ecoliChrSet = ['NC_000913.3']

datasets_order = ["NA19240", "APL", "K562", "HL60"]

ToolNameList = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon']
# TODO: check all consistent, now just for plot figure
ToolNameList = ['Nanopolish', 'Megalodon', 'DeepSignal', 'Guppy', 'Tombo', 'METEORE', 'DeepMod']

# TODO: change to nanopolish, megalodon and deepsignal order later
Top3ToolNameList = ['Nanopolish', 'Megalodon', 'DeepSignal']

ToolEncodeList = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod.C', 'DeepMod.Cluster', 'Megalodon', 'Megalodon.ZW', 'Guppy', 'Guppy.ZW', 'Guppy.gcf52ref', 'METEORE']

BGTruthEncodeList = ['bismark', 'encode']  # 'bed',

ToolsColorList = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442"]
ToolsColorList = ["#56B4E9", "#CC79A7", "#999999", "#009E73", "#E69F00", "#0072B2", "#D55E00"]

referenceGenomeFile = '/projects/li-lab/Ziwei/Nanopore/data/reference/hg38.fa'
referenceGenomeFile = "/projects/li-lab/Nanopore_compare/nf_input/reference_genome/hg38/hg38.fasta"
# These two files are defined from Reference Genome
# Singletons:       XXXXXCGXXXXX            >=k bp distance CpGs
# Nonsingletons:    XXXXXCGXXXCGXXXXCGXX    <k bp for pair of neighbors of CpGs
singletonsFile = "hg38_singletons_10bp.bed"
nonsingletonsFile = "hg38_nonsingletons_10bp.bed"
singletonFileExtStr = "_10bp.bed"

# singletonsFile = "hg38_singletons.bed"
# nonsingletonsFile = "hg38_nonsingletons.bed"
# singletonFileExtStr = ".bed"

narrowCoordNameList = ['x.x.Genome-wide', singletonsFile, nonsingletonsFile, "ONT.hg38.cpgIslandExt.bed", "ONT.hg38.cpgShoresExt.bed", "ONT.hg38.cpgShelvesExt.bed", "ONT.hg38.exonFeature.bed", "ONT.hg38.geneFeature.bed", "ONT.hg38.intergenic.bed", "ONT.hg38.intronFeature.bed", "ONT.hg38.promoterFeature.flank_100.bed",
        "ONT.hg38.promoterFeature.flank_1000.bed",
        "ONT.hg38.promoterFeature.flank_200.bed", "ONT.hg38.promoterFeature.flank_2000.bed", "ONT.hg38.promoterFeature.flank_500.bed", "ONT.hg38.promoterFeature.flank_750.bed"]

# None means no coordinate used, i.e. Genome-wide
narrowCoordFileList = [None] + [os.path.join(data_base_dir, 'genome-annotation', cofn) for cofn in narrowCoordNameList[1:]]

# CG density bed regions file basename
cg_density_coord_name_list = ["hg38.gc5Base_merged.cg_bin20.bed.gz", "hg38.gc5Base_merged.cg_bin40.bed.gz", "hg38.gc5Base_merged.cg_bin60.bed.gz", "hg38.gc5Base_merged.cg_bin80.bed.gz", "hg38.gc5Base_merged.cg_bin100.bed.gz"]
# CG density bed regions file full path
cg_density_file_list = [os.path.join(data_base_dir, 'genome-annotation', cofn) for cofn in cg_density_coord_name_list]

# Repetitive regions
rep_coord_name_list = ["hg38.repetitive.rep_All.bed.gz", "hg38.repetitive.rep_SINE.bed.gz", "hg38.repetitive.rep_LINE.bed.gz", "hg38.repetitive.rep_LTR.bed.gz", "hg38.repetitive.rep_Simple_repeat.bed.gz",
        "hg38.repetitive.rep_DNA.bed.gz", "hg38.repetitive.rep_Low_complexity.bed.gz", "hg38.repetitive.rep_Satellite.bed.gz", "hg38.repetitive.rep_Retroposon.bed.gz", "hg38.repetitive.rep_Others.bed.gz"
        ]
rep_file_list = [os.path.join(data_base_dir, 'genome-annotation', cofn) for cofn in rep_coord_name_list]

# Map each bed file name to a standard name, used by curve data plotting
location_filename_to_abbvname = {
        'x.x.Genome-wide'                          : 'Genome-wide',
        singletonsFile                             : 'Singletons',
        nonsingletonsFile                          : 'Non-singletons',
        "ONT.hg38.cpgIslandExt.bed"                : "CpG Island",
        "ONT.hg38.cpgShoresExt.bed"                : 'CpG Shores',
        "ONT.hg38.cpgShelvesExt.bed"               : 'CpG Shelves',
        "ONT.hg38.exonFeature.bed"                 : 'Exons',
        "ONT.hg38.geneFeature.bed"                 : 'GeneFeature',
        "ONT.hg38.intergenic.bed"                  : 'Intergenic',
        "ONT.hg38.intronFeature.bed"               : 'Introns',
        "ONT.hg38.promoterFeature.flank_100.bed"   : 'Promoter_flank100',
        "ONT.hg38.promoterFeature.flank_1000.bed"  : 'Promoter_flank1000',
        "ONT.hg38.promoterFeature.flank_200.bed"   : 'Promoter_flank200',
        "ONT.hg38.promoterFeature.flank_2000.bed"  : 'Promoters',  # 2k bp promoters used
        "ONT.hg38.promoterFeature.flank_500.bed"   : 'Promoter_flank500',
        "ONT.hg38.promoterFeature.flank_750.bed"   : 'Promoter_flank750',
        'hg38_nonsingletons.concordant.bed'        : 'Concordant',
        'hg38_nonsingletons.discordant.bed'        : 'Discordant',

        'hg38.gc5Base_merged.cg_bin20.bed.gz'      : 'CG_20',
        'hg38.gc5Base_merged.cg_bin40.bed.gz'      : 'CG_40',
        'hg38.gc5Base_merged.cg_bin60.bed.gz'      : 'CG_60',
        'hg38.gc5Base_merged.cg_bin80.bed.gz'      : 'CG_80',
        'hg38.gc5Base_merged.cg_bin100.bed.gz'     : 'CG_100',

        'hg38.repetitive.rep_All.bed.gz'           : 'reg_All',
        'hg38.repetitive.rep_SINE.bed.gz'          : 'rep_SINE',
        'hg38.repetitive.rep_LINE.bed.gz'          : 'rep_LINE',
        'hg38.repetitive.rep_LTR.bed.gz'           : 'rep_LTR',
        'hg38.repetitive.rep_Simple_repeat.bed.gz' : 'rep_Simple_repeat',
        'hg38.repetitive.rep_DNA.bed.gz'           : 'rep_DNA',
        'hg38.repetitive.rep_Low_complexity.bed.gz': 'rep_Low_complexity',
        'hg38.repetitive.rep_Satellite.bed.gz'     : 'rep_Satellite',
        'hg38.repetitive.rep_Retroposon.bed.gz'    : 'rep_Retroposon',
        'hg38.repetitive.rep_Others.bed.gz'        : 'rep_Others',
        }

# Tagname list of coordinate file list
narrowCoordFileTag = ['Genome-wide'] + [location_filename_to_abbvname[cofn] for cofn in narrowCoordNameList[1:]]
cgCoordFileTag = [location_filename_to_abbvname[cofn] for cofn in cg_density_coord_name_list]
repCoordFileTag = [location_filename_to_abbvname[cofn] for cofn in rep_coord_name_list]

# which column of performance table is extracted and returned
perf_report_columns = ['Dataset', 'Tool', 'Location', 'Accuracy', "Macro-F1", 'ROC-AUC', 'Average-Precision', "Macro-Precision", "Macro-Recall", "Micro-F1", "Micro-Precision", "Micro-Recall", 'F1_5mC', 'F1_5C', 'Precision_5mC', 'Precision_5C', 'Recall_5mC', 'Recall_5C', 'mCsites_called', 'Csites_called', 'mCsites', 'Csites', 'referenceCpGs', 'prefix',
        'coord']

# Rename raw name of region file (third section of .) to print name, such as 'ONT.hg38.cpgIslandExt.bed'
location_name_map_raw_to_standard = {
        'cpgIslandExt'       : 'CpG Island',
        'discordant'         : 'Discordant',
        'concordant'         : 'Concordant',
        'cpgShoresExt'       : 'CpG Shores',
        'cpgShelvesExt'      : 'CpG Shelves',
        'exonFeature'        : 'Exons',
        'intergenic'         : 'Intergenic',
        'intronFeature'      : 'Introns',
        'promoterFeature2000': 'Promoters',  # We use 2k bp promoter bed region
        'geneFeature'        : 'GeneFeature',
        'cg_bin20'           : 'CG_20',  # CG density 20%
        'cg_bin40'           : 'CG_40',
        'cg_bin60'           : 'CG_60',
        'cg_bin80'           : 'CG_80',
        'cg_bin100'          : 'CG_100'
        }

locations_category = ["Genome-wide", "CpG Island", "Promoters", "Exons", "Intergenic", "Introns", "CpG Shores", "CpG Shelves", "GeneFeature"]
locations_singleton = ["Singletons", "Non-singletons", "Discordant", "Concordant"]  # TODO: Nonsingletons

# New introduced regions for CG density and repetitive
locations_new = ["CG_20", "CG_40", "CG_60", "CG_80", "CG_100"] + ["rep_All", "rep_SINE", "rep_LINE", "rep_LTR", "rep_Simple_repeat", "rep_DNA", "rep_Low_complexity", "rep_Satellite", "rep_Retroposon", "rep_Others"]

locations_order = ["Genome-wide", "Singletons", "Non-singletons", "Discordant", "Concordant", "CpG Island", "Promoters", "Exons", "Intergenic", "Introns", "CpG Shores", "CpG Shelves", "GeneFeature"] + ["CG_20", "CG_40", "CG_60", "CG_80", "CG_100"] + ["rep_All", "rep_SINE", "rep_LINE", "rep_LTR", "rep_Simple_repeat", "rep_DNA", "rep_Low_complexity",
        "rep_Satellite", "rep_Retroposon", "rep_Others"]


def get_tool_name(encode_name):
    if encode_name.find('.') != -1:  # cut of DeepMod.C or DeepMod.Cluster
        return encode_name[:encode_name.find('.')]
    return encode_name


def rename_location_from_coordinate_name(df):
    """
    Rename and change raw values of report df to more meaning full for display
    :param df:
    :return:
    """
    df = df.replace(to_replace=singletonsFile, value="x.x.Singletons")
    df = df.replace(to_replace=nonsingletonsFile, value="x.x.Non-singletons")
    df['coord'] = df['coord'].str.replace("promoterFeature.flank_", "promoterFeature")  # merge third section together for promoter file

    # coord file like: x.x.Singletons, x.x.Non-singletons, HL60_RRBS_2Reps.hg38_nonsingletons.discordant.bed
    # we select the third of . split sections
    df["Location"] = df["coord"].str.split(".", n=3, expand=True)[2]
    df['Location'] = df['Location'].replace(location_name_map_raw_to_standard)
    return df


def save_done_file(outdir, filename="DONE.txt"):
    """
    Save a done file for finish flag
    :param outdir:
    :param filename:
    :return:
    """
    outfn = os.path.join(outdir, filename)
    with open(outfn, "w") as outf:
        outf.write("DONE\n")


if __name__ == '__main__':

    pass
