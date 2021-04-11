"""
Define global settings of  variables and functions

"""

import importlib
import os

import matplotlib.colors as mcolors

import nanocompare.legacy.performance_plots as pp
from nanocompare.global_config import data_base_dir, results_base_dir

importlib.reload(pp)

# chr 1-22 X and Y
humanChrSet = [f'chr{k}' for k in range(1, 23)] + ['chrX', 'chrY']

ToolNameList = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon']
Top3ToolNameList = ['DeepSignal', 'Nanopolish', 'Megalodon']

ToolEncodeList = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'DeepMod.C', 'DeepMod.Cluster', 'Megalodon']
BGTruthEncodeList = ['bismark', 'encode']  # 'bed',
# ToolsColorList = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal'])

ToolsColorList = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442"]

referenceGenomeFile = '/projects/li-lab/Ziwei/Nanopore/data/reference/hg38.fa'

# These two files are defined from Reference Genome
# Singletons:       XXXXXCGXXXXX            >=k bp CpGs
# Nonsingletons:    XXXXXCGXXXCGXXXXCGXX    <k bp for pair of neighbors of CpGs
# singletonsFile = "hg38_singletons_10bp.bed"
# nonsingletonsFile = "hg38_nonsingletons_10bp.bed"
# singletonFileExtStr = "_10bp.bed"

singletonsFile = "hg38_singletons.bed"
nonsingletonsFile = "hg38_nonsingletons.bed"
singletonFileExtStr = ".bed"

narrowCoordNameList = ['x.x.GenomeWide', singletonsFile, nonsingletonsFile, "ONT.hg38.cpgIslandExt.bed", "ONT.hg38.cpgShoresExt.bed", "ONT.hg38.cpgShelvesExt.bed", "ONT.hg38.exonFeature.bed", "ONT.hg38.geneFeature.bed", "ONT.hg38.intergenic.bed", "ONT.hg38.intronFeature.bed", "ONT.hg38.promoterFeature.flank_100.bed",
        "ONT.hg38.promoterFeature.flank_1000.bed",
        "ONT.hg38.promoterFeature.flank_200.bed", "ONT.hg38.promoterFeature.flank_2000.bed", "ONT.hg38.promoterFeature.flank_500.bed", "ONT.hg38.promoterFeature.flank_750.bed"]

# Map each bed file name to a standard name, used by curve data plotting
location_filename_to_abbvname = {
        'x.x.GenomeWide'                         : 'Genomewide',
        singletonsFile                           : 'Singleton',
        nonsingletonsFile                        : 'Nonsingleton',
        "ONT.hg38.cpgIslandExt.bed"              : "CpG Island",
        "ONT.hg38.cpgShoresExt.bed"              : 'CpG Shores',
        "ONT.hg38.cpgShelvesExt.bed"             : 'CpG Shelves',
        "ONT.hg38.exonFeature.bed"               : 'Exons',
        "ONT.hg38.geneFeature.bed"               : 'GeneFeature',
        "ONT.hg38.intergenic.bed"                : 'Intergenic',
        "ONT.hg38.intronFeature.bed"             : 'Introns',
        "ONT.hg38.promoterFeature.flank_100.bed" : 'Promoter_flank100',
        "ONT.hg38.promoterFeature.flank_1000.bed": 'Promoter_flank1000',
        "ONT.hg38.promoterFeature.flank_200.bed" : 'Promoter_flank200',
        "ONT.hg38.promoterFeature.flank_2000.bed": 'Promoters',  # 2k bp promoters used
        "ONT.hg38.promoterFeature.flank_500.bed" : 'Promoter_flank500',
        "ONT.hg38.promoterFeature.flank_750.bed" : 'Promoter_flank750',
        # 'hg38_singletons.absolute.bed'           : 'Absolute',
        'hg38_nonsingletons.concordant.bed'      : 'Concordant',
        'hg38_nonsingletons.discordant.bed'      : 'Discordant'
        }

# None means no coordinate used, i.e. Genome-wide
narrowCoordFileList = [None] + [os.path.join(data_base_dir, 'genome-annotation', cofn) for cofn in narrowCoordNameList[1:]]

narrowCoordFileTag = ['Genomewide'] + [location_filename_to_abbvname[cofn] for cofn in narrowCoordNameList[1:]]

# which column of performance table is extracted and returned
perf_report_columns = ['Dataset', 'Tool', 'Location', 'Accuracy', "Macro-F1", 'ROC-AUC', 'Average-Precision', "Macro-Precision", "Macro-Recall", "Micro-F1", "Micro-Precision", "Micro-Recall", 'F1_5mC', 'F1_5C', 'Precision_5mC', 'Precision_5C', 'Recall_5mC', 'Recall_5C', 'mCsites_called', 'Csites_called', 'mCsites', 'Csites', 'referenceCpGs', 'prefix',
        'coord']

# Rename raw name of retion to print name
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
        # 'absolute'           : 'Absolute',
        'geneFeature'        : 'GeneFeature'}

locations_category = ["Genome-wide", "CpG Island", "Promoters", "Exons", "Intergenic", "Introns", "CpG Shores", "CpG Shelves", "GeneFeature"]
locations_singleton = ["Singletons", "Non-singletons", "Discordant", "Concordant"] #TODO: Nonsingletons

runPrefixDefault = {
        'K562_WGBS'   : os.path.join(results_base_dir, 'MethPerf-K562_WGBS'),
        'HL60_RRBS'   : os.path.join(results_base_dir, 'MethPerf-HL60_RRBS'),
        'APL_RRBS'    : os.path.join(results_base_dir, 'MethPerf-APL_RRBS'),
        'NA19240_RRBS': os.path.join(results_base_dir, 'MethPerf-NA19240_RRBS'),
        }


def ConvertRGB2sth(r, g, b):
    tmp_r = (r / 255.0)
    tmp_g = (g / 255.0)
    tmp_b = (b / 255.0)

    return tmp_r, tmp_g, tmp_b


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def get_tool_name(toolname):
    if toolname.find('.') != -1:  # cut of DeepMod.C or DeepMod.Cluster
        return toolname[:toolname.find('.')]
    return toolname


def rename_location_from_coordinate_name(df):
    """
    Rename and change raw values of report df to more meaning full for display
    :param df:
    :return:
    """
    # Rename column names
    # df = df.rename(columns=raw_to_standard_perf_colname)

    # Replace coordinate name, and define unified Location column
    # df = df.replace(to_replace="False", value="x.x.Genome-wide")
    df = df.replace(to_replace=singletonsFile, value="x.x.Singletons")
    df = df.replace(to_replace=nonsingletonsFile, value="x.x.Nonsingletons")
    df['coord'] = df['coord'].str.replace("promoterFeature.flank_", "promoterFeature")

    # coord file like: x.x.Singletons, x.x.Non-singletons, HL60_RRBS_2Reps.hg38_nonsingletons.discordant.bed
    # we select the third of . split sections
    df["Location"] = df["coord"].str.split(".", n=3, expand=True)[2]
    df['Location'] = df['Location'].replace(location_name_map_raw_to_standard)
    return df


# Color used of correlation plots
agressiveHot = make_colormap([ConvertRGB2sth(255, 255, 255), ConvertRGB2sth(255, 237, 8), 0.05, ConvertRGB2sth(255, 181, 0), ConvertRGB2sth(218, 33, 0), 0.3, ConvertRGB2sth(218, 33, 0), ConvertRGB2sth(0, 0, 0)])

if __name__ == '__main__':

    pass
