"""
Define global settings of  variables and functions

"""

import importlib
import os

import matplotlib.colors as mcolors

import nanocompare.legacy.performance_plots as pp
from global_config import pkl_base_dir

importlib.reload(pp)

# Project base dir for Nanocompare
nanocompare_basedir = "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare"

# Tool lists and abbreviations
tools = ['DeepSignal_calls', 'Tombo_calls', 'DeepMod_calls', 'Nanopolish_calls']
tools_abbr = ['DeepSignal', 'Tombo', 'DeepMod', 'Nanopolish']

# map_tools_to_abbr = {tools[0]: tools_abbr[0], tools[1]: tools_abbr[1], tools[2]: tools_abbr[2], tools[3]: tools_abbr[3]}

# tool -> abbr
dict_tools_to_abbr = {tools[i]: tools_abbr[i] for i in range(len(tools))}

# Locations of two types
locations_category = ["GW", "cpgIslandExt", "promoters_500bp", "exonFeature", "intergenic", "intronFeature"]
locations_singleton = ["singletons", "nonsingletons", "discordant", "concordant"]

locations_category2 = ["Genome-wide", "CpG Island", "Promoters", "Exons", "Intergenic", "Introns"]
locations_singleton2 = ["Singletons", "Non-singletons", "Discordant", "Concordant"]

# Plot performance order
perf_order = ['F1_5mC', 'F1_5C']

# Correlation ploting input tsv's important fields
cor_tsv_fields = ["DeepSignal_freq", "Tombo_freq", "Nanopolish_freq", "DeepMod_freq", "DeepMod_clust_freq", "BSseq"]

cor_tsv_fields_abbr = ["DeepSignal", "Tombo", "Nanopolish", "DeepMod", "DeepMod_clust", "BSseq"]

gbtruth_filedict = {
        'NA19240': os.path.join('/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/joined_reps/RRBS/extractBismark', 'NA19240_joined_RRBS.Read_R1.Rep_1_trimmed_bismark_bt2.bismark.cov.gz'),
        'K562'   : os.path.join('/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/joined_reps/RRBS/extractBismark', 'K562_joined_RRBS.Read_R1.Rep_4_trimmed_bismark_bt2.bismark.cov.gz'),
        'HL60'   : os.path.join('/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/joined_reps/RRBS/extractBismark', 'HL60_joined_RRBS.Read_R1.Rep_3_trimmed_bismark_bt2.bismark.cov.gz'),  # /projects/li-lab/wangj/AML-oxbs/AML-bs_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
        'APL'    : '/projects/li-lab/wangj/AML-oxbs/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz'
        }

singletonsFile = "hg38_singletons.bed"
nonsingletonsFile = "hg38_nonsingletons.bed"

narrowCoord = [False, singletonsFile, nonsingletonsFile, "ONT.hg38.cpgIslandExt.bed", "ONT.hg38.cpgShoresExt.bed", "ONT.hg38.cpgShelvesExt.bed", "ONT.hg38.exonFeature.bed", "ONT.hg38.geneFeature.bed", "ONT.hg38.intergenic.bed", "ONT.hg38.intronFeature.bed", "ONT.hg38.promoterFeature.flank_100.bed", "ONT.hg38.promoterFeature.flank_1000.bed",
        "ONT.hg38.promoterFeature.flank_200.bed", "ONT.hg38.promoterFeature.flank_2000.bed", "ONT.hg38.promoterFeature.flank_500.bed", "ONT.hg38.promoterFeature.flank_750.bed"]

narrowCoord = [False] + [os.path.join(pkl_base_dir, 'nanocompare', 'coordinate', cofn) for cofn in narrowCoord[1:]]

important_region_bed_fns = [narrowCoord[-2], narrowCoord[6], narrowCoord[9], narrowCoord[8], narrowCoord[3]]


def dict_cor_tsv_to_abbr():
    dict_cor = {cor_tsv_fields[i]: cor_tsv_fields_abbr[i] for i in range(len(cor_tsv_fields))}
    return dict_cor


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


# greenApple = make_colormap([ConvertRGB2sth(255, 255, 255), ConvertRGB2sth(173, 255, 50), 0.05, ConvertRGB2sth(173, 255, 50), ConvertRGB2sth(78, 121, 32), 0.3, ConvertRGB2sth(78, 121, 32), ConvertRGB2sth(0, 0, 0)])
# G2O = make_colormap([ConvertRGB2sth(0, 158, 115), ConvertRGB2sth(0, 0, 0), 0.5, ConvertRGB2sth(0, 0, 0), ConvertRGB2sth(230, 159, 0)])

# Color used of correlation plots
agressiveHot = make_colormap([ConvertRGB2sth(255, 255, 255), ConvertRGB2sth(255, 237, 8), 0.05, ConvertRGB2sth(255, 181, 0), ConvertRGB2sth(218, 33, 0), 0.3, ConvertRGB2sth(218, 33, 0), ConvertRGB2sth(0, 0, 0)])


def map_from_tool_to_abbr(tool):
    """
    From 'DeepSignal_calls' to 'DeepSignal'
    :param tool:
    :return:
    """
    return dict_tools_to_abbr[tool]


if __name__ == '__main__':

    pass
