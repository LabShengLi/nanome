#!/usr/bin/env python3

"""
Define names and global variables
"""
import json
import os

import pandas as pd

from nanocompare.global_config import data_base_dir, set_log_debug_level

nanome_version = "1.3.4"

# define the small error of 0 and 1, for fully-meth and unmeth eval
epslong = 1e-5

# can be 1.0, or 0.9 for the level of fully-meth definition
fully_meth_level = 1.0

# chr 1-22 X and Y
humanChrSet = [f'chr{k}' for k in range(1, 23)] + ['chrX', 'chrY']

ecoliChrSet = ['NC_000913.3']

datasets_order = ["NA12878", "NA19240", "APL", "K562", "HL60"]

ToolNameList = ['Nanopolish', 'Megalodon', 'DeepSignal', 'Guppy', 'Tombo', 'METEORE', 'DeepMod']

# read format for tools
ToolEncodeList = ['Nanopolish', 'Megalodon', 'Megalodon.ZW', 'DeepSignal',
                  'Guppy', 'Guppy.ZW', 'Guppy.gcf52ref', 'Tombo',
                  'METEORE', 'DeepMod', 'DeepMod.C', 'DeepMod.Cluster', ]

# read format for bs-seq
BGTruthEncodeList = ['bismark', 'encode']

ToolsColorList = ["#56B4E9", "#CC79A7", "#999999", "#009E73", "#E69F00", "#0072B2", "#D55E00"]

# default reference genome file location
referenceGenomeFile = "/projects/li-lab/Nanopore_compare/nf_input/reference_genome/hg38/hg38.fasta"

singletonsFile = "hg38_singletons_10bp.bed.gz"
nonsingletonsFile = "hg38_nonsingletons_10bp.bed.gz"
singletonFileExtStr = "_10bp.bed.gz"

narrowCoordNameList = ['x.x.Genome-wide',
                       singletonsFile, nonsingletonsFile,
                       "ONT.hg38.promoterFeature.flank_2000.bed.gz", "ONT.hg38.exonFeature.bed.gz",
                       "ONT.hg38.intronFeature.bed.gz", "ONT.hg38.intergenic.bed.gz",
                       "ONT.hg38.cpgIslandExt.bed.gz", "ONT.hg38.cpgShoresExt.bed.gz",
                       "ONT.hg38.cpgShelvesExt.bed.gz"]

# None means no coordinate used, i.e. Genome-wide
narrowCoordFileList = [None] + [os.path.join(data_base_dir, 'genome-annotation', cofn) for cofn in
                                narrowCoordNameList[1:]]

# CG density bed regions file basename
cg_density_coord_name_list = ["hg38.gc5Base_merged.cg_bin20.bed.gz", "hg38.gc5Base_merged.cg_bin40.bed.gz",
                              "hg38.gc5Base_merged.cg_bin60.bed.gz", "hg38.gc5Base_merged.cg_bin80.bed.gz",
                              "hg38.gc5Base_merged.cg_bin100.bed.gz"]
# CG density bed regions file full path
cg_density_file_list = [os.path.join(data_base_dir, 'genome-annotation', cofn) for cofn in cg_density_coord_name_list]

# Repetitive regions
rep_coord_name_list = ["hg38.repetitive.rep_SINE.bed.gz", "hg38.repetitive.rep_LINE.bed.gz",
                       "hg38.repetitive.rep_LTR.bed.gz", "hg38.repetitive.rep_DNA.bed.gz",
                       "hg38.repetitive.rep_Others.bed.gz"]
rep_file_list = [os.path.join(data_base_dir, 'genome-annotation', cofn) for cofn in rep_coord_name_list]

region_name_to_fn_dict = {
    'Genome-wide': narrowCoordFileList[0], 'Singletons': narrowCoordFileList[1],
    'Non-singletons': narrowCoordFileList[2],
    'Promoters': narrowCoordFileList[3], 'Exons': narrowCoordFileList[4], 'Introns': narrowCoordFileList[5],
    'Intergenic': narrowCoordFileList[6], 'CpG islands': narrowCoordFileList[7], 'CpG shores': narrowCoordFileList[8],
    'CpG shelves': narrowCoordFileList[9],
    'CG_20': cg_density_file_list[0], 'CG_40': cg_density_file_list[1],
    'CG_60': cg_density_file_list[2], 'CG_80': cg_density_file_list[3],
    'CG_100': cg_density_file_list[4],
    'rep_SINE': rep_file_list[0], 'rep_LINE': rep_file_list[1],
    'rep_LTR': rep_file_list[2], 'rep_DNA': rep_file_list[3],
    'rep_Others': rep_file_list[4]
}

# List of all start is 0-based bed files
list_base0_bed_basefn = ["ONT.hg38.cpgIslandExt.bed.gz", "ONT.hg38.cpgShoresExt.bed.gz",
                         "ONT.hg38.cpgShelvesExt.bed.gz"] + \
                        ["hg38.gc5Base_merged.cg_bin20.bed.gz", "hg38.gc5Base_merged.cg_bin40.bed.gz",
                         "hg38.gc5Base_merged.cg_bin60.bed.gz", "hg38.gc5Base_merged.cg_bin80.bed.gz",
                         "hg38.gc5Base_merged.cg_bin100.bed.gz"] + \
                        ["hg38.repetitive.rep_SINE.bed.gz", "hg38.repetitive.rep_LINE.bed.gz",
                         "hg38.repetitive.rep_LTR.bed.gz", "hg38.repetitive.rep_DNA.bed.gz",
                         "hg38.repetitive.rep_Others.bed.gz"]
enable_base_detection_bedfile = True
# enable_base_detection_bedfile = False


# Map each bed file name to a standard name, used by curve data plotting
location_filename_to_abbvname = {
    'x.x.Genome-wide': 'Genome-wide',
    singletonsFile: 'Singletons',
    nonsingletonsFile: 'Non-singletons',
    "ONT.hg38.cpgIslandExt.bed.gz": "CpG islands",
    "ONT.hg38.cpgShoresExt.bed.gz": 'CpG shores',
    "ONT.hg38.cpgShelvesExt.bed.gz": 'CpG shelves',
    "ONT.hg38.exonFeature.bed.gz": 'Exons',
    "ONT.hg38.intergenic.bed.gz": 'Intergenic',
    "ONT.hg38.intronFeature.bed.gz": 'Introns',
    "ONT.hg38.promoterFeature.flank_2000.bed.gz": 'Promoters',

    'hg38_nonsingletons.concordant.bed.gz': 'Concordant',
    'hg38_nonsingletons.discordant.bed.gz': 'Discordant',

    'hg38.gc5Base_merged.cg_bin20.bed.gz': 'CG_20',
    'hg38.gc5Base_merged.cg_bin40.bed.gz': 'CG_40',
    'hg38.gc5Base_merged.cg_bin60.bed.gz': 'CG_60',
    'hg38.gc5Base_merged.cg_bin80.bed.gz': 'CG_80',
    'hg38.gc5Base_merged.cg_bin100.bed.gz': 'CG_100',

    'hg38.repetitive.rep_SINE.bed.gz': 'rep_SINE',
    'hg38.repetitive.rep_LINE.bed.gz': 'rep_LINE',
    'hg38.repetitive.rep_LTR.bed.gz': 'rep_LTR',
    'hg38.repetitive.rep_DNA.bed.gz': 'rep_DNA',
    'hg38.repetitive.rep_Others.bed.gz': 'rep_Others',
}

# Tagname list of coordinate file list
narrowCoordFileTag = ['Genome-wide'] + [location_filename_to_abbvname[cofn] for cofn in narrowCoordNameList[1:]]
cgCoordFileTag = [location_filename_to_abbvname[cofn] for cofn in cg_density_coord_name_list]
repCoordFileTag = [location_filename_to_abbvname[cofn] for cofn in rep_coord_name_list]

# which column of performance table is extracted and returned
perf_report_columns = ['Dataset', 'Tool', 'Location', 'Accuracy', "Macro-F1", 'ROC-AUC', 'Average-Precision',
                       "Macro-Precision", "Macro-Recall", "Micro-F1", "Micro-Precision", "Micro-Recall", 'F1_5mC',
                       'F1_5C', 'Precision_5mC', 'Precision_5C', 'Recall_5mC', 'Recall_5C', 'mCsites_called',
                       'Csites_called', 'mCsites', 'Csites', 'referenceCpGs', 'prefix',
                       'coord']

# Rename raw name of region file (third section of .) to print name, such as 'ONT.hg38.cpgIslandExt.bed'
location_name_map_raw_to_standard = {
    'cpgIslandExt': 'CpG islands',
    'discordant': 'Discordant',
    'concordant': 'Concordant',
    'cpgShoresExt': 'CpG shores',
    'cpgShelvesExt': 'CpG shelves',
    'exonFeature': 'Exons',
    'intergenic': 'Intergenic',
    'intronFeature': 'Introns',
    'promoterFeature2000': 'Promoters',  # We use 2k bp promoter bed region
    'geneFeature': 'GeneFeature',
    'cg_bin20': 'CG_20',  # CG density 20%
    'cg_bin40': 'CG_40',
    'cg_bin60': 'CG_60',
    'cg_bin80': 'CG_80',
    'cg_bin100': 'CG_100'
}

locations_category = ["Promoters", "Exons", "Introns",
                      "Intergenic", "CpG islands",
                      "CpG shores", "CpG shelves"]

locations_singleton = ["Genome-wide", "Singletons", "Non-singletons",
                       "Concordant", "Discordant"]

# New introduced regions for CG density and repetitive
locations_new = ["CG_20", "CG_40", "CG_60", "CG_80", "CG_100"] + \
                ["rep_SINE", "rep_LINE", "rep_LTR", "rep_DNA", "rep_Others"]

locations_order = ["Genome-wide", "Singletons", "Non-singletons", "Concordant", "Discordant"] + \
                  ["Promoters", "Exons", "Introns", "Intergenic", "CpG islands", "CpG shores", "CpG shelves"] + \
                  ["CG_20", "CG_40", "CG_60", "CG_80", "CG_100"] + \
                  ["rep_SINE", "rep_LINE", "rep_LTR", "rep_DNA", "rep_Others"]


def get_tool_name(encode_name):
    if encode_name.find('.') != -1:  # cut of DeepMod.C or DeepMod.Cluster
        return encode_name[:encode_name.find('.')]
    return encode_name


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


default_config_name = 'nanome_genome_annotation.csv'


def load_genome_annotation_config(verbose=False):
    ## config file can be located at pwd, then the module dir
    search_list = [os.getcwd(), os.path.dirname(__file__)]

    for bdir in search_list:
        config_filepath = os.path.join(bdir, default_config_name)
        if os.path.exists(config_filepath):
            break
    if not os.path.exists(config_filepath):
        raise Exception(f"Can not find genome annotaion config file:{default_config_name} from {search_list}")
    if verbose:
        print(f"Load config from {config_filepath}", flush=True)

    df = pd.read_csv(config_filepath)
    ret1 = dict()  # tagname-> (fn, 0-1 format, )
    ret2 = dict()  # fn-> (tagname, 0-1 format, strand-sensi, )
    for index, row in df.iterrows():
        ret1[str(row['tagname']).strip()] = (str(row['filename']).strip(), int(row['format-0/1']),)
        ret2[str(row['filename']).strip()] = (str(row['tagname']).strip(), int(row['format-0/1']),
                                              str(row['strand-sensitive']).strip().upper() == 'Y',)
    if verbose:
        print(f"Config file loaded, results are below:", flush=True)
        print(json.dumps(ret2, indent=4), flush=True)
    return ret1, ret2


genome_wide_tagname = 'Genome-wide'
sing_tagname = 'Singletons'
nonsing_tagname = 'Non-singletons'
concord_tagname = 'Concordant'
discord_tagname = 'Discordant'
# key -> value, key=tagname, value = (filename, 0/1 format)
# from singletons, non-singletons, to genetic/intergenic, cg-density and repetitive regions
region_tagname_dict, region_filename_dict = load_genome_annotation_config()

if __name__ == '__main__':
    set_log_debug_level()
    print(region_tagname_dict)
    print(region_filename_dict)
