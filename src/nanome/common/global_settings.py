#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : global_settings.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Define names and global variables
"""
import json
import os
from pathlib import Path

import pandas as pd
from nanome.common.global_config import set_log_debug_level, current_time_str

NANOME_VERSION = "2.0.6"

# define the small error of 0 and 1, for fully-meth and unmeth eval
EPSLONG = 1e-5

# dataframe loading chuncksize default
CHUNKSIZE = 500000

# can be 1.0, or 0.9 for the level of fully-meth definition
FULLY_METH_LEVEL = 1.00

# chr 1-22 X and Y
HUMAN_CHR_SET = [f'chr{k}' for k in range(1, 23)] + ['chrX', 'chrY']

ECOLI_CHR_SET = ['NC_000913.3']

datasets_order = ["NA12878", "NA19240", "APL", "K562", "HL60"]

ToolNameList = ['Nanopolish', 'Megalodon', 'DeepSignal', 'Guppy', 'Tombo', 'METEORE', 'DeepMod']

# read/site format encode name for tools
ToolEncodeList = ['Nanopolish', 'Megalodon', 'Megalodon.ZW', 'DeepSignal',
                  'Guppy', 'Guppy.ZW', 'Guppy.gcf52ref', 'Tombo',
                  'METEORE', 'DeepMod', 'DeepMod.C', 'DeepMod.Cluster',
                  'NANOME', 'UNIREAD', 'UNISITE']

# format for bs-seq
BGTruthEncodeList = ['bismark', 'encode', 'UNISITE']

ToolsColorList = ["#56B4E9", "#CC79A7", "#999999", "#009E73", "#E69F00", "#0072B2", "#D55E00"]

# default reference genome file location
reference_genome_hg38_fn = "/projects/li-lab/Nanopore_compare/nf_input/reference_genome/hg38/hg38.fasta"

enable_base_detection_bedfile = True
# enable_base_detection_bedfile = False


# Map each bed file name to a standard name, used by curve data plotting
location_filename_to_abbvname = {}

# which column of performance table is extracted and returned
perf_report_columns = ['Dataset', 'Tool', 'Location', 'Accuracy', "Macro-F1", 'ROC-AUC', 'Average-Precision',
                       "Macro-Precision", "Macro-Recall", "Micro-F1", "Micro-Precision", "Micro-Recall", 'F1_5mC',
                       'F1_5C', 'Precision_5mC', 'Precision_5C', 'Recall_5mC', 'Recall_5C', 'mCsites_called',
                       'Csites_called', 'mCsites', 'Csites', 'referenceCpGs', 'prefix',
                       'coord']

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

## Commonly used tagname for regions
genome_wide_tagname = 'Genome-wide'
sing_tagname = 'Singletons'
nonsing_tagname = 'Non-singletons'
concord_tagname = 'Concordant'
discord_tagname = 'Discordant'

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

    ret1 = dict()  # tagname-> (fn, 0-1 format, )
    ret2 = dict()  # fn-> (tagname, 0-1 format, strand-sensi, )

    try:
        df = pd.read_csv(config_filepath)
        for index, row in df.iterrows():
            ret1[str(row['tagname']).strip()] = (str(row['filename']).strip(), int(row['format-0/1']),)
            ret2[str(row['filename']).strip()] = (str(row['tagname']).strip(), int(row['format-0/1']),
                                                  str(row['strand-sensitive']).strip().upper() == 'Y',)
    except:
        print(f"ERROR: occur error when reading {config_filepath}", flush=True)
    if verbose:
        print(f"Config file loaded, results are below:")
        print(json.dumps(ret2, indent=4), flush=True)
    return ret1, ret2


# key -> value, key=tagname, value = (filename, 0/1 format)
# from singletons, non-singletons, to genetic/intergenic, cg-density and repetitive regions
region_tagname_dict, region_filename_dict = load_genome_annotation_config()


def get_tool_name(encode_name):
    if encode_name.find('.') != -1:  # cut of DeepMod.C or DeepMod.Cluster
        return encode_name[:encode_name.find('.')]
    return encode_name


def save_done_file(outdir, filename=None):
    """
    Save a done file for finish flag
    :param outdir:
    :param filename:
    :return:
    """
    time_tag = current_time_str()
    if filename == None:
        filename = f"DONE_{time_tag}.txt"
    outfn = os.path.join(outdir, filename)
    with open(outfn, "w") as outf:
        outf.write(f"DONE at {time_tag}\n")


if __name__ == '__main__':
    set_log_debug_level()
    ## load_genome_annotation_config(True)
    # save_done_file('.')
