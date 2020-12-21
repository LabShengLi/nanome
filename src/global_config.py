"""
=================================================================
Global variable initializations
=================================================================

Sample usage:
from global_config import *

Main global variables

logger:  log object
dataset_base_dir:
results_base_dir:
pkl_base_dir:   pickle file dir
log_base_dir:   log file dir
pic_base_dir:

All specific used data file name, such as
"""

import datetime
import logging
import os
import sys

import matplotlib.pyplot as plt

init_log_level_prj = logging.INFO

if sys.platform == 'linux':  # Linux helix dir config
    # logger.debug("Running on Linux")
    print("Running on Linux")
    results_dir = "/projects/li-lab/yang/results"  # temp output base
    project_base_dir = "/projects/li-lab/yang/workspace/nano-compare"  # project base
    data_base_dir = os.path.join(project_base_dir, 'data')  # all used data base
    src_base_dir = os.path.join(project_base_dir, 'src')  # source code base
    pkl_base_dir = '/projects/liuya/results/pkl'  # will deprecated later

    # sspairr_path = "/projects/liuya/workspace/R/smooth_scatter/ss_pair.R"

    import platform

    if not platform.node().startswith('helix'):
        # TODO helix do not have this settings
        os.environ['R_HOME'] = "/home/liuya/anaconda3/envs/nmf/lib/R"

elif sys.platform == 'darwin':  # Mac dir config
    # logger.debug("Running on MacOS")
    print("Running on MacOS")
    dataset_base_dir = "/Users/liuya/dataset"
    results_dir = "/Users/liuya/results"
    project_base_dir = "/Users/liuya/PycharmProjects/tcgajax"
    # sspairr_path = "/Users/liuya/Box/dev/R/smooth_scatter/ss_pair.R"
else:
    raise Exception("Unsupported running system.")

# R script for smooth scatter location
# r_smooth_scatter_pair_path = os.path.join(project_base_dir, "rscript", "ss_pair.R")

# ensure TCGA data are located in this base folder for pipeline


# all pictures will be produced at this folder

today_str = datetime.date.today().strftime("%Y-%m-%d")

# ensure pkl is located at this folder, or pkl will be produced at this folder (depends)

# PKL_DIR=RESULTS_DIR/pkl
# LOG_DIR=RESULTS_DIR/log
# PIC_DIR=RESULTS_DIR/YYYY-MM-DD automatic create each day
# pkl_base_dir = os.path.join(results_base_dir, "pkl")
# pkl_base_dir = os.path.join('/projects/liuya/results', "pkl")

log_base_dir = os.path.join(results_dir, "log")
pic_base_dir = os.path.join(results_dir, today_str)

logger = logging.getLogger()  # 不加名称设置root logger
logger.setLevel(logging.DEBUG)


# Ensure dir is created if it is not exist
def ensure_dir(dir_name):
    if not os.path.exists(dir_name):
        os.umask(0)
        os.makedirs(dir_name, exist_ok=True)
        logger.debug("create dir [{}]".format(dir_name))


# create 3 folders if needed
ensure_dir(pic_base_dir)
ensure_dir(log_base_dir)


# Must be after log_base_dir var defined
def init_logging():
    """
    Init both stdout and log file located at LOG_DIR/log-YYYY-MM-DD.txt
    Sample: /Users/liuya/results/log/log-11-30.txt
    :return:
    """
    formatter_log_results = logging.Formatter('%(asctime)s - [%(filename)s:%(lineno)d] - %(levelname)s: %(message)s')  # datefmt='%Y-%m-%d %H:%M:%S'

    # 使用StreamHandler输出到屏幕
    console_handler = logging.StreamHandler()
    console_handler.setLevel(init_log_level_prj)
    console_handler.setFormatter(formatter_log_results)

    # 使用FileHandler输出到文件
    file_handler = logging.FileHandler(os.path.join(log_base_dir, 'log-{}.txt'.format(today_str)))
    file_handler.setLevel(init_log_level_prj)
    file_handler.setFormatter(formatter_log_results)

    # 添加两个Handler
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    # suppress some libs' logging results
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('numexpr').setLevel(logging.WARNING)
    plt.rcParams.update({'figure.max_open_warning': 0})
    return file_handler, console_handler


file_handler, console_handler = init_logging()


# logger.debug("init_logging ok.")


def set_log_debug_level():
    """
    Set debug logger level to file and std out
    :return:
    """
    file_handler.setLevel(logging.DEBUG)
    console_handler.setLevel(logging.DEBUG)


def set_log_info_level():
    """
    Set info logger level to file and std out
    :return:
    """
    file_handler.setLevel(logging.INFO)
    console_handler.setLevel(logging.INFO)


def set_log_error_level():
    file_handler.setLevel(logging.ERROR)
    console_handler.setLevel(logging.ERROR)


#
# # where to find the related files (DataFrame of pickle)
# SURVIVAL_ALL_NAME = "survival_all.pkl"
# SUBTYPE_ALL_NAME = "subtype_info_all.pkl"
# SURVIVAL_SUBTYPE_NAME = "surv_subtype.pkl"
#
# PINFO_NAME = "pinfo_all.pkl"
#
# # globally used column name of pinfo
# PINFO_COL_METH_NAME = "het_methy_value"
# PINFO_COL_MUT_NAME = "het_mut_freq_value"
# PINFO_COL_MUT_LOG_NAME = "het_mut_freq_value_log"
# PINFO_COL_TUMOR_TYPE = "tumor.type"
# PINFO_COL_TUMOR_SUBTYPE_FROM_YUE = "subtype"  # subtype too little, from Yue, same with Haitham
# PINFO_COL_TUMOR_SUBTYPE_FROM_2018_SOURCE = "tumor subtype"  # subtype too little, 'tumor subtype' is too much from
# # https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=info:doi/10.1371/journal.pgen.1007669.s008
#
# PINFO_COL_TUMOR_SUBTYPE = PINFO_COL_TUMOR_SUBTYPE_FROM_YUE
# PINFO_COL_SURV_TIME = "OS.time"
# PINFO_COL_SURV_EVENT = "OS"
#
# PINFO_COL_ABPURITY = 'ABSOLUTE purity'
# PINFO_COL_TCGACUPURITY = 'TCGA purity (curated)'
# PINFO_COL_CPEPURITY = 'CPE'
# PINFO_COL_ESTIMATEPURITY = 'ESTIMATE'
# PINFO_COL_IHCPURITY = 'IHC'
#
# # purity used in pipeline
# PINFO_COL_PURITY = PINFO_COL_CPEPURITY
#
# PINFO_COL_AVG_METH = "avg_meth_val"
#
# PINFO_COL_NEW_PURITY = "Purity"
# PINFO_COL_NEW_SUBTYPE = "Subtype"
#
# # set the threshold of filtering low purity samples
# PURITY_FILTER_THRESHOLD = 0.6
#
# FSMETHOD_RFS = "rfs"
# FSMETHOD_XGBOOST = "xgboost"
# FSMETHOD_CHI2 = "chi2"
#
# FSMETHOD_LIST = [FSMETHOD_RFS, FSMETHOD_XGBOOST, FSMETHOD_CHI2]
#
# NEW_SURVIVAL_CSV = "1-s2.0-S0092867418302290-mmc1.csv"
#
# # from Yue's survival source
# SURVIVAL_MEASURE_DSS = "DSS"
#
# # from journal source
# SURVIVAL_MEASURE_PFI = "PFI"
# SURVIVAL_MEASURE_PFS = "PFS"
# SURVIVAL_MEASURE_DSS1 = "DSS1"
#
# SURVIVAL_MEASURE_LIST = [SURVIVAL_MEASURE_DSS, SURVIVAL_MEASURE_DSS1, SURVIVAL_MEASURE_PFI, SURVIVAL_MEASURE_PFS]
#
#
# class SurvivalMeasure(Enum):
#     """
#     Three types of different array data, Gene53 is the selected CpGs based on Gene53
#     """
#     DSS = 1
#     PFI = 2
#     PFS = 3
#     DSS1 = 4
#     OS = 5
#     EFS = 6
#
#
# CUTOFF_METHOD_MEDIAN = "median"
# CUTOFF_METHOD_STRATIFY = "stratify"
#
# CUTOFF_METHOD_LIST = [CUTOFF_METHOD_MEDIAN, CUTOFF_METHOD_STRATIFY]
#
# PLOT_SCATTER_MARKSIZE = 12
# PLOT_SCATTER_MARKALPHA = 0.6
#
# # read file names convension
# meth_dir_name = "meth"
# nmf_dir_name = "nmf"
#
# METHYLATION450_FN = "meth.txt"
#
# DNAM_AGE_COE_FN = '13059_2013_3156_MOESM3_ESM.csv'
#
# TCGA_METH_FMT = "{dsname}_meth.pkl"
#
# CV_SUMMARY_FMT = "{dsname}_cross_validation_summary_nmf_k_{startk}_{endk}_step_{stepk}"
# CV_GRIDCV_FMT = "{dsname}_cross_validation_gridcv_results_k_{startk}_{endk}_step_{stepk}"
#
# NMF450K_DE_X_PKL_FMT = "450k_de_x_{dsname}.pkl"
# NMF450K_DE_W_K_PKL_FMT = "450k_de_w_{dsname}_k{k}.pkl"
# NMF450K_DE_H_K_PKL_FMT = "450k_de_h_{dsname}_k{k}.pkl"
# NMF450K_DE_MODEL_PKL_FMT = "450k_de_model_{dsname}_k{k}.joblib"
#
#
# # data set need to be analyzed
# # dataset = ["BRCA", "KIRC", "KIRP", "LGG"]
# # dataset_series = pd.Series(dataset)
#
# class DataType(Enum):
#     """
#     Three types of different array data, Gene53 is the selected CpGs based on Gene53
#     """
#     All_450K = 1
#     Gene53 = 2
#     Gene538 = 3
#     CKB = 4
#     OncoKB = 5
#     OncoPathway = 6
#
#
# class HetType(Enum):
#     """
#     Two types of heterogeneity computation, entorpy or epipolymorphism
#     """
#     ENTROPY = 1
#     EPIPOLYMORPHISM = 2
#     PDR = 3


def current_time_str():
    time_str = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f")
    time_str = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")

    return time_str
