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

if sys.platform == 'linux':  # Linux HPC dir config
    # print("Running on Linux")
    results_dir = "/projects/li-lab/yang/results"  # temp output base
    project_base_dir = "/projects/li-lab/yang/workspace/nano-compare"  # project base
    data_base_dir = os.path.join(project_base_dir, 'data')  # all used data base
    src_base_dir = os.path.join(project_base_dir, 'src')  # source code base
    pkl_base_dir = '/projects/liuya/results/pkl'  # will deprecated later

    results_base_dir = '/projects/li-lab/Nanopore_compare/result'
    cache_dir = '/fastscratch/liuya/nanocompare/cache_dir' # cache readed object to pkl

    import platform

    if not platform.node().startswith('helix'):
        # TODO helix do not have this settings
        os.environ['R_HOME'] = "/home/liuya/anaconda3/envs/nmf/lib/R"

elif sys.platform == 'darwin':  # Mac dir config
    print("Running on MacOS")
    dataset_base_dir = "/Users/liuya/dataset"
    results_dir = "/Users/liuya/results"
    project_base_dir = "/Users/liuya/PycharmProjects/tcgajax"
    # sspairr_path = "/Users/liuya/Box/dev/R/smooth_scatter/ss_pair.R"
else:
    raise Exception("Unsupported running system.")

today_str = datetime.date.today().strftime("%Y-%m-%d")

log_base_dir = os.path.join(results_dir, "log")
pic_base_dir = os.path.join(results_dir, today_str)

logger = logging.getLogger()
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

formatter_log_results = logging.Formatter('%(asctime)s - [%(filename)s:%(lineno)d] - %(levelname)s: %(message)s')  # datefmt='%Y-%m-%d %H:%M:%S'


def add_logging_file(fn):
    file_handler = logging.FileHandler(fn)
    file_handler.setFormatter(formatter_log_results)
    logger.addHandler(file_handler)
    file_handler.setLevel(logging.DEBUG)


# Must be after log_base_dir var defined
def init_logging():
    """
    Init both stdout and log file located at LOG_DIR/log-YYYY-MM-DD.txt
    Sample: /Users/liuya/results/log/log-11-30.txt
    :return:
    """
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


def current_time_str():
    # time_str = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f")
    time_str = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
    return time_str
