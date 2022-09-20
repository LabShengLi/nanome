#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : cs_agg_site.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Evaluate consensus model
"""
import argparse
import os.path
import warnings

import pandas as pd
from scipy import stats
from scipy.stats import PearsonRConstantInputWarning
from sklearn.metrics import mean_squared_error

from nanome.common.global_config import logger, set_log_debug_level, set_log_info_level
from nanome.common.global_settings import CHUNKSIZE
from nanome.xgboost.ml_common import top3_tools, region_order

SITE_COLUMNS = ['Chr', 'Pos', 'Strand']
READ_COLUMNS = ['ID'] + SITE_COLUMNS

# top3_tools = ['megalodon', 'nanopolish', 'deepsignal']
# region_order = ['Genome-wide', 'Discordant', 'Concordant', 'Singleton', 'Nonsingleton']

dna_seq_order = ['A', 'C', 'G', 'T']

# will be infered later
k_mer_k = 17


# cutoff_llr_tools = {'nanopolish': 2, 'megalodon': math.log(4)}


def report_pcc(tool, y_test, y_pred, region_name="Genome-wide", dsname="NA12878"):
    """
    report performance of tools at region for dataset name
    Args:
        tool:
        y_test:
        y_pred:
        y_score:
        region_name:
        dsname:

    Returns:

    """
    ## evaluate  model
    try:  # too few samples will fail
        # with warnings.catch_warnings(): # not function
        warnings.filterwarnings('ignore', category=PearsonRConstantInputWarning)
        coe, pval = stats.pearsonr(y_test, y_pred)
    except:
        coe, pval = None, None

    mse = mean_squared_error(y_test, y_pred)

    ret = {
        'Dataset': dsname,
        'Tool': tool,
        'Region': region_name,
        'PCC': coe,
        'Pvalue': pval,
        'MSE': mse,
        '#Sites': len(y_test),
    }
    return ret


def parse_arguments():
    parser = argparse.ArgumentParser(prog='cs_agg_site (NANOME)', description='Consensus model site-level eval on data')
    parser.add_argument('-i', nargs='+', required=True,
                        help='input data file')
    parser.add_argument('--dsname', type=str, default="NA12878",
                        help='dataset name, default is NA12878')
    parser.add_argument('--model-name', nargs='+', type=str, default="xgboost",
                        help='model name: rf, xgboost, etc.')
    parser.add_argument('-o', type=str, required=True,
                        help='output file dir')
    parser.add_argument('--bs-cov', type=int, default=5,
                        help='bs-seq coverage cutoff, default is 5')
    parser.add_argument('--tool-cov', type=int, default=1,
                        help='ONT tool coverage cutoff, default is 1')
    parser.add_argument('--eval-type', type=str, default='site-level',
                        help='evaluation type, i.e., site-level')
    parser.add_argument('--test-lines', type=int, default=None,
                        help='test top N rows, such as 10000, default is None')
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE,
                        help=f'chunk size for load large data, default is {CHUNKSIZE}')
    parser.add_argument('--save-data', type=str, default=None,
                        help='if save prediction outputs')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    ## input arguments
    args = parse_arguments()
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()
    logger.debug(f"args={args}")

    dflist = []
    for infn in args.i:
        df = pd.read_csv(infn, index_col=None, header=0, nrows=args.test_lines)
        dflist.append(df)
    site_df = pd.concat(dflist)

    dataset = []
    for tool in top3_tools + args.model_name:
        for region in [None] + list(site_df['Region'].unique()):
            # logger.debug(region)
            y_true = site_df['Freq']
            y_pred = site_df[tool]
            if region is not None:
                region_index = (site_df['Region'] == region)
                y_true = y_true[region_index]
                y_pred = y_pred[region_index]
            ret = report_pcc(tool, y_true, y_pred,
                             region_name=region if region is not None else 'Genome-wide',
                             dsname=args.dsname)
            # logger.debug(ret)
            dataset.append(ret)

    ret_df = pd.DataFrame.from_dict(dataset)
    ret_df['Region'] = pd.Categorical(ret_df['Region'],
                                      categories=region_order,
                                      ordered=True)
    ret_df.sort_values(by=['Dataset', 'Region', 'PCC'], ascending=[True, True, False], inplace=True)
    logger.debug(ret_df)

    outfn = os.path.join(args.o, f"Consensus_site_level_eval_{args.dsname}.csv")
    ret_df.to_csv(outfn, index=False)
    logger.info(f"save to {outfn}")
    logger.info(f"## consensus site-level eval DONE")
