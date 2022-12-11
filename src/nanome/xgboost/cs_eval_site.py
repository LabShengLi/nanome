#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : cs_eval_read.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Evaluate consensus model
"""
import argparse
import os.path
import warnings

import joblib
import math
import pandas as pd
from scipy import stats
from scipy.stats import PearsonRConstantInputWarning
from sklearn.metrics import mean_squared_error
from tqdm import tqdm

from nanome.common.global_config import logger, set_log_debug_level, set_log_info_level
from nanome.common.global_settings import CHUNKSIZE
from nanome.xgboost.ml_common import READS_COLUMN_LIST, prob_to_llr_2, top3_tools, region_order


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
    parser = argparse.ArgumentParser(prog='cs_eval (NANOME)', description='Consensus model train on data')
    parser.add_argument('-i', nargs='+', required=True,
                        help='input data file')
    parser.add_argument('--dsname', type=str, default="NA12878",
                        help='dataset name, default is NA12878')
    parser.add_argument('--model-name', nargs='+', type=str, default="xgboost",
                        help='model name: rf, xgboost, etc.')
    parser.add_argument('--model-file', nargs='+', type=str, default="xgboost.pkl",
                        help='model file')
    parser.add_argument('-o', type=str, required=True,
                        help='output file dir')
    parser.add_argument('--processors', type=int, default=1,
                        help='number of processors, default is 1')
    parser.add_argument('--bs-cov', type=int, default=5,
                        help='bs-seq coverage cutoff, default is 5')
    parser.add_argument('--tool-cov', type=int, default=1,
                        help='ONT tool coverage cutoff, default is 1')
    parser.add_argument('--eval-type', type=str, default='site-level',
                        help='evaluation type, i.e., site-level')
    parser.add_argument('--model-base-dir', type=str, default='.',
                        help="model file's base dir")
    parser.add_argument('--test-lines', type=int, default=None,
                        help='test top N rows, such as 10000, default is None')
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE,
                        help=f'chunk size for load large data, default is {CHUNKSIZE}')
    parser.add_argument('--save-data', type=str, default=None,
                        help='if save prediction outputs')
    parser.add_argument('--force-llr2', help="if convert megalodon llr to llr2", action='store_true')
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

    ## load models
    model_list = {}
    for k in range(len(args.model_name)):
        model_name = args.model_name[k]
        model_file = args.model_file[k]
        infn = os.path.join(args.model_base_dir, model_file)
        # logger.debug(f"load modelname={model_name}, model file: {infn}")
        model_cls = joblib.load(infn)
        model_list[model_name] = model_cls
    logger.debug(f"num of tools = {len(model_list)}")

    ## predict on each input
    site_df_list = []
    for infn in tqdm(args.i[:]):
        # logger.debug(f"Processing input: {infn}")
        df_iter = pd.read_csv(infn, sep='\t', header=0, index_col=False, iterator=True,
                              chunksize=args.chunksize, nrows=args.test_lines)
        df_list = []
        for chunck_df in df_iter:
            df1 = chunck_df[
                READS_COLUMN_LIST + top3_tools + ['k_mer'] +
                ['Freq', 'Coverage', 'Region']].copy()
            df1.dropna(subset=top3_tools, inplace=True, how='any')
            df1.dropna(subset=['Freq', 'Coverage', 'Region'], inplace=True, how='any')
            df1['k_mer'].fillna('N' * 17, inplace=True)
            # df1.info()
            # logger.debug(df1['k_mer'])
            df_list.append(df1)
        df = pd.concat(df_list)

        df.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)
        df = df[df['Coverage'] >= args.bs_cov]

        if args.force_llr2:  # convert ln to log2
            df['megalodon'] = df['megalodon'] / math.log(2)
        df.reset_index(drop=True, inplace=True)
        # logger.debug(f"df={df.shape}")

        llr2_df1 = df[READS_COLUMN_LIST + top3_tools + ['Freq', 'Coverage', 'Region']].copy()
        top3Features = df[top3_tools]
        # logger.debug(top3Features.shape)
        seqFeatures = df['k_mer'].apply(lambda x: pd.Series(list(x))).copy()
        seqFeatures.columns = [f"DNASeq_{k}" for k in range(len(seqFeatures.columns))]
        # logger.debug(seqFeatures.shape)

        X12 = pd.concat([top3Features, seqFeatures], axis=1)
        # logger.debug(X12.shape)

        for model_name in model_list:
            mm = model_list[model_name]
            if not model_name.endswith('_seq'):
                y_score = pd.DataFrame(mm.predict_proba(top3Features), index=top3Features.index)[1]
            else:
                y_score = pd.DataFrame(mm.predict_proba(X12), index=X12.index)[1]
            # logger.debug(f"y_score.shape={y_score.shape}, model_name={model_name}")

            llr2_df1[model_name] = y_score.apply(prob_to_llr_2)
        # logger.debug(llr2_df1)
        ## convert llr2 to pred
        for tool in top3_tools + args.model_name:
            llr2_df1[tool] = llr2_df1[tool].apply(lambda x: 1 if x >= 0 else 0)
        # logger.debug(llr2_df1)
        # llr2_df1.info()

        ##     df_gp = df[["Chr", "Pos", "Strand", "ID"] + tool_list + ["Freq"]].groupby(by=["Chr", "Pos", "Strand"]).agg(agg_func)

        agg_func = {tool: 'mean' for tool in top3_tools + args.model_name}
        agg_func.update({'Freq': 'first', 'Coverage': 'first', 'Region': 'first', 'ID': 'count'})
        site_df1 = llr2_df1.groupby(by=['Chr', 'Pos', 'Strand']).agg(agg_func)
        site_df1 = site_df1[site_df1['ID'] >= args.tool_cov]
        site_df1.reset_index(drop=False, inplace=True)
        # logger.debug(site_df1)

        site_df_list.append(site_df1)
    site_df = pd.concat(site_df_list)
    site_df.reset_index(drop=True, inplace=True)

    if args.save_data is not None:
        outfn = os.path.join(args.o, args.save_data)
        site_df.to_csv(outfn, index=False)
        logger.info(f"save to {outfn}")
    # logger.debug(site_df)
    # site_df.info()

    ## evaluate site-level performance
    # logger.debug(list(site_df['Region'].unique()))
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
