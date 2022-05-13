#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_predict.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Predict NANOME consensus results
"""

import argparse
import numpy as np
import os.path
import sys
from functools import reduce

import joblib
import pandas as pd
from tqdm import tqdm

from nanome.common.eval_common import load_tool_read_level_unified_as_df
from nanome.common.global_config import set_log_debug_level, set_log_info_level, logger
from nanome.common.global_settings import CHUNKSIZE, NANOME_VERSION
from nanome.xgboost.xgboost_common import SITES_COLUMN_LIST, READS_COLUMN_LIST, nanome_model_dict, \
    xgboost_mode_base_dir, nanome_model_tool_list_dict


def parse_arguments():
    parser = argparse.ArgumentParser(prog='xgboost_predict (NANOME)', description='XGBoost predict for data')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', nargs='+', required=True,
                        help='input data for predicting')
    parser.add_argument('-m', type=str, required=True,
                        help=f'model file, existing model list: {",".join(list(nanome_model_dict.keys()))}')
    parser.add_argument('--dsname', type=str, required=True,
                        help='dataset name')
    parser.add_argument('-o', type=str, required=True,
                        help='output file name')
    parser.add_argument('-t', nargs='+', help='tools used for prediction, default is: [megalodon, deepsignal]',
                        default=['megalodon', 'deepsignal'])
    parser.add_argument('--random-state', type=int, default=42,
                        help='random state, default is 42')
    parser.add_argument('--processors', type=int, default=8,
                        help='num of processors, default is 8')
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE,
                        help=f'chunk size for load large data, default is {CHUNKSIZE}')
    parser.add_argument('--contain-na', help="if make prediction on NA values", action='store_true')
    parser.add_argument('--tsv-input', help="if input is tsv for tools' read-level format, or else is combined input",
                        action='store_true')
    parser.add_argument('--chrs', nargs='+', help='chromosomes used', default=None)
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()
    logger.debug(f"args={args}")

    import sklearn;
    sklearn.show_versions()

    if args.m in nanome_model_dict:
        infn = os.path.join(xgboost_mode_base_dir, nanome_model_dict[args.m])
    else:
        infn = args.m
    logger.info(f"XGBoost model file: {infn}")
    xgboost_cls = joblib.load(infn)
    try:
        logger.debug(f"Model info: xgboost_cls={xgboost_cls}")
        logger.debug(f"best_params={xgboost_cls.best_params_}")
    except:
        logger.debug(f"WARNING: print params encounter problem")

    if args.tsv_input:
        if len(args.i) != 1:
            raise Exception(f"Not support args.i={args.i}, in tsv_input mode")
        df_tsvfile = pd.read_csv(args.i[0], header=None, sep='\t')
        tool_list = list(df_tsvfile[0])
        dflist = []
        for index, row in tqdm(df_tsvfile.iterrows()):
            if row[1] in ['None', 'NA', 'NULL']:
                empty_frame = {'Chr': pd.Series(dtype='str'), "ID": pd.Series(dtype='str'),
                               "Pos": pd.Series(dtype='int'), "Strand": pd.Series(dtype='str'),
                               row[0]: pd.Series(dtype='float')}
                df = pd.DataFrame(empty_frame)
            else:
                df = load_tool_read_level_unified_as_df(row[1], toolname=row[0], filterChrSet=args.chrs,
                                                        chunksize=args.chunksize)
            # logger.debug(f"df={df}, df_type={df.info()}")
            dflist.append(df)
        if args.contain_na:
            datadf = reduce(
                lambda left, right: pd.merge(left, right, how='outer', on=["ID", "Chr", "Pos", "Strand"]),
                dflist)
        else:
            datadf = reduce(
                lambda left, right: pd.merge(left, right, how='inner', on=["ID", "Chr", "Pos", "Strand"]),
                dflist)
        # release memory
        dflist = None
        datadf.drop_duplicates(subset=["ID", "Chr", "Pos", "Strand"], inplace=True)
        if len(datadf) <= 0:
            key_list = ['Chr', "ID", "Pos", "Strand"] + tool_list + ['Prediction', 'Prob_methylation']
            empty_frame = {keystr: [] for keystr in key_list}
            outdf = pd.DataFrame(empty_frame)
            outdf.to_csv(args.o, sep='\t', index=False)
            logger.info(f"make no predictions")
            logger.error(f"The combined results are empty, for tool_list={tool_list}, fn_list={df_tsvfile[1]}")
            sys.exit(0)

        logger.debug(f"tool_list={tool_list}")
        logger.debug(f"datadf={datadf}")
    else:  # combined input as default input
        dflist = []
        for infn in tqdm(args.i):
            if args.chrs is not None and len(args.chrs) >= 1:
                iter_df = pd.read_csv(infn, header=0, index_col=False, sep=",", iterator=True,
                                      chunksize=args.chunksize)
                datadf1 = pd.concat([chunk[chunk['Chr'].isin(args.chrs)] for chunk in iter_df])
            else:
                datadf1 = pd.read_csv(infn, header=0, sep=',', index_col=False)
            if not args.contain_na:  ## remove any NAs
                datadf1.dropna(subset=args.t, inplace=True)
            else:
                datadf1.dropna(subset=args.t, inplace=True, how='all')
            datadf1.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)
            dflist.append(datadf1)
        datadf = pd.concat(dflist)

        tool_list = list(args.t)
        datadf = datadf[list(datadf.columns[0:4]) + args.t]
        datadf.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)
        logger.debug(f"tool_list={tool_list}")
        logger.debug(f"datadf={datadf}")

    if datadf is not None:
        datadf.reset_index(inplace=True, drop=True)
    else:
        raise Exception(f"datadf can not be None")

    if args.m in nanome_model_tool_list_dict:
        # Check if model's tool is in
        nanome_model_tool_list = nanome_model_tool_list_dict[args.m]
        logger.info(f"XGBoost model {args.m} use tool order: nanome_model_tool_list={nanome_model_tool_list}")
        for t1 in nanome_model_tool_list:
            if t1 not in tool_list:
                raise Exception(
                    f"Input tool file not contain {t1}, tool_list={tool_list}, nanome_model_tool_list={nanome_model_tool_list}")
        # Rearrange tool column as model used order, remove unused order
        datadf = datadf[READS_COLUMN_LIST + nanome_model_tool_list]
        # Update tool_list to model's input tool list
        tool_list = nanome_model_tool_list
    else:
        # The order will be ensured by user
        logger.info(
            f"The order user input is: tool_list={tool_list}, please make sure you provde the same order of XGBoost model used.")

    ## Output read-level, site-level distributions
    logger.debug(f"Read stats: total={len(datadf):,}")

    # logger.debug(f"datadf={datadf}")
    # logger.debug(f"datadf={datadf.info()}")
    sitedf = datadf[SITES_COLUMN_LIST].drop_duplicates()
    logger.debug(f"Site stats: total={len(sitedf):,}")
    sitedf = None

    logger.debug(f"Start predict by XGBoost......")
    predX = datadf.loc[:, tool_list].astype(float)
    prediction = pd.DataFrame(xgboost_cls.predict(predX))
    prediction.rename(columns={0: "Prediction"}, inplace=True)

    prediction_prob = pd.DataFrame(xgboost_cls.predict_proba(predX))[[1]]
    prediction_prob.rename(columns={1: "Prob_methylation"}, inplace=True)

    nanome_df = pd.concat([datadf, prediction, prediction_prob], axis=1)
    nanome_df = nanome_df[READS_COLUMN_LIST + tool_list + ["Prediction", "Prob_methylation"]]
    logger.debug(f"nanome_df={nanome_df}")

    nanome_df.to_csv(args.o, sep='\t', index=False)
    logger.info(f"make predictions:{len(nanome_df):,}")
    logger.info(f"save to {args.o}")
    logger.info(f"### Done for model:{args.m} predict")
