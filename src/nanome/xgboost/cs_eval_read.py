#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : cs_eval_read.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Evaluate consensus model at read-level
"""
import argparse
import os.path

import joblib
import math
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix, \
    classification_report
from tqdm import tqdm

from nanome.common.global_config import logger, set_log_debug_level, set_log_info_level
from nanome.common.global_settings import CHUNKSIZE
from nanome.xgboost.ml_common import READS_COLUMN_LIST, prob_to_llr_2, top3_tools, region_order, tool_feature_dict


def report_model_performance(model_name, y_test, y_pred, y_score, region_name="Genome-wide", dsname="NA12878"):
    """
    report performance of tools at region for dataset name
    Args:
        model_name:
        y_test:
        y_pred:
        y_score:
        region_name:
        dsname:

    Returns:

    """
    ## evaluate  model
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    try:
        roc_auc = roc_auc_score(y_test, y_score)
    except:
        roc_auc = None
    conf_matrix = confusion_matrix(y_test, y_pred)
    class_report = classification_report(y_test, y_pred)

    ret = {
        'Dataset': dsname,
        'Tool': model_name,
        'Region': region_name,
        'F1': f1,
        'ROC_AUC': roc_auc,
        'Accuracy': accuracy,
        'Precision': precision,
        'Recall': recall,
        '#Preds': len(y_test),
        '#5C_Preds': (y_test == 0).sum(),
        '#5mC_Preds': (y_test == 1).sum(),
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
    parser.add_argument('--eval-type', type=str, default='read-level',
                        help='evaluation type, read-level or site-level')
    parser.add_argument('--model-base-dir', type=str, default='.',
                        help="model file's base dir")
    parser.add_argument('--test-lines', type=int, default=None,
                        help='test top N rows, such as 10000, default is None')
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE,
                        help=f'chunk size for load large data, default is {CHUNKSIZE}')
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
    llr2_df_list = []
    for infn in tqdm(args.i[:]):
        # logger.debug(f"Processing input: {infn}")
        df_iter = pd.read_csv(infn, sep='\t', header=0, index_col=False, iterator=True,
                              chunksize=args.chunksize, nrows=args.test_lines)
        df_list = []
        for chunck_df in df_iter:
            df1 = chunck_df[
                READS_COLUMN_LIST + top3_tools + ['k_mer'] +
                ['Freq', 'Coverage', 'Label', 'Region']]
            df_list.append(df1)
        df = pd.concat(df_list)

        df.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)
        df = df[df['Coverage'] >= args.bs_cov]
        df.dropna(subset=top3_tools, inplace=True, how='any')
        if args.force_llr2:  # convert ln to log2
            df['megalodon'] = df['megalodon'] / math.log(2)
        df.reset_index(drop=True, inplace=True)
        # logger.debug(f"df={df.shape}")

        llr2_df1 = df[top3_tools + ['Label', 'Region']].copy()

        seqFeatures = df['k_mer'].apply(lambda x: pd.Series(list(x))).copy()
        seqFeatures.columns = [f"DNASeq_{k}" for k in range(len(seqFeatures.columns))]
        # logger.debug(seqFeatures.shape)

        for model_name in model_list:
            mm = model_list[model_name]

            if '_basic' in model_name:
                X1 = df[tool_feature_dict['basic']]
            elif '_megalodon_deepsignal' in model_name:
                X1 = df[tool_feature_dict['megalodon_deepsignal']]
            else:
                raise Exception(f"Not support model_name={model_name}")

            if model_name.endswith('_seq'):
                X2 = seqFeatures
            else:
                X2 = None

            if X2 is not None:
                X12 = pd.concat([X1, X2], axis=1)
            else:
                X12 = X1

            y_score = pd.DataFrame(mm.predict_proba(X12), index=X12.index)[1]
            llr2_df1[model_name] = y_score.apply(prob_to_llr_2)
        # logger.debug(llr2_df1)
        llr2_df_list.append(llr2_df1)
    llr2_df = pd.concat(llr2_df_list)
    llr2_df.reset_index(drop=True, inplace=True)
    # logger.debug(llr2_df)
    llr2_df.info()

    ## evaluate read-level performance
    dataset = []
    for tool in top3_tools + args.model_name:
        # logger.debug(tool)
        # logger.debug(llr2_df['Region'].unique())
        for region in [None] + list(llr2_df['Region'].unique()):
            logger.debug(region)
            y_true = llr2_df['Label']
            y_pred = llr2_df[tool].apply(lambda x: 1 if x >= 0 else 0)
            y_score = llr2_df[tool]
            if region is not None:
                region_index = (llr2_df['Region'] == region)
                y_true = y_true[region_index]
                y_pred = y_pred[region_index]
                y_score = y_score[region_index]
            ret = report_model_performance(tool, y_true, y_pred, y_score,
                                           region_name=region if region is not None else 'Genome-wide',
                                           dsname=args.dsname)
            logger.debug(ret)
            dataset.append(ret)

    ret_df = pd.DataFrame.from_dict(dataset)
    ret_df['Region'] = pd.Categorical(ret_df['Region'],
                                      categories=region_order,
                                      ordered=True)
    ret_df.sort_values(by=['Dataset', 'Region', 'F1'], ascending=[True, True, False], inplace=True)
    logger.debug(ret_df)

    outfn = os.path.join(args.o, f"Consensus_read_level_eval_{args.dsname}.csv")
    ret_df.to_csv(outfn, index=None)
    logger.info(f"## consensus read-level eval DONE")
