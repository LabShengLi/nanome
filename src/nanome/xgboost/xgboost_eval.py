#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Evaluate NANOME consensus model on data (0% or 100% BS-seq sites), joined reads
"""

import argparse
import os.path
from collections import defaultdict

import joblib
import math
import pandas as pd
from sklearn.metrics import roc_auc_score, f1_score, recall_score, precision_score, accuracy_score
from tqdm import tqdm

from nanome.common.eval_common import freq_to_label
from nanome.common.global_config import set_log_debug_level, set_log_info_level, logger
from nanome.common.global_settings import CHUNKSIZE, NANOME_VERSION, EPSLONG
from nanome.xgboost.xgboost_common import SITES_COLUMN_LIST, READS_COLUMN_LIST, nanome_model_dict, xgboost_mode_base_dir


def report_read_level_performance(df):
    logger.debug("Evaluate on data...")
    logger.debug(f"Total predictions:{len(df):,}")

    numSites = len(df.drop_duplicates(subset=SITES_COLUMN_LIST))
    logger.debug(f"Total CpGs:{numSites:,}")

    y_truth = df['Freq'].apply(freq_to_label).astype(int)
    dataset = defaultdict(list)
    for tool in all_tools:
        y_score = df[tool]
        y_pred = df[tool].apply(lambda x: 1 if x > 0 else 0)

        accuracy = accuracy_score(y_truth, y_pred)
        precision = precision_score(y_truth, y_pred)
        recall = recall_score(y_truth, y_pred)
        f1 = f1_score(y_truth, y_pred)
        roc_auc = roc_auc_score(y_truth, y_score)
        dataset['Tool'].append(tool)
        dataset['Accuracy'].append(accuracy)
        dataset['Precision'].append(precision)
        dataset['Recall'].append(recall)
        dataset['F1'].append(f1)
        dataset['ROC_AUC'].append(roc_auc)
        dataset['#Bases'].append(numSites)
        dataset['#Pred'].append(len(y_truth))
    outdf = pd.DataFrame.from_dict(dataset)
    outfn = args.o.replace('.csv.gz', '_evaluation.csv')
    outdf.to_csv(outfn)
    logger.info(outdf)
    logger.info(f"save to {outfn}")


def prob_to_log(prob):
    """
    convert probability of 5mC to log-likelyhood
    Args:
        prob:

    Returns:

    """
    return math.log2((prob + EPSLONG) / (1 - prob + EPSLONG))


def parse_arguments():
    parser = argparse.ArgumentParser(prog='xgboost_eval (NANOME)', description='XGBoost eval on data')
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
    parser.add_argument('--fully-meth-threshold', type=float, default=1.0,
                        help='fully methylated threshold (e.g., 0.9), default is 1.0')
    parser.add_argument('--bsseq-cov', type=int, default=5,
                        help='coverage cutoff for BS-seq data, default is 5')
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
    dir_name = os.path.dirname(args.o)

    if args.m in nanome_model_dict:
        infn = os.path.join(xgboost_mode_base_dir, nanome_model_dict[args.m])
    else:
        infn = args.m
    logger.debug(f"Model file: {infn}")
    xgboost_cls = joblib.load(infn)
    try:
        logger.debug(f"Model info: xgboost_cls={xgboost_cls}")
        logger.debug(f"best_params={xgboost_cls.best_params_}")
    except:
        logger.debug(f"WARNNING: print params encounter problem")

    dflist = []
    for infn in tqdm(args.i):
        if args.chrs is not None and len(args.chrs) >= 1:
            iter_df = pd.read_csv(infn, header=0, index_col=False, sep=",", iterator=True,
                                  chunksize=args.chunksize)
            datadf1 = pd.concat([chunk[chunk['Chr'].isin(args.chrs)] for chunk in iter_df])
        else:
            datadf1 = pd.read_csv(infn, header=0, sep=',', index_col=False)

        if 'guppy' in datadf1:
            datadf1.drop('guppy', axis=1, inplace=True)
        datadf1.dropna(subset=['Freq', 'Coverage'], inplace=True)
        datadf1 = datadf1[datadf1['Coverage'] >= args.bsseq_cov]
        datadf1 = datadf1[(datadf1['Freq'] <= EPSLONG) | (datadf1['Freq'] >= args.fully_meth_threshold - EPSLONG)]

        tool_list_input = list(datadf1.columns[4:datadf1.columns.get_loc('Freq')])
        datadf1.dropna(subset=tool_list_input, inplace=True)

        if 'nanopolish' in datadf1:
            datadf1 = datadf1[(datadf1['nanopolish'] >= 2.0) | (datadf1['nanopolish'] <= -2.0)]

        if 'megalodon' in datadf1:
            datadf1 = datadf1[(datadf1['megalodon'] >= 2.0) | (datadf1['megalodon'] <= -2.0)]

        datadf1.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)
        datadf1.reset_index(inplace=True, drop=True)
        logger.debug(f"pred_tool={args.t}")
        logger.debug(f"datadf={datadf1}")

        ## Output read-level, site-level distributions
        logger.debug(f"Read stats: total={len(datadf1):,}")

        sitedf = datadf1[SITES_COLUMN_LIST].drop_duplicates()
        logger.debug(f"Site stats: total={len(sitedf):,}")
        sitedf = None

        logger.debug(f"Start predict by XGBoost......")
        model_name = args.m.replace('NA12878_', '')
        all_tools = tool_list_input + [model_name]

        ## XGBoost predictions
        predX = datadf1.loc[:, args.t].astype(float)
        prediction = pd.DataFrame(xgboost_cls.predict(predX))
        prediction.rename(columns={0: "Prediction"}, inplace=True)
        prediction_prob = pd.DataFrame(xgboost_cls.predict_proba(predX))[[1]]
        prediction_prob.rename(columns={1: "Prob_methylation"}, inplace=True)
        xgboost_score = prediction_prob["Prob_methylation"].apply(prob_to_log)
        xgboost_score.rename(model_name, inplace=True)

        nanome_df1 = pd.concat([datadf1, xgboost_score], axis=1)
        nanome_df1 = nanome_df1[READS_COLUMN_LIST + all_tools + ['Freq', 'Coverage']]
        logger.debug(f"nanome_df1={nanome_df1}")

        outfn = os.path.join(dir_name, os.path.basename(infn).replace('.csv.gz', f'_pred_{model_name}.csv.gz'))
        nanome_df1.to_csv(outfn, index=False)
        logger.info(f"make predictions:{len(nanome_df1):,}")
        logger.info(f"save to {outfn}")
        dflist.append(nanome_df1)

    nanome_df = pd.concat(dflist)
    nanome_df.to_csv(args.o, index=False)
    logger.info(f"save to {args.o}")

    report_read_level_performance(nanome_df)
    logger.info(f"### Done for model:{args.m} predict")
