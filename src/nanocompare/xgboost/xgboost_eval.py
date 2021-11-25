#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score, f1_score, roc_auc_score, \
    precision_score, recall_score
from tqdm import tqdm

from nanocompare.eval_common import freq_to_label, tool_pred_class_label
from nanocompare.global_config import set_log_debug_level, set_log_info_level, logger
from nanocompare.global_settings import nanome_model_dict, CHUNKSIZE, NANOME_VERSION, EPSLONG
from nanocompare.xgboost.xgboost_common import load_xgboost_model

XGBOOST_MODEL_NAME = 'XGBoost'


def pred_on_data(cls, X_test, y_test):
    ret_conf_matrix = {}
    ret_ypred = {}
    ret_yscore = {}

    ## make prediction on test data
    y_pred1 = cls.predict(X_test)
    y_pred_prob1 = pd.DataFrame(cls.predict_proba(X_test))[1].values
    ret_conf_matrix[XGBOOST_MODEL_NAME] = confusion_matrix(y_test, y_pred1, labels=[0, 1])
    ret_ypred[XGBOOST_MODEL_NAME] = y_pred1
    ret_yscore[XGBOOST_MODEL_NAME] = y_pred_prob1

    ## evaluate for base model performance
    for tool in list(X_test.columns):
        y_score_tool = X_test[tool].fillna(value=0)
        y_pred_tool = y_score_tool.apply(tool_pred_class_label)
        ret_conf_matrix[tool] = confusion_matrix(y_test, y_pred_tool, labels=[0, 1])
        ret_ypred[tool] = y_pred_tool.values
        ret_yscore[tool] = y_score_tool.values
    return ret_conf_matrix, ret_ypred, ret_yscore


def eval_model_on_data(cls, df):
    if len(tool_list) != 2:
        raise Exception(f"Not supported tools={tool_list}")

    ## drop NAs for BS-seq info
    df.dropna(subset=['Freq', 'Coverage'], inplace=True)
    ## drop NAs for both tools
    df.dropna(subset=tool_list, how='all', inplace=True)
    ## filter bsseq cov
    df['Coverage'] = df['Coverage'].astype(int)
    df = df[df['Coverage'] >= args.bsseq_cov]
    ## Apply fully-meth cutoff for BS-seq
    df = df[(df['Freq'] <= EPSLONG) | (df['Freq'] >= args.fully_meth_threshold - EPSLONG)]

    if len(df) < 1:
        return None

    if 'Truth_label' not in df.columns:
        ## create truth label
        df['Truth_label'] = df['Freq'].apply(freq_to_label, args=(args.fully_meth_threshold, EPSLONG))

    X_test = df[tool_list]
    y_test = df['Truth_label']
    conf_matrix3, y_pred3, y_score3 = pred_on_data(cls, X_test, y_test)
    return conf_matrix3, y_pred3, y_score3, y_test.values


def parse_arguments():
    parser = argparse.ArgumentParser(prog='xgboost_predict (NANOME)', description='XGBoost predict for data')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', nargs='+', required=True,
                        help='input data for evaluation')
    parser.add_argument('-m', type=str, required=True,
                        help=f'model file, existing model list: {",".join(list(nanome_model_dict.keys()))}')
    parser.add_argument('--dsname', type=str, required=True,
                        help='dataset name')
    parser.add_argument('-o', type=str, required=True,
                        help='output dir')
    parser.add_argument('--sep', type=str, default='\t',
                        help=f'seperator of input file, default is TAB')
    parser.add_argument('--fully-meth-threshold', type=float, default=1.0,
                        help='fully methylated threshold (e.g., 0.9), default is 1.0')
    parser.add_argument('--bsseq-cov', type=int, default=5,
                        help='coverage cutoff for BS-seq data, default is 5')
    parser.add_argument('--processors', type=int, default=8,
                        help='num of processors, default is 8')
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE,
                        help=f'chunk size for load large data, default is {CHUNKSIZE}')
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

    xgboost_clf = load_xgboost_model(args.m)

    tool_list = None

    ret_list = []
    for infn in args.i:
        logger.info(f"Analyse file:{infn}")
        df = pd.read_csv(infn, index_col=None, sep=args.sep, chunksize=args.chunksize)
        for chunk in tqdm(df):
            # logger.debug(chunk)
            if tool_list is None:
                freq_column_index = chunk.columns.get_loc('Freq')
                tool_list = list(chunk.columns[4:freq_column_index])
            ret = eval_model_on_data(xgboost_clf, chunk)
            ret_list.append(ret)

    ## combine confusion matrix
    combine_conf_matrix = {}
    y_pred = defaultdict(list)
    y_score = defaultdict(list)
    y_truth = []
    all_tools = [XGBOOST_MODEL_NAME] + tool_list
    for ret in ret_list:
        if ret is None:
            continue
        conf_matrix3 = ret[0]
        y_pred3 = ret[1]
        y_score3 = ret[2]
        y_truth1 = ret[3]
        y_truth.extend(y_truth1)
        for tool in conf_matrix3:
            if tool not in combine_conf_matrix:
                combine_conf_matrix[tool] = np.zeros((2, 2), dtype=int)
            combine_conf_matrix[tool] = combine_conf_matrix[tool] + conf_matrix3[tool]
            y_pred[tool].extend(y_pred3[tool])
            y_score[tool].extend(y_score3[tool])

    for k, tool in enumerate(all_tools):
        logger.info(f"""
Tool={tool}, 
conf_matrix={combine_conf_matrix[tool]},
classification_report={classification_report(y_truth, y_pred[tool])},
accuracy={accuracy_score(y_truth, y_pred[tool]):.3f}, F1={f1_score(y_truth, y_pred[tool]):.3f}, roc_auc={roc_auc_score(y_truth, y_score[tool]):.3f}, Precision={precision_score(y_truth, y_pred[tool]):.3f}, Recall={recall_score(y_truth, y_pred[tool]):.3f}
        """)

    logger.info("Done for XGBoost eval")
