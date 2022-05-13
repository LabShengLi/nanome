#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_train.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome
"""
Train model on data using RF or XGBoost
"""
import argparse
import os.path
import sys
from collections import defaultdict
from functools import reduce
from warnings import simplefilter

import math
from tqdm import tqdm

from nanome.common.eval_common import freq_to_label
from nanome.xgboost.xgboost_common import TRUTH_LABEL_COLUMN, default_xgboost_params, gridcv_xgboost_params, \
    gridcv_rf_params, SITES_COLUMN_LIST, default_rf_params

simplefilter(action='ignore', category=FutureWarning)
simplefilter(action='ignore', category=UserWarning)

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, roc_auc_score, precision_score, recall_score, \
    classification_report
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split, StratifiedKFold
from xgboost import XGBClassifier

from nanome.common.global_config import set_log_debug_level, set_log_info_level, logger
from nanome.common.global_settings import NANOME_VERSION, EPSLONG


def report_performance(y_test, y_pred, y_score, toolname=None):
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_score)
    logger.info(
        f"tool={toolname} accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}, #Pred={len(y_test):,}")
    return accuracy, precision, recall, f1, roc_auc


def evaluation_on_na_test(clf, nadf, x_col_list=None, y_col_name=None):
    """
    Evaluation on NA predictions
    Args:
        clf:
        nadf:
        x_col_list:
        y_col_name:

    Returns:

    """
    nadf.reset_index(drop=True, inplace=True)
    logger.debug(f"Eval on NAs   Total predictions:{len(nadf):,}")

    numSites = len(nadf.drop_duplicates(subset=SITES_COLUMN_LIST))
    logger.debug(f"Eval on NAs   Total CpGs:{numSites:,}")
    X_test = nadf[x_col_list]
    y_test = nadf[y_col_name]

    y_pred = pd.DataFrame(clf.predict(X_test))[0]
    y_pred_prob = pd.DataFrame(clf.predict_proba(X_test))[1]

    dataset = []
    for tool in x_col_list:
        mask = nadf[tool].notna()
        tool_notna_num = len(nadf[mask].drop_duplicates(subset=SITES_COLUMN_LIST))
        y_tool_score = nadf[mask][tool]
        y_tool_pred = y_tool_score.apply(lambda x: 1 if x > 0 else 0)
        logger.info(f"For non-NA of {tool} in NA_DF2, total pred={len(y_test[mask]):,}, sites={tool_notna_num:,}")
        accuracy, precision, recall, f1, roc_auc = report_performance(y_test[mask], y_tool_pred, y_tool_score,
                                                                      toolname=tool)
        ret = {'Tool': tool, 'Accuracy': accuracy, 'Precision': precision,
               'Recall': recall, 'F1': f1, 'ROC_AUC': roc_auc,
               '#Bases': tool_notna_num, '#Pred': len(y_tool_score)
               }
        dataset.append(ret)

        accuracy, precision, recall, f1, roc_auc = report_performance(y_test[mask], y_pred[mask], y_pred_prob[mask],
                                                                      toolname=args.model_type)
        ret = {'Tool': args.model_type, 'Accuracy': accuracy, 'Precision': precision,
               'Recall': recall, 'F1': f1, 'ROC_AUC': roc_auc,
               '#Bases': tool_notna_num, '#Pred': len(y_tool_score)
               }
        dataset.append(ret)
    outdf = pd.DataFrame(dataset)
    outfn = args.o + f'_{"_".join(train_tool_list)}' + f'_{args.model_type}_niter{args.niter}_NA_eval.csv'
    outdf.to_csv(outfn)
    logger.info(f"save to {outfn}")


def evaluation_on_test(clf, testdf, x_col_list=None, y_col_name=None, all_tools=None):
    """
    Evaluate trained classifier on non-NA data
    Args:
        clf:
        testdf:
        x_col_list:
        y_col_name:
        all_tools:

    Returns:

    """
    logger.debug("Evaluate on test data...")
    testdf.reset_index(drop=True, inplace=True)
    logger.debug(f"Total predictions:{len(testdf):,}")

    numSites = len(testdf.drop_duplicates(subset=SITES_COLUMN_LIST))
    logger.debug(f"Total CpGs:{numSites:,}")
    X_test = testdf[x_col_list]
    y_test = testdf[y_col_name]

    ## make prediction on test data
    y_pred = clf.predict(X_test)
    y_pred_prob = pd.DataFrame(clf.predict_proba(X_test))[1]

    dataset = defaultdict(list)
    ## evaluate XGBoost model
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred_prob)

    conf_matrix = confusion_matrix(y_test, y_pred)
    logger.info(f"\nconf_matrix=\n{conf_matrix}")
    class_report = classification_report(y_test, y_pred)
    logger.info(f"\nclass_report=\n{class_report}")

    dataset['Tool'].append(args.model_type)
    dataset['Accuracy'].append(accuracy)
    dataset['Precision'].append(precision)
    dataset['Recall'].append(recall)
    dataset['F1'].append(f1)
    dataset['ROC_AUC'].append(roc_auc)
    dataset['#Bases'].append(numSites)
    dataset['#Pred'].append(len(y_test))

    logger.info(
        f"tool=NANOME({args.model_type}) accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}, #Pred={len(y_pred):,}, #Bases={numSites:,}")

    ## evaluate for base model performance

    for tool in all_tools:
        y_score_tool = testdf[tool]
        y_pred_tool = testdf[tool].apply(lambda x: 1 if x >= 0 else 0)
        accuracy = accuracy_score(y_test, y_pred_tool)
        precision = precision_score(y_test, y_pred_tool)
        recall = recall_score(y_test, y_pred_tool)
        f1 = f1_score(y_test, y_pred_tool)
        roc_auc = roc_auc_score(y_test, y_score_tool)
        dataset['Tool'].append(tool)
        dataset['Accuracy'].append(accuracy)
        dataset['Precision'].append(precision)
        dataset['Recall'].append(recall)
        dataset['F1'].append(f1)
        dataset['ROC_AUC'].append(roc_auc)
        dataset['#Bases'].append(numSites)
        dataset['#Pred'].append(len(y_pred_tool))

        logger.info(
            f"tool={tool} accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}, #Pred={len(y_pred_tool):,}, #Bases={numSites:,}")
    outdf = pd.DataFrame.from_dict(dataset)
    outfn = args.o + f"_{args.model_type}_niter{args.niter}_evaluation_on_test_across_tools.csv"
    outdf.to_csv(outfn, index=False)
    print(f"save to {outfn}")


def train_classifier_model(datadf, nadf=None, train_tool_list=None):
    """
    Train xgboost/RF model, tune the params
    Args:
        datadf:

    Returns:

    """
    ## read level to site level df
    sitedf = datadf[SITES_COLUMN_LIST + [TRUTH_LABEL_COLUMN]].drop_duplicates()
    sitedf.reset_index(inplace=True, drop=True)

    if not math.isclose(args.train_size, 1.0):  # not equal 1.0
        ## split sites into train and test
        siteX = sitedf.loc[:, SITES_COLUMN_LIST]
        siteY = sitedf.loc[:, TRUTH_LABEL_COLUMN]
        X_train_site, X_test_site, y_train_site, y_test_site = \
            train_test_split(siteX, siteY, test_size=1 - args.train_size, random_state=args.random_state,
                             stratify=siteY)

        train_test_array = [len(y_train_site), len(y_test_site)]
        logger.debug(
            f"""Split, site level report:
            Train:Test={train_test_array / np.sum(train_test_array)}

            Train data:Sites={len(y_train_site):,}
            class_distribution=\n{y_train_site.value_counts()}
            class_freq=\n{y_train_site.value_counts(normalize=True)}

            Test data:Sites={len(y_test_site):,}
            class_distribution=\n{y_test_site.value_counts()}
            class_freq=\n{y_test_site.value_counts(normalize=True)}
            """)

        ## select predictions (read) based on sites split
        traindf = datadf.merge(X_train_site, on=SITES_COLUMN_LIST, how='inner')
        testdf = datadf.merge(X_test_site, on=SITES_COLUMN_LIST, how='inner')

        ## general info of train and test for non-NA data
        logger.info(general_info_df(traindf, TRUTH_LABEL_COLUMN, 'Train DF (non-NA data)'))
        logger.info(general_info_df(testdf, TRUTH_LABEL_COLUMN, 'Test DF (non-NA data)'))

    else:
        traindf = datadf
        testdf = None

    if nadf is not None:
        sitenadf = nadf[SITES_COLUMN_LIST + [TRUTH_LABEL_COLUMN]].drop_duplicates()
        sitenadf.reset_index(inplace=True, drop=True)
        logger.debug(
            f"NADF  Sites={len(sitenadf):,}, class_distribution=\n{sitenadf[TRUTH_LABEL_COLUMN].value_counts()}")
        if not math.isclose(args.train_size, 1.0):
            ## split sites into train and test
            sitenaX = sitenadf.loc[:, SITES_COLUMN_LIST]
            sitenaY = sitenadf.loc[:, TRUTH_LABEL_COLUMN]
            naX_train_site, naX_test_site, nay_train_site, nay_test_site = \
                train_test_split(sitenaX, sitenaY, test_size=1 - args.train_size, random_state=args.random_state,
                                 stratify=sitenaY)
            logger.debug(
                f"NA: Total CpGs={len(sitenadf):,}, train CpGs={len(naX_train_site):,}, test CpGs={len(naX_test_site):,}")
            train_nadf = nadf.merge(naX_train_site, on=SITES_COLUMN_LIST, how='inner')
            test_nadf = nadf.merge(naX_test_site, on=SITES_COLUMN_LIST, how='inner')
            ## general info of train and test for any-NA data
            logger.info(general_info_df(train_nadf, TRUTH_LABEL_COLUMN, 'Train na_DF (any-NA data)'))
            logger.info(general_info_df(test_nadf, TRUTH_LABEL_COLUMN, 'Test na_DF (any-NA data)'))
        else:
            train_nadf = nadf
            test_nadf = None

        traindf = pd.concat([traindf, train_nadf])
        # testdf = pd.concat([testdf, test_nadf])

    ## save train and test data
    if args.gen_data is not None:
        outdir = os.path.dirname(args.o)

        outfn = os.path.join(outdir,
                             f"{args.gen_data}_train_data_trainsize{args.train_size:.2f}_{args.model_type}_data.csv.gz")
        traindf.to_csv(outfn, index=False)
        logger.info(f"save to {outfn}")

        outfn = os.path.join(outdir,
                             f"{args.gen_data}_test_data_trainsize{args.train_size:.2f}_{args.model_type}_data.csv.gz")
        testdf.to_csv(outfn, index=False)
        logger.info(f"save to {outfn}")

    traindf.reset_index(drop=True, inplace=True)

    if testdf is not None:
        testdf.reset_index(drop=True, inplace=True)

    X_train = traindf[train_tool_list]
    y_train = traindf[TRUTH_LABEL_COLUMN]

    logger.debug(
        f"""Split, read level report:
    Train data:Reads_pred={len(y_train):,}
    class_distribution=\n{y_train.value_counts()}
    class_freq=\n{y_train.value_counts(normalize=True)}
    """)

    ## min max scaller for RF model
    # X_train = MinMaxScaler().fit_transform(X_train)

    ## Training process
    if args.model_type.startswith('XGBoost'):  # ['XGBoost', 'XGBoost_NA', 'XGBoostNA2T', 'XGBoostNA3T', ]
        ## train model using CV and search best params
        default_xgboost_params.update({'random_state': args.random_state})
        clf_model = XGBClassifier(**default_xgboost_params)
        gridcv_params = gridcv_xgboost_params
        logger.debug(f"\n\nDefault params={default_xgboost_params}\n\nSearch gridcv_params={gridcv_params}")
    elif args.model_type == 'RF':
        default_rf_params.update({'random_state': args.random_state})
        clf_model = RandomForestClassifier(**default_rf_params)
        gridcv_params = gridcv_rf_params
        logger.debug(f"\n\nDefault params={default_rf_params}\n\nSearch gridcv_params={gridcv_params}")
    else:
        raise Exception(f"args.model_type={args.model_type} is not supported")

    clf = RandomizedSearchCV(clf_model, gridcv_params, n_iter=args.niter,
                             random_state=args.random_state, n_jobs=args.processors,
                             cv=StratifiedKFold(n_splits=args.cv, shuffle=True,
                                                random_state=args.random_state),
                             scoring='f1', verbose=10, refit=True,
                             return_train_score=True)
    clf.fit(X_train, y_train)
    sys.stdout.flush()
    sys.stderr.flush()
    logger.info(f"model_type=NANOME({args.model_type})")
    logger.info(f"best_params={clf.best_params_}")
    logger.info(f"best_score={clf.best_score_}")

    ## save cv results
    perf_df = pd.DataFrame(clf.cv_results_)
    outfn = args.o + f'_{"_".join(train_tool_list)}' + f'_{args.model_type}_niter{args.niter}_cv_results.xlsx'
    perf_df.to_excel(outfn)
    logger.info(f"save to {outfn}")

    ## save model
    outfn = args.o + f'_{"_".join(train_tool_list)}' + f'_{args.model_type}_niter{args.niter}_model.pkl'
    joblib.dump(clf, outfn)
    logger.info(f"save model to {outfn}")

    ## Evaluate on overlapped predictions for CpGs by all tools
    if testdf is not None:
        testdf.dropna(inplace=True)
        testdf.reset_index(drop=True, inplace=True)
        logger.debug(f"Total overlaped predictions: {len(testdf):,}")
        y_test = testdf[TRUTH_LABEL_COLUMN]

        logger.debug(
            f"""Split, read level report:
                Test data:Reads_pred={len(y_test):,}
                class_distribution=\n{y_test.value_counts()}
                class_freq=\n{y_test.value_counts(normalize=True)}
                """)

    ## evaluate model
    if clf is not None and testdf is not None:
        evaluation_on_test(clf, testdf, x_col_list=train_tool_list, y_col_name=TRUTH_LABEL_COLUMN, all_tools=all_tools)
        if nadf is not None and test_nadf is not None:
            evaluation_on_na_test(clf, test_nadf, x_col_list=train_tool_list, y_col_name=TRUTH_LABEL_COLUMN)


def general_info_df(df, label, dfname=None):
    """
    Output the general info of pred df
    Args:
        df:
        label:

    Returns:

    """
    outstr = '\n\n###### General info of data ######\n'
    if dfname is not None:
        outstr += f"dataset name: {dfname}\n\n"
    outstr += '------------------\n'
    outstr += f'Total predictions: {len(df):,}\n\n'
    outstr += f'Predictions distribution:\n{df[label].value_counts()}\n{df[label].value_counts(normalize=True)}\n\n'
    outstr += '------------------\n'
    sitedf = df[["Chr", "Pos", "Strand", TRUTH_LABEL_COLUMN]].drop_duplicates()
    outstr += f'Total sites: {len(sitedf):,}\n\n'
    outstr += f'Sites distribution:\n{sitedf[label].value_counts()}\n{sitedf[label].value_counts(normalize=True)}\n\n'
    outstr += '####################################\n\n'
    return outstr


def parse_arguments():
    parser = argparse.ArgumentParser(prog='xgboost_train (NANOME)', description='XGBoost train on data')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', nargs='+', help='input data (can be multiple files by chrs) for training', required=True)
    parser.add_argument('-t', nargs='+',
                        help='trained on specific tools, such as megalodon, deepsignal, nanopolish, and guppy',
                        required=True)
    parser.add_argument('-o', type=str, required=True,
                        help='output trained model file name, suggest suffixed with .pkl')
    parser.add_argument('-m', type=str, default=None,
                        help='import XGBoost model from file')
    parser.add_argument('--santity-meteore-fn', type=str, default=None,
                        help='sanity check old meteore results')
    parser.add_argument('--model-type', type=str, default='XGBoost',
                        help='model type for training, can be XGBoost, or RF, or XGBoost_NA')
    parser.add_argument('--bsseq-cov', type=int, default=5,
                        help='coverage cutoff for BS-seq data, default is 5')
    parser.add_argument('--random-state', type=int, default=42,
                        help='random state, default is 42')
    parser.add_argument('--cv', type=int, default=5,
                        help='cv folds, default is 5')
    parser.add_argument('--niter', type=int, default=30,
                        help='num of iterations for RandomeSearchCV, default is 30')
    parser.add_argument('--processors', type=int, default=8,
                        help='num of processors, default is 8')
    parser.add_argument('--train-size', type=float, default=0.5,
                        help='test data ratio: 0.0-1.0, default is 0.5')
    parser.add_argument('--fully-meth-threshold', type=float, default=1.0,
                        help='fully methylated threshold (e.g., 0.9), default is 1.0')
    parser.add_argument('--gen-data', type=str, default=None,
                        help='generate train and test data if specified its name, such as APL, NA12878')
    parser.add_argument('--eval-only', help="if only output evaluation results, not training", action='store_true')
    parser.add_argument('--test', type=int, help="if only test on some number of rows", default=None)
    parser.add_argument('--save-df', help="if save the data into dataframe", action='store_true')
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

    train_tool_list = list(args.t)
    logger.debug(f"train_tool_list={train_tool_list}")

    ## ID,Chr,Pos,Strand,nanopolish,megalodon,deepsignal,guppy,meteore,Freq,Coverage
    dtype = {'ID': str, 'Chr': str, 'Pos': np.int64, 'Strand': str, 'nanopolish': float,
             'megalodon': float, 'deepsignal': float, 'guppy': float, 'meteore': float,
             'Freq': float, 'Coverage': float}

    if math.isclose(args.train_size, 1.0):
        logger.info(f"Full data set is used for training/cross-validation the model")

    datadf_list = []
    nadatadf_list = []
    for infn in tqdm(args.i):
        datadf1 = pd.read_csv(infn, dtype=dtype, header=0, index_col=False, nrows=args.test)
        datadf1['Pos'] = datadf1['Pos'].astype(np.int64)
        # drop guppy due to there is not clear read level outputs for it
        datadf1.drop('guppy', axis=1, inplace=True)
        all_tools = list(datadf1.columns[4:datadf1.columns.get_loc('Freq')])
        datadf1.drop_duplicates(subset=["ID", "Chr", "Pos", "Strand"], inplace=True)
        datadf1.dropna(subset=["Freq", "Coverage"], inplace=True)
        datadf1 = datadf1[datadf1['Coverage'] >= args.bsseq_cov]
        datadf1 = datadf1[(datadf1['Freq'] <= EPSLONG) | (datadf1['Freq'] >= args.fully_meth_threshold - EPSLONG)]
        ## add label for distribution report
        datadf1[TRUTH_LABEL_COLUMN] = datadf1['Freq'].apply(freq_to_label).astype(int)
        datadf1.reset_index(drop=True, inplace=True)

        ## non-NA train data for all tools for evaluation
        datadf2 = datadf1.dropna(subset=all_tools)
        datadf_list.append(datadf2)

        ## NA train data at any one or more tools
        mask_list = [datadf1[tool].isna() for tool in args.t]
        mask_na_index = reduce(lambda left, right: (left | right), mask_list)
        datadf3 = datadf1[mask_na_index]
        nadatadf_list.append(datadf3)

    logger.debug(f"all_tools={all_tools}")
    outdir = os.path.dirname(args.o)
    os.makedirs(outdir, exist_ok=True)

    datadf = pd.concat(datadf_list)
    datadf.reset_index(drop=True, inplace=True)
    datadf_list = None  # save memory
    logger.debug(f"DF1 pred (non-NA)={len(datadf):,}")
    if args.save_df:
        outfn = args.o + f'_joined_nonNA_df.csv.gz'
        datadf.to_csv(outfn, index=False)
        logger.debug(f"save to {outfn}")

    nadatadf = pd.concat(nadatadf_list)
    nadatadf.reset_index(drop=True, inplace=True)
    nadatadf_list = None  # save memory
    logger.debug(f"DF2 pred (any-NA)={len(nadatadf):,}")
    if args.save_df:
        outfn = args.o + f'_joined_NA_df.csv.gz'
        nadatadf.to_csv(outfn, index=False)
        logger.debug(f"save to {outfn}")

    ## Output read-level, site-level distributions for DF1 and DF2
    logger.info(general_info_df(datadf, TRUTH_LABEL_COLUMN, 'DF1 (non-NA data)'))
    logger.info(general_info_df(nadatadf, TRUTH_LABEL_COLUMN, 'DF2 (any-NA data)'))

    train_classifier_model(datadf, nadatadf, train_tool_list=args.t)

    logger.info(f"### train model={args.model_type} on train_tools={args.t} program DONE")
