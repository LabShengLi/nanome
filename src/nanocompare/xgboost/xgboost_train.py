#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_train.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome
"""
Train model on data using RF or XGBoost
"""
import argparse
import os.path
import sys
from collections import defaultdict
from warnings import simplefilter

import math

from nanocompare.eval_common import freq_to_label
from nanocompare.xgboost.xgboost_common import TRUTH_LABEL_COLUMN, default_xgboost_params, gridcv_xgboost_params, \
    gridcv_rf_params, SITES_COLUMN_LIST, default_rf_params, meteore_deepsignal_megalodon_model, load_meteore_model

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

from nanocompare.global_config import set_log_debug_level, set_log_info_level, logger
from nanocompare.global_settings import NANOME_VERSION, EPSLONG


def report_performance(y_test, y_pred, y_score, toolname=None):
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_score)
    logger.info(
        f"tool={toolname} accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}")
    return accuracy, precision, recall, f1, roc_auc


def evaluation_on_na_test(clf, nadf, x_col_list=None, y_col_name=None):
    nadf.reset_index(drop=True, inplace=True)
    logger.debug(f"Eval on NAs   Total predictions:{len(nadf):,}")
    logger.debug(f"Eval on NAs   Total CpGs:{len(nadf.drop_duplicates(subset=SITES_COLUMN_LIST)):,}")
    X_test = nadf[x_col_list]
    y_test = nadf[y_col_name]

    y_pred = pd.DataFrame(clf.predict(X_test))[0]
    y_pred_prob = pd.DataFrame(clf.predict_proba(X_test))[1]

    for tool in x_col_list:
        mask = nadf[tool].notna()
        y_tool_score = nadf[mask][tool]
        y_tool_pred = y_tool_score.apply(lambda x: 1 if x > 0 else 0)
        logger.info(f"Total subsets={len(y_test[mask]):,}")
        report_performance(y_test[mask], y_tool_pred, y_tool_score, toolname=tool)
        report_performance(y_test[mask], y_pred[mask], y_pred_prob[mask], toolname=args.model_type)


def evaluation_on_test(clf, testdf, x_col_list=None, y_col_name=None, all_tools=None):
    """
    Evaluate trained classifier on test dataframe
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
    logger.debug(f"Total CpGs:{len(testdf.drop_duplicates(subset=SITES_COLUMN_LIST)):,}")
    X_test = testdf[x_col_list]
    y_test = testdf[y_col_name]

    ## Scale features before make prediction
    # X_test = MinMaxScaler().fit_transform(X_test)

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

    # cpg_xgboost = len(y_pred)
    # logger.info(f"Number of predictions (read-level): XGBoost={cpg_xgboost:,}")
    # for tool in tool_list:
    #     logger.info(f"Number of predictions (read-level): {tool}={X_test[tool].notna().sum():,}")
    # logger.info(f"Number of predictions (read-level): METEORE={len(X_test.dropna(subset=tool_list)):,}")

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

    logger.info(
        f"tool=NANOME({args.model_type}) accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}")

    ## evaluate for base model performance

    for tool in all_tools + ['METEORE']:
        y_pred_tool = testdf[f"{tool}_pred"]
        if tool != 'METEORE':
            y_score_tool = testdf[tool]
        else:
            y_score_tool = testdf['METEORE_prob']
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

        logger.info(
            f"tool={tool} accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}")
        if tool == 'METEORE':
            # refit and re eval
            meteore_model = load_meteore_model(meteore_deepsignal_megalodon_model)
            X_test_scaled = testdf[['deepsignal_scale', 'megalodon_scale']]
            tmp_y_pred = meteore_model.predict(X_test_scaled)
            tmp_y_pred_prob = pd.DataFrame(meteore_model.predict_proba(X_test_scaled))[1]
            ## evaluate model
            accuracy = accuracy_score(y_test, tmp_y_pred)
            precision = precision_score(y_test, tmp_y_pred)
            recall = recall_score(y_test, tmp_y_pred)
            f1 = f1_score(y_test, tmp_y_pred)
            roc_auc = roc_auc_score(y_test, tmp_y_pred_prob)
            dataset['Tool'].append(f"{tool}_online_pred")
            dataset['Accuracy'].append(accuracy)
            dataset['Precision'].append(precision)
            dataset['Recall'].append(recall)
            dataset['F1'].append(f1)
            dataset['ROC_AUC'].append(roc_auc)
            logger.info(
                f"tool={tool} accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}.  (Refit min-max preprocessing)")
            conf_matrix = confusion_matrix(y_test, tmp_y_pred)
            logger.info(f"\nconf_matrix=\n{conf_matrix}")
    outdf = pd.DataFrame.from_dict(dataset)
    outfn = args.o + f"_{args.model_type}_evaluation_on_test_across_tools.csv"
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
    logger.info(f"Sites={len(sitedf)}, class_distribution=\n{sitedf[TRUTH_LABEL_COLUMN].value_counts()}")

    if not math.isclose(args.train_size, 1.0):
        ## split sites into train and test
        siteX = sitedf.loc[:, SITES_COLUMN_LIST]
        siteY = sitedf.loc[:, TRUTH_LABEL_COLUMN]
        X_train_site, X_test_site, y_train_site, y_test_site = \
            train_test_split(siteX, siteY, test_size=1 - args.train_size, random_state=args.random_state,
                             stratify=siteY)

        train_test_array = [len(y_train_site), len(y_test_site)]
        logger.info(
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
    else:
        traindf = datadf
        testdf = None

    if nadf is not None:
        sitenadf = nadf[SITES_COLUMN_LIST + [TRUTH_LABEL_COLUMN]].drop_duplicates()
        sitenadf.reset_index(inplace=True, drop=True)
        logger.info(
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
            logger.debug(
                f"NA: Total Predictions={len(nadf):,}, train preds={len(train_nadf):,}, test preds={len(test_nadf):,}")
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

    logger.info(
        f"""Split, read level report:
    Train data:Reads_pred={len(y_train):,}
    class_distribution=\n{y_train.value_counts()}
    class_freq=\n{y_train.value_counts(normalize=True)}
    """)

    ## min max scaller for RF model
    # X_train = MinMaxScaler().fit_transform(X_train)

    ## Training process
    if args.model_type in ['XGBoost', 'XGBoost_NA', 'XGBoostNA2T', 'XGBoostNA3T', ]:
        ## train model using CV and search best params
        default_xgboost_params.update({'random_state': args.random_state})
        clf_model = XGBClassifier(**default_xgboost_params)
        gridcv_params = gridcv_xgboost_params
        logger.info(f"\n\nDefault params={default_xgboost_params}\n\nSearch gridcv_params={gridcv_params}")
    elif args.model_type == 'RF':
        default_rf_params.update({'random_state': args.random_state})
        clf_model = RandomForestClassifier(**default_rf_params)
        gridcv_params = gridcv_rf_params
        logger.info(f"\n\nDefault params={default_rf_params}\n\nSearch gridcv_params={gridcv_params}")
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
    outfn = args.o + f'_{"_".join(train_tool_list)}' + f'_{args.model_type}_cv_results.xlsx'
    perf_df.to_excel(outfn)
    logger.info(f"save to {outfn}")

    ## save model
    outfn = args.o + f'_{"_".join(train_tool_list)}' + f'_{args.model_type}_model.pkl'
    joblib.dump(clf, outfn)
    logger.info(f"save model to {outfn}")

    ## Evaluate on overlapped predictions for CpGs by all tools (Top4 + XGBoost/RF)
    if testdf is not None:
        testdf.dropna(inplace=True)
        testdf.reset_index(drop=True, inplace=True)
        logger.debug(f"Total overlaped predictions: {len(testdf):,}")
        y_test = testdf[TRUTH_LABEL_COLUMN]

        logger.info(
            f"""Split, read level report:
                Test data:Reads_pred={len(y_test):,}
                class_distribution=\n{y_test.value_counts()}
                class_freq=\n{y_test.value_counts(normalize=True)}
                """)

    # if args.eval_only and False:
    #     ## import model
    #     if args.m is not None:
    #         logger.debug(f"Model file: {args.m}")
    #         xgboost_cls = joblib.load(args.m)
    #         try:
    #             logger.debug(f"Model info: xgboost_cls={xgboost_cls}")
    #             logger.debug(f"best_params={xgboost_cls.best_params_}")
    #         except:
    #             logger.debug(f"WARNNING: print params encounter problem, may be scikit learn confilicts issue")
    #     else:
    #         raise Exception(f"No model input file specified")
    #
    #     ## evaluation on test data
    #     evaluation_on_test(xgboost_cls, X_test, y_test)
    #     return

    ## evaluate model
    if clf is not None and testdf is not None:
        evaluation_on_test(clf, testdf, x_col_list=train_tool_list, y_col_name=TRUTH_LABEL_COLUMN, all_tools=all_tools)
        if nadf is not None and test_nadf is not None:
            evaluation_on_na_test(clf, test_nadf, x_col_list=train_tool_list, y_col_name=TRUTH_LABEL_COLUMN)


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
    parser.add_argument('--na-data', nargs='+', default=None,
                        help='data sets that include NAs')
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

    tool_list = list(args.t)
    logger.debug(f"tool_list={tool_list}")

    if math.isclose(args.train_size, 1.0):
        logger.info(f"Full data set is used for training/cross-validation the model")

    datadf_list = []
    for infn in args.i:
        datadf1 = pd.read_csv(infn, index_col=False, nrows=args.test)
        datadf1['Pos'] = datadf1['Pos'].astype(np.int64)
        # drop guppy due to there is not clear read level outputs for it
        datadf1.drop('guppy', axis=1, inplace=True)
        datadf1.drop_duplicates(subset=["ID", "Chr", "Pos", "Strand"], inplace=True)
        datadf1.dropna(subset=["Freq", "Coverage"], inplace=True)
        datadf1.dropna(subset=tool_list, inplace=True)
        datadf_list.append(datadf1)
    datadf = pd.concat(datadf_list)
    datadf_list = None  # save memory
    datadf1 = None  # save memory
    logger.debug(f"datadf={datadf}")

    if args.na_data is not None:
        datadf_list = []
        for infn in args.na_data:
            datadf1 = pd.read_csv(infn, index_col=False, nrows=args.test)
            datadf1['Pos'] = datadf1['Pos'].astype(np.int64)
            datadf1.drop('guppy', axis=1, inplace=True)
            datadf1.drop_duplicates(subset=["ID", "Chr", "Pos", "Strand"], inplace=True)
            datadf1.dropna(subset=["Freq", "Coverage"], inplace=True)
            datadf1 = datadf1[datadf1['Coverage'] >= 5]
            datadf1 = datadf1[(datadf1['Freq'] <= EPSLONG) | (datadf1['Freq'] >= 1.0 - EPSLONG)]
            datadf1[TRUTH_LABEL_COLUMN] = datadf1['Freq'].apply(freq_to_label).astype(int)

            datadf1.dropna(subset=tool_list, how='all', inplace=True)
            datadf_list.append(datadf1)
        nadf = pd.concat(datadf_list)
        logger.debug(
            f"NA  nadf={nadf}\n\n total NA preds={len(nadf):,}, nanopolish={nadf['nanopolish'].notna().sum():,}, megalodon={nadf['megalodon'].notna().sum():,}, deepsignal={nadf['deepsignal'].notna().sum():,}")
        nadf_sites = nadf[SITES_COLUMN_LIST].drop_duplicates()
        logger.debug(f"NA  sites={len(nadf_sites):,}")
    else:
        nadf = None

    all_tools = list(datadf.columns[4:datadf.columns.get_loc('Freq')])
    logger.debug(f"all_tools={all_tools}")

    ## Output read-level, site-level distributions
    logger.info(f"Read stats: total={len(datadf):,}, distribution=\n{datadf[TRUTH_LABEL_COLUMN].value_counts()}")

    sitedf = datadf[["Chr", "Pos", "Strand", TRUTH_LABEL_COLUMN]].drop_duplicates()
    logger.info(f"Site stats: total={len(sitedf):,}, distribution=\n{sitedf[TRUTH_LABEL_COLUMN].value_counts()}")
    sitedf = None

    # ## use all-value data, NA for tool1, and tool2 value data
    # logger.debug(f"datadf={datadf}")
    #
    # num_all_value_read_level = len(datadf.dropna(subset=tool_list))
    # num_na_value_for_tool = {}
    # for tool in tool_list:
    #     num_na_value_for_tool[tool] = len(datadf[datadf[tool].isna()])
    # logger.info(
    #     f"NA value stats for read level: all-values={num_all_value_read_level:,}, NA-values={num_na_value_for_tool}")
    # logger.info(
    #     f"Assert all-values+NA-values=total: {num_all_value_read_level + sum(list(num_na_value_for_tool.values()))} == {len(datadf)}")

    train_classifier_model(datadf, nadf, train_tool_list=args.t)

    logger.info(f"### train model={args.model_type} on tools={args.t} program DONE")
