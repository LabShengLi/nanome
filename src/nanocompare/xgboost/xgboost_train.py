#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_train.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome
import argparse
import os.path
import sys

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, roc_auc_score, precision_score, recall_score, \
    classification_report
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split, StratifiedKFold
from xgboost import XGBClassifier

from nanocompare.eval_common import tool_pred_class_label, freq_to_label
from nanocompare.global_config import set_log_debug_level, set_log_info_level, logger
from nanocompare.global_settings import EPSLONG, NANOME_VERSION

default_params = {
    'objective': 'binary:logistic',
    'booster': 'gbtree',
    'eval_metric': 'mlogloss',
    'use_label_encoder': False,
    'nthread': -1,
    'learning_rate': 0.1,
    'n_estimators': 100,
    'max_depth': 6,
    'subsample': 1,
    'colsample_bytree': 1,
    'reg_alpha': 0,
    'reg_lambda': 1,
}

gridcv_search_params = {
    'learning_rate': [0.01, 0.05, 0.1, 0.2],
    'n_estimators': [50, 100, 200, 500],
    'max_depth': [3, 6, 9],
    'subsample': [0.6, 0.8, 1],
    'colsample_bytree': [0.6, 0.8, 1],
    'reg_alpha': [0.5, 0, 2, 5],
    'reg_lambda': [0.5, 1, 2],
}

APL_best_params = {
    'subsample': [0.8], 'reg_lambda': [1],
    'reg_alpha': [0], 'n_estimators': [200],
    'max_depth': [6], 'learning_rate': [0.05],
    'colsample_bytree': [1]
}

NA12878_best_params = {
    'subsample': [0.8], 'reg_lambda': [1],
    'reg_alpha': [0], 'n_estimators': [200],
    'max_depth': [6], 'learning_rate': [0.05],
    'colsample_bytree': [1]
}


# gridcv_search_params = APL_best_params
# gridcv_search_params = NA12878_best_params


def train_xgboost_model(datadf):
    """
    Train xgboost model, tune the params
    Args:
        datadf:

    Returns:

    """
    ## read level to site level df
    sitedf = datadf[["Chr", "Pos", "Strand", "Truth_label"]].drop_duplicates()
    sitedf.reset_index(inplace=True, drop=True)
    logger.info(f"Sites={len(sitedf)}, class_distribution=\n{sitedf['Truth_label'].value_counts()}")

    ## split sites into train and test
    siteX = sitedf.loc[:, ["Chr", "Pos", "Strand"]]
    siteY = sitedf.loc[:, "Truth_label"]
    X_train_site, X_test_site, y_train_site, y_test_site = \
        train_test_split(siteX, siteY, test_size=args.test_size, random_state=args.random_state, stratify=siteY)

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
    traindf = datadf.merge(X_train_site, on=["Chr", "Pos", "Strand"], how='inner')
    X_train = traindf[tool_list]
    y_train = traindf['Truth_label']

    testdf = datadf.merge(X_test_site, on=["Chr", "Pos", "Strand"], how='inner')
    X_test = testdf[tool_list]
    y_test = testdf['Truth_label']

    ## save train and test data
    if args.gen_data is not None:
        outdir = os.path.dirname(args.o)

        outfn = os.path.join(outdir, f"{args.gen_data}_NANOME_train{1 - args.test_size:.2f}_data.csv.gz")
        traindf.to_csv(outfn, index=False)
        logger.info(f"save to {outfn}")

        outfn = os.path.join(outdir, f"{args.gen_data}_NANOME_test{args.test_size:.2f}_data.csv.gz")
        testdf.to_csv(outfn, index=False)
        logger.info(f"save to {outfn}")
        pass

    train_test_array = [len(y_train), len(y_test)]
    logger.info(
        f"""Split, read level report:
            Train:Test={train_test_array / np.sum(train_test_array)}

            Train data:Reads_pred={len(y_train):,}
            class_distribution=\n{y_train.value_counts()}
            class_freq=\n{y_train.value_counts(normalize=True)}

            Test data:Reads_pred={len(y_test):,}
            class_distribution=\n{y_test.value_counts()}
            class_freq=\n{y_test.value_counts(normalize=True)}
            """)

    logger.info(f"\n\nDefault params={default_params}\n\nSearch parameters={gridcv_search_params}")

    ## train model using CV and search best params
    xgb_model = XGBClassifier(**default_params)
    # clf = GridSearchCV(xgb_model, search_cv_params, n_jobs=args.processors,
    #                    cv=StratifiedKFold(n_splits=args.cv, shuffle=True,
    #                                       random_state=args.random_state),
    #                    scoring='f1', verbose=10, refit=True, return_train_score=True)

    clf = RandomizedSearchCV(xgb_model, gridcv_search_params, n_iter=args.niter,
                             random_state=args.random_state, n_jobs=args.processors,
                             cv=StratifiedKFold(n_splits=args.cv, shuffle=True,
                                                random_state=args.random_state),
                             scoring='f1', verbose=10, refit=True,
                             return_train_score=True)
    clf.fit(X_train, y_train)
    sys.stdout.flush()
    sys.stderr.flush()

    logger.info(f"best_params={clf.best_params_}")
    logger.info(f"best_score={clf.best_score_}")

    ## save cv results
    perf_df = pd.DataFrame(clf.cv_results_)
    outfn = args.o.replace('.pkl', '') + '_cv_results.xlsx'
    perf_df.to_excel(outfn)
    logger.info(f"save to {outfn}")

    ## save model
    joblib.dump(clf, args.o)
    logger.info(f"save model to {args.o}")

    ## make prediction on test data
    y_pred = clf.predict(X_test)
    y_pred_prob = pd.DataFrame(clf.predict_proba(X_test))[[1]]

    ## evaluate XGBoost model
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred_prob)

    cpg_xgboost = len(y_pred)
    logger.info(f"Number of predictions (read-level): XGBoost={cpg_xgboost:,}")
    for tool in tool_list:
        logger.info(f"Number of predictions (read-level): {tool}={X_test[tool].notna().sum():,}")
    logger.info(f"Number of predictions (read-level): METEORE={len(X_test.dropna(subset=tool_list)):,}")

    conf_matrix = confusion_matrix(y_test, y_pred)
    logger.info(f"\nconf_matrix=\n{conf_matrix}")

    class_report = classification_report(y_test, y_pred)
    logger.info(f"\nclass_report=\n{class_report}")

    logger.info(
        f"XGBoost accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}")

    ## evaluate for base model performance
    for tool in tool_list:
        y_pred_tool = X_test[tool].apply(tool_pred_class_label)
        y_score_tool = X_test[tool].fillna(value=0)

        ## DeepSignal may be np.nan, since its original results contains NAN
        accuracy = accuracy_score(y_test, y_pred_tool)
        precision = precision_score(y_test, y_pred_tool)
        recall = recall_score(y_test, y_pred_tool)
        f1 = f1_score(y_test, y_pred_tool)
        roc_auc = roc_auc_score(y_test, y_score_tool)
        logger.info(
            f"{tool} accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}")


def parse_arguments():
    parser = argparse.ArgumentParser(prog='xgboost_train (NANOME)', description='XGBoost train on data')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', type=str, required=True,
                        help='input data for training')
    parser.add_argument('-o', type=str, required=True,
                        help='output trained model file name, suggest suffixed with .pkl')
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
    parser.add_argument('--test-size', type=float, default=0.5,
                        help='test data ratio: 0.0-1.0, default is 0.5')
    parser.add_argument('--fully-meth-threshold', type=float, default=1.0,
                        help='fully methylated threshold, default is 1.0')
    parser.add_argument('--gen-data', type=str, default=None,
                        help='generate train and test data if specified its name, such as APL, NA12878')
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

    datadf = pd.read_csv(args.i, sep='\t', index_col=False)
    tool_list = list(datadf.columns[4:-2])
    logger.debug(f"tool_list={tool_list}")

    datadf['Pos'] = datadf['Pos'].astype(np.int64)
    datadf.drop_duplicates(subset=["ID", "Chr", "Pos", "Strand"], inplace=True)
    datadf.dropna(subset=["Freq", "Coverage"], inplace=True)
    datadf.dropna(subset=tool_list, inplace=True, how='all')
    logger.debug(f"datadf={datadf}")

    ## Apply coverage cutoff for BS-seq
    datadf = datadf[datadf['Coverage'] >= args.bsseq_cov]
    logger.debug(f"After cov cutoff={args.bsseq_cov}, len={len(datadf)}")

    ## Apply fully-meth cutoff for BS-seq
    datadf = datadf[(datadf['Freq'] <= EPSLONG) | (datadf['Freq'] >= args.fully_meth_threshold - EPSLONG)]
    logger.debug(
        f"After select fully meth and unmeth, fully_meth_threshold={args.fully_meth_threshold:.2f}, len={len(datadf)}")

    ## Add truth_label
    datadf['Truth_label'] = datadf['Freq'].apply(freq_to_label, args=(args.fully_meth_threshold, EPSLONG))

    ## Output read-level, site-level distributions
    logger.info(f"Read stats: total={len(datadf):,}, distribution=\n{datadf['Truth_label'].value_counts()}")

    sitedf = datadf[["Chr", "Pos", "Strand", "Truth_label"]].drop_duplicates()
    logger.info(f"Site stats: total={len(sitedf):,}, distribution=\n{sitedf['Truth_label'].value_counts()}")
    sitedf = None

    ## use all-value data, NA for tool1, and tool2 value data
    logger.debug(f"datadf={datadf}")

    num_all_value_read_level = len(datadf.dropna(subset=tool_list))
    num_na_value_for_tool = {}
    for tool in tool_list:
        num_na_value_for_tool[tool] = len(datadf[datadf[tool].isna()])
    logger.info(
        f"NA value stats for read level: all-values={num_all_value_read_level:,}, NA-values={num_na_value_for_tool}")
    logger.info(
        f"Assert all-values+NA-values=total: {num_all_value_read_level + sum(list(num_na_value_for_tool.values()))} == {len(datadf)}")

    train_xgboost_model(datadf)

    logger.info(f"### xgboost DONE")
