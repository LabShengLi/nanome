#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : cs_train.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Train and evaluate consensus model on different regions, using xgboost, or RF base model, weight learning.

Base model: xgboost, rf
Model name:
    basic (top3 tools),
    basic_w (top3 tools with weight learning),
    basic_w_seq (top3 tools + DNA seq with weight learning)

    megalodon_deepsignal_w
    megalodon_deepsignal_w_seq
"""
import argparse
import os.path
import sys
from collections import defaultdict

import joblib
import math
import numpy as np
import pandas as pd
from sklearn.compose import make_column_transformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix, \
    classification_report
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from xgboost import XGBClassifier

from nanome.common.global_config import logger, set_log_debug_level, set_log_info_level
from nanome.common.global_settings import prob_to_llr2
from nanome.xgboost.ml_common import SITES_COLUMN_LIST, READS_COLUMN_LIST, top3_tools

default_rf_params = {
    'n_jobs': -1
}

gridcv_rf_params = {
    'cls__n_estimators': [5, 10, 20, 40, 80, 100, 150, 200],
    'cls__max_depth': [2, 3, 6, 9, 11, 15, 20],
    'cls__oob_score': [True, False],
    'cls__min_samples_leaf': [1, 5, 10, 50, 100],
    'cls__criterion': ['gini', 'entropy'],
    'cls__max_features': ['auto', 'sqrt', 'log2'],
}

default_xgboost_params = {
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

# gridcv_xgboost_params = {
#     'cls__learning_rate': [0.01, 0.05, 0.1, 0.2, 0.3],
#     'cls__n_estimators': [5, 10, 20, 40, 60],
#     'cls__max_depth': [3, 6, 9, 12],
#     'cls__subsample': [0.4, 0.6, 0.8, 1],
#     'cls__colsample_bytree': [0.4, 0.6, 0.8, 1],
#     'cls__reg_alpha': [0.5, 0, 2, 5],
#     'cls__reg_lambda': [0.5, 1, 2],
# }

gridcv_xgboost_params = {
    'cls__learning_rate': [0.01, 0.05, 0.1, 0.2],
    'cls__n_estimators': [10, 20, 40, 60],
    'cls__max_depth': [3, 6, 9],
    'cls__subsample': [0.6, 0.8, 1],
    'cls__colsample_bytree': [0.4, 0.6, 0.8, 1],
    'cls__reg_alpha': [0.5, 0, 1],
    'cls__reg_lambda': [0.5, 1, 2],
}

tool_list = ['megalodon', 'nanopolish', 'deepsignal', 'nanome3t', 'nanome2t', 'meteore']

exclude_report_tools = ['nanome3t', 'nanome2t', 'meteore']

region_order = ['Genome-wide', 'Discordant', 'Concordant', 'Singleton']

dna_seq_order = ['A', 'C', 'G', 'T']

# will be infered later
k_mer_k = 17

# cutoff_llr_tools = {'nanopolish': 2, 'megalodon': math.log(4),
#                     'xgboost_basic': 2, 'xgboost_basic_w': 2, 'xgboost_basic_w_seq': 2,
#                     'rf_basic': 2, 'rf_basic_w': 2, 'rf_basic_w_seq': 2,
#                     }


def describe_data(df):
    """
    Describe the dataset
    Args:
        df:

    Returns:

    """
    logger.info("#######################################")
    logger.info("########## Dataset description ########")
    ret = pd.concat([df['Label'].value_counts(), df['Label'].value_counts(normalize=True)], axis=1)
    logger.info(f"\n#Pred={len(df):,}, Distribution={ret}\n")

    ret = pd.concat([df['Region'].value_counts(), df['Region'].value_counts(normalize=True)], axis=1)
    logger.info(f"\nRegion Distribution={ret}\n")

    label_vector = df.drop_duplicates(subset=SITES_COLUMN_LIST)['Label']
    ret = pd.concat([label_vector.value_counts(), label_vector.value_counts(normalize=True)], axis=1)
    logger.info(f"\n#Base={len(label_vector):,}, Distribution={ret}\n")
    logger.info("#######################################")


def compute_sample_weights(regionList):
    """
    Compute weights based on region label list
    Args:
        regionList: region label for each sample, len=sample length, type is Series

    Returns:
        sample_weights as Series

    """
    total = len(regionList)
    vc = regionList.value_counts()
    num_class = len(vc)

    class_weight = total / (num_class * vc)
    logger.debug(f"vc={vc}")
    logger.debug(f"class_weight={class_weight}")

    sample_weight = regionList.apply(lambda x: class_weight[x])
    logger.debug(f"sample_weight={sample_weight}")
    return sample_weight


def get_data(infn, model_name="basic", input_tools=top3_tools):
    """
    Input data format:
    <class 'pandas.core.frame.DataFrame'>
    Int64Index: 100000 entries, 0 to 99999
    Data columns (total 22 columns):
     #   Column          Non-Null Count   Dtype
    ---  ------          --------------   -----
     0   Chr             100000 non-null  object
     1   Pos             100000 non-null  int64
     2   Strand          100000 non-null  object
     3   Pos_deprecated  100000 non-null  int64
     4   ID              100000 non-null  object
     5   read_strand     100000 non-null  object
     6   k_mer           100000 non-null  object
     7   signal_means    100000 non-null  object
     8   signal_stds     100000 non-null  object
     9   signal_lens     100000 non-null  object
     10  raw_signals     100000 non-null  object
     11  methy_label     100000 non-null  int64
     12  nanome3t        100000 non-null  float64
     13  nanome2t        100000 non-null  float64
     14  megalodon       100000 non-null  float64
     15  nanopolish      100000 non-null  float64
     16  deepsignal      100000 non-null  float64
     17  meteore         100000 non-null  float64
     18  Freq            100000 non-null  float64
     19  Coverage        100000 non-null  int64
     20  Label           100000 non-null  int64
     21  Region          100000 non-null  object
    dtypes: float64(7), int64(5), object(10)
    memory usage: 17.5+ MB

    Args:
        infn:

    Returns:
        X       tools +/- DNAseq
        y       label
        toolDF  contains 'Region'

    """
    ## decide including usecols
    if "_seq" in model_name.lower():
        usecols = READS_COLUMN_LIST + tool_list + ['k_mer'] + ['Region', 'Label']
    else:
        usecols = READS_COLUMN_LIST + tool_list + ['Region', 'Label']

    if type(infn) == list:
        dflist = []
        for infn1 in infn:
            dflist.append(pd.read_csv(infn1, sep='\t', header=0, nrows=args.test_lines, usecols=usecols))
        df = pd.concat(dflist)
        dflist = None  # free memory
    else:
        df = pd.read_csv(infn, sep='\t', header=0, nrows=args.test_lines, usecols=usecols)
    logger.debug(f"len={len(df):,}")
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna()
    df.reset_index(drop=True, inplace=True)
    logger.debug(f"df={df}")
    logger.debug(f"df.dtypes={df.dtypes}")
    logger.debug(f"len={len(df):,}")

    describe_data(df)

    y = df['Label'].copy()
    toolDF = df[top3_tools + ['Region']].copy()

    if "_seq" in model_name.lower():
        # top3 + DNA seq as feature, i.e., ACCACACCCGGCTAATT
        global k_mer_k
        k_mer_k = len(df['k_mer'].iloc[0])
        logger.debug(f"Detect k_mer_k={k_mer_k}")

        X1 = df[input_tools].copy()
        X2 = df['k_mer'].apply(lambda x: pd.Series(list(x))).copy()
        X2.columns = [f"DNASeq_{k}" for k in range(len(X2.columns))]
        global cat_columns
        cat_columns = list(X2.columns)
        logger.debug(f"X2={X2}")
        X = pd.concat([X1, X2], axis=1)
        logger.debug(f"X1.shape={X1.shape}, X2.shape={X2.shape}, X.shape={X.shape}, X.dtypes={X.dtypes}")
    else:
        X = df[input_tools].copy()

    X.reset_index(drop=True, inplace=True)
    y.reset_index(drop=True, inplace=True)
    toolDF.reset_index(drop=True, inplace=True)
    logger.debug(f"X.shape={X.shape}, y.shape={y.shape}, toolDF.shape={toolDF.shape}")
    logger.debug(f"X={X}, y={y}, toolDF={toolDF}")

    return X, y, toolDF


def train_model(X, y, model_name="basic", region_vector=None, scoring='f1'):
    logger.info(f"Start train, base_model={args.base_model}")

    ## two types of features
    numeric_columns = list(X.select_dtypes(include=[np.number]).columns)
    cat_columns = list(X.select_dtypes(include=['object']).columns)
    logger.debug(f"numeric_columns={numeric_columns}, cat_columns={cat_columns}")

    ## compute weights
    if "_w" in model_name.lower() or "_w_seq" in model_name.lower():
        ## compute sample weights
        sample_weight = compute_sample_weights(region_vector)
    else:
        sample_weight = None

    ## compose models
    simple_imp = SimpleImputer(missing_values=np.nan, strategy='constant')
    ohe = OneHotEncoder(categories=[dna_seq_order] * k_mer_k, handle_unknown='ignore', sparse=False)
    if args.base_model.lower() == 'rf':
        # Make column transformer with one-hot encoder and simple imputer
        column_transform = make_column_transformer(
            (simple_imp, numeric_columns),
            (ohe, cat_columns),
            remainder='passthrough')

        default_rf_params.update({'random_state': args.random_state})
        logger.debug(f"\n\nDefault params={default_rf_params}")
        clf_model = RandomForestClassifier(**default_rf_params)
        gridcv_params = gridcv_rf_params
    elif args.base_model.lower() == 'xgboost':
        # Make column transformer with one-hot encoder and simple imputer
        column_transform = make_column_transformer(
            (ohe, cat_columns),
            remainder='passthrough')

        default_xgboost_params.update({'random_state': args.random_state})
        logger.debug(f"\n\nDefault params={default_xgboost_params}")
        clf_model = XGBClassifier(**default_xgboost_params)
        gridcv_params = gridcv_xgboost_params
    else:
        raise Exception(f"params base_model={args.base_model} not support")
    logger.debug(f"\nSearch gridcv_params={gridcv_params}")

    # Pipeline
    pipe = Pipeline([('prep', column_transform), ('cls', clf_model)])
    logger.debug(f"pipeline={pipe}")

    mm = RandomizedSearchCV(pipe, gridcv_params, n_iter=args.niter,
                            random_state=args.random_state, n_jobs=args.processors,
                            cv=StratifiedKFold(n_splits=args.cv, shuffle=True,
                                               random_state=args.random_state),
                            scoring=scoring, verbose=10, refit=True,
                            return_train_score=True)

    # mm = HalvingRandomSearchCV(pipe, gridcv_params, n_iter=args.niter,
    #                         random_state=args.random_state, n_jobs=args.processors,
    #                         cv=StratifiedKFold(n_splits=args.cv, shuffle=True,
    #                                            random_state=args.random_state),
    #                         scoring='f1', verbose=10, refit=True,
    #                         return_train_score=True)

    if sample_weight is not None:
        logger.info(f"Fit with sample_weight={sample_weight}")
        mm.fit(X, y, cls__sample_weight=list(sample_weight))
    else:
        logger.info(f"Fit without sample weight")
        mm.fit(X, y)
    sys.stdout.flush()
    sys.stderr.flush()
    logger.debug(f"best_params={mm.best_params_}")
    logger.debug(f"best_score={mm.best_score_}")

    ## save cv results
    perf_df = pd.DataFrame(mm.cv_results_)
    outfn = f'NA12878_{train_chrs}_{report_mm_name}_cv_results.xlsx'
    perf_df.to_excel(os.path.join(args.o, outfn))
    logger.info(f"save to {outfn}")

    ## save model
    outfn = f'NA12878_{train_chrs}_{report_mm_name}_model.pkl'
    joblib.dump(mm, os.path.join(args.o, outfn))
    logger.info(f"save model to {outfn}")

    return mm


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

    if args.show_confusion_matrix:
        logger.debug(f"\n\n########################################")
        logger.debug(f"\nreport_mm_name={report_mm_name}: conf_matrix=\n{conf_matrix}")
        logger.debug(f"\nreport_mm_name={report_mm_name}: class_report=\n{class_report}")
        logger.debug(
            f"report_mm_name={report_mm_name}: accuracy={accuracy:.3f}, precision={precision:.3f}, recall={recall:.3f}, f1={f1:.3f}, roc_auc={roc_auc:.3f}, #Pred={len(y_pred):,}")

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


def model_predict(mm, X_test):
    """
    Return prediction class and probabilities [0,1]
    Args:
        mm:
        X_test:

    Returns:

    """
    y_pred = pd.Series(mm.predict(X_test), index=X_test.index)
    y_score = pd.DataFrame(mm.predict_proba(X_test), index=X_test.index)[1]
    return y_pred, y_score


def eval_tools(y_test, toolDF=None, regionList=None, dsname=None):
    # used for aggregate all test data
    ret_pred = defaultdict(list)
    dataset = []

    for region in ['Genome-wide'] + list(regionList.unique()):
        for tool in toolDF.columns:
            y_score = toolDF[tool]
            y_pred = toolDF[tool].apply(lambda x: 1 if x > 0 else 0)

            if region == 'Genome-wide':
                y_score1 = y_score
                y_pred1 = y_pred
                y_test1 = y_test
            else:
                region_index = (regionList == region)
                y_score1 = y_score[region_index]
                y_pred1 = y_pred[region_index]
                y_test1 = y_test[region_index]
            ret_pred[(tool + "_pred", region)] = list(y_pred1)
            ret_pred[(tool + "_score", region)] = list(y_score1)
            ret_pred[('y_test', region)] = list(y_test1)

            ret = report_model_performance(tool, y_test1, y_pred1, y_score1, region_name=region, dsname=dsname)
            dataset.append(ret)
    perfDF = pd.DataFrame.from_dict(dataset)
    logger.info(f"perfDF={perfDF}")

    return perfDF, ret_pred


def validate_args(args):
    if len(args.train) != len(args.train_chr):
        raise Exception(f"Params length for train != train_chr: train={args.train}, train_chr={args.train_chr}")
    if len(args.test) != len(args.test_chr):
        raise Exception(f"Params length for test != test_chr: test={args.test}, test_chr={args.test_chr}")


def parse_arguments():
    parser = argparse.ArgumentParser(prog='cs_train (NANOME)', description='Consensus model train on data')
    parser.add_argument('--train', nargs='+', required=True,
                        help='train data file')
    parser.add_argument('--train-chr', nargs='+', required=True,
                        help='train chr file')
    parser.add_argument('--test', nargs='+', required=True,
                        help='test data file')
    parser.add_argument('--test-chr', nargs='+', required=True,
                        help='train chr file')
    parser.add_argument('--input-tools', nargs='+', type=str, default=['megalodon', 'nanopolish', 'deepsignal'],
                        help='input features for train, default is megalodon, nanopolish, and deepsignal')
    parser.add_argument('--dsname', type=str, default="NA12878",
                        help='dataset name, default is NA12878')
    parser.add_argument('--model-name', type=str, default="basic",
                        help='model name: basic, etc.')
    parser.add_argument('--base-model', type=str, default="xgboost",
                        help='base model name: rf, xgboost, etc.')
    parser.add_argument('-o', type=str, required=True,
                        help='output file dir')
    parser.add_argument('--niter', type=int, default=20,
                        help='number of iterations for random CV, default is 20')
    parser.add_argument('--cv', type=int, default=3,
                        help='number of CV, default is 3')
    parser.add_argument('--scoring', type=str, default='f1',
                        help='optimized score name, i.e., f1, roc_auc, etc., default is f1')
    parser.add_argument('--random-state', type=int, default=42,
                        help='random state 42')
    parser.add_argument('--processors', type=int, default=1,
                        help='number of processors, default is 1')
    parser.add_argument('--test-lines', type=int, default=None,
                        help='test top N rows, such as 10000, default is None')
    parser.add_argument('--show-confusion-matrix', help="if output verbose info", action='store_true')
    parser.add_argument('--apply-cutoff', help="if apply default cutoff of tools", action='store_true')
    parser.add_argument('--apply-cutoff-train', help="if apply default cutoff of tools before train",
                        action='store_true')
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

    validate_args(args)
    train_chrs = '_'.join(args.train_chr)
    report_mm_name = args.base_model + "_" + args.model_name
    logger.info(
        f"Train model name: {args.model_name}_{train_chrs}, base_mode={args.base_model}, report_mm_name={report_mm_name}")

    ## Get train data
    X, y, toolDF = get_data(args.train, model_name=args.model_name, input_tools=args.input_tools)

    ## apply cutoff before train
    if args.apply_cutoff_train:
        raise Exception("Under development")
        ## apply cutoff of tools, filter X_test, y_test, etc.
        logger.debug(f"cutoff_llr_tools={cutoff_llr_tools}")
        logger.debug(f"Apply cutoff to train data, before cutoff, len={len(toolDF)}")
        keep_index = None
        for tool in cutoff_llr_tools:
            if tool not in toolDF:
                continue
            logger.debug(f"Apply LLR cutoff={cutoff_llr_tools[tool]:.2f} for tool={tool}")
            if keep_index is None:
                keep_index = toolDF[tool].abs() >= cutoff_llr_tools[tool]
            else:
                keep_index = keep_index & (toolDF[tool].abs() >= cutoff_llr_tools[tool])
        X = X[keep_index].reset_index(drop=True).copy()
        y = y[keep_index].reset_index(drop=True).copy()
        toolDF = toolDF[keep_index].reset_index(drop=True).copy()
        logger.debug(f"Apply cutoff, after cutoff, len={len(toolDF)}")

    ## train model
    mm = train_model(X, y, model_name=args.model_name, region_vector=toolDF['Region'], scoring=args.scoring)

    ## Eval the trained model
    total_ret_pred = None
    outdflist = []
    for infn, chr in zip(args.test, args.test_chr):
        logger.info(f"Evaluate chr={chr}, infn={infn}")
        X_test, y_test, toolDF = get_data(infn, model_name=args.model_name, input_tools=args.input_tools)

        dsname = '_'.join([args.dsname, chr])  # NA12878_chr1
        y_pred, y_score = model_predict(mm, X_test)
        y_llr2 = y_score.apply(prob_to_llr2)
        y_llr2.rename(report_mm_name, inplace=True)
        logger.debug(f"y_pred={y_pred}, y_score={y_score}, y_llr2={y_llr2}")
        toolDF = pd.concat([toolDF, y_llr2], axis=1)
        # logger.debug(f"toolDF={toolDF}")

        regionList = toolDF['Region']
        toolDF = toolDF[top3_tools + [report_mm_name]]
        # logger.debug(f"toolDF={toolDF}, regionList={regionList}")

        if args.apply_cutoff:
            ## apply cutoff of tools, filter X_test, y_test, etc.
            logger.debug(f"cutoff_llr_tools={cutoff_llr_tools}")
            logger.debug(f"Apply cutoff, before cutoff, len={len(toolDF)}")
            keep_index = None
            for tool in cutoff_llr_tools:
                if tool not in toolDF:
                    continue
                logger.debug(f"Apply LLR cutoff={cutoff_llr_tools[tool]:.2f} for tool={tool}")
                if keep_index is None:
                    keep_index = toolDF[tool].abs() >= cutoff_llr_tools[tool]
                else:
                    keep_index = keep_index & (toolDF[tool].abs() >= cutoff_llr_tools[tool])
            y_test = y_test[keep_index].reset_index(drop=True).copy()
            toolDF = toolDF[keep_index].reset_index(drop=True).copy()
            regionList = regionList[keep_index].reset_index(drop=True).copy()
            logger.debug(f"Apply cutoff, after cutoff, len={len(toolDF)}")
            # logger.debug(f"toolDF={toolDF}, regionList={regionList}, y_test={y_test}")

        if len(y_test) > 0:
            df, ret_pred = eval_tools(y_test, toolDF=toolDF, regionList=regionList, dsname=dsname)
            outdflist.append(df)

            if total_ret_pred is None:
                total_ret_pred = ret_pred
            else:
                for key in ret_pred:
                    total_ret_pred[key].extend(ret_pred[key])
        else:
            logger.error(f"After cutoff, no data exist for chr={chr}, infn={infn}")
    outdf = pd.concat(outdflist)

    ## aggregate all test data predictions, report the performance
    dataset_all = []
    for tool in tool_list + [report_mm_name]:
        for region in region_order:
            if (f"{tool}_pred", region) not in total_ret_pred:
                continue
            y_test = pd.Series(total_ret_pred[('y_test', region)])
            y_pred = total_ret_pred[(f"{tool}_pred", region)]
            y_score = total_ret_pred[(f"{tool}_score", region)]

            if len(y_test) > 0:
                ret = report_model_performance(tool, y_test, y_pred, y_score, region_name=region,
                                               dsname="NA12878_test_all")
                dataset_all.append(ret)
            else:
                logger.error(f"No data for agg all, tool={tool}, region={region}")

    if len(dataset_all) < 1:
        logger.error("No perf results, bye bye")
        sys.exit(0)
    perf_on_all_test = pd.DataFrame.from_dict(dataset_all)
    logger.debug(perf_on_all_test)

    outdf = pd.concat([outdf, perf_on_all_test])

    ## Save the performance report
    outfn = os.path.join(args.o, f'NA12878_{train_chrs}_{report_mm_name}_performance_report_v1.csv')
    outdf.to_csv(outfn, index=False)
    logger.info(f"save to {outfn}")

    ## Report selected tools
    outdf1 = outdf[~outdf['Tool'].isin(exclude_report_tools)].copy()
    outdf1['Region'] = pd.Categorical(outdf1['Region'],
                                      categories=region_order,
                                      ordered=True)
    outdf1.sort_values(by=['Dataset', 'Region', args.scoring.upper()], ascending=[True, True, False], inplace=True)
    outfn = os.path.join(args.o, f'NA12878_{train_chrs}_{report_mm_name}_performance_report_v2.csv')
    outdf1.to_csv(outfn, index=False)
    logger.info(f"perfDF1={outdf1}")
    logger.info(f"save to {outfn}")

    logger.info(f"## train and test for report_mm_name={report_mm_name} DONE")
