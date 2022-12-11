#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : cs_predict.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Predict NANOME consensus results
"""

import argparse
import os.path
import sys

import joblib
import numpy as np
import pandas as pd
import sklearn
from tqdm import tqdm

from nanome.common.eval_common import load_tool_read_level_unified_as_df
from nanome.common.global_config import set_log_debug_level, set_log_info_level, logger
from nanome.common.global_settings import CHUNKSIZE, NANOME_VERSION
from nanome.xgboost.ml_common import SITES_COLUMN_LIST, READS_COLUMN_LIST, nanome_model_dict, \
    xgboost_mode_base_dir, prob_to_llr_2, top3_tools, k_mer


class _Getch:
    """Gets a single character from standard input.  Does not echo to the screen."""

    def __init__(self):
        try:
            self.impl = _GetchUnix()
        except ImportError:
            self.impl = _GetchWindows()

    def __call__(self):
        return self.impl()


class _GetchUnix:
    def __init__(self):
        pass

    def __call__(self):
        import sys, tty, termios
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch


class _GetchWindows:
    def __init__(self):
        pass

    def __call__(self):
        import msvcrt
        return msvcrt.getch()


def split_k_mer_str(instr):
    if instr is None or (isinstance(instr, float) and np.isnan(instr)):
        return pd.Series(['N'] * k_mer)
    if len(str(instr)) != k_mer:
        raise Exception(f"instr={instr}, len={len(instr)} is not correct, since k_mer={k_mer}")
    return pd.Series(list(str(instr)))


def dummy_empty_tool_df(tool):
    empty_frame = {'Chr': pd.Series(dtype='str'), "ID": pd.Series(dtype='str'),
                   "Pos": pd.Series(dtype='int'), "Strand": pd.Series(dtype='str'),
                   tool.lower(): pd.Series(dtype='float')}
    df = pd.DataFrame(empty_frame)
    return df


def dummy_empty_seq_df():
    empty_frame = {"ID": pd.Series(dtype='str'),
                   'Chr': pd.Series(dtype='str'),
                   "Pos": pd.Series(dtype='int'),
                   "Strand": pd.Series(dtype='str'),
                   "Seq": pd.Series(dtype='str')}
    df = pd.DataFrame(empty_frame)
    return df


def model_pred(datadf, model, outfn, header=True, mode='w',
               interactive=False):
    """
    Make xgboost prediction on DF
    Args:
        datadf:
        model:
        tool_list:
        outfn:
        reads_column_list:
        header:
        mode:
        interactive:

    Returns:

    """

    X1 = datadf[top3_tools].copy().astype(float)

    if args.model_specific.startswith('rf_'):
        # deal with inf value, NA value will be handeled by model's impute preprocessor
        X1.replace([np.inf, -np.inf], np.nan, inplace=True)

    if 'Seq' in datadf:
        X2 = datadf['Seq'].apply(split_k_mer_str).copy()
        X2.columns = [f"DNASeq_{k}" for k in range(len(X2.columns))]
        inputX = pd.concat([X1, X2], axis=1)
    else:
        inputX = X1

    prediction = pd.DataFrame(model_cls.predict(inputX))
    prediction.rename(columns={0: "Prediction"}, inplace=True)

    prediction_prob = pd.DataFrame(model.predict_proba(inputX))[[1]]
    prediction_prob.rename(columns={1: "Prob_methylation"}, inplace=True)

    outdf = pd.concat([datadf, prediction, prediction_prob], axis=1)
    outdf = outdf[feature_out + ["Prediction", "Prob_methylation"]]
    logger.debug(f"outdf={outdf}")

    outdf.to_csv(outfn, sep='\t', index=False, header=header, mode=mode)
    logger.info(f"make predictions:{len(outdf):,}")
    logger.info(f"save to {outfn}")

    if interactive:
        ## Get single character from input keyboard
        getch = _Getch()

        print("\n\nNow you are in interactive mode to check inference, use key q/Q to quit, any other key to continue")
        outdf['LLR'] = outdf['Prob_methylation'].apply(prob_to_llr_2)
        print('\t'.join(list(outdf.columns)))
        for index, row in outdf.iterrows():
            row_str_list = [str(k) for k in row]
            print('\t'.join(row_str_list))
            value = getch()
            if value.lower() == 'q':
                break


def parse_arguments():
    parser = argparse.ArgumentParser(prog='cs_predict (NANOME)', description='Consensus model predict for data')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', nargs='+', default=None,
                        help='input tsv combined data for predicting')
    parser.add_argument('--nanopolish', type=str, default=None,
                        help='input nanopolish unified read-level file')
    parser.add_argument('--megalodon', type=str, default=None,
                        help='input megalodon unified read-level file')
    parser.add_argument('--deepsignal', type=str, default=None,
                        help='input deepsignal unified read-level file')
    parser.add_argument('--feature', type=str, default=None,
                        help='input feature file for DNAseq')
    parser.add_argument('--feature-readids-col', nargs='+', default=[0, 1, 2, 4],
                        help='column index for ID, Chr, Pos and Strand')
    parser.add_argument('--feature-readids-col-order', nargs='+', default=[3, 0, 1, 2],
                        help='column index order for ID, Chr, Pos and Strand')
    parser.add_argument('--feature-seq-col', type=int, default=6,
                        help='column index for DNA seq feature')
    parser.add_argument('--model_specific', type=str, default='xgboost_basic_seq_w',
                        help='specific model info')
    parser.add_argument('-m', type=str, required=True,
                        help=f'model file, existing model list: {",".join(list(nanome_model_dict.keys()))}')
    parser.add_argument('--dsname', type=str, required=True,
                        help='dataset name')
    parser.add_argument('-o', type=str, required=True,
                        help='output file name')
    parser.add_argument('-t', nargs='+', help='tools used for prediction, default is None',
                        default=None)
    parser.add_argument('--random-state', type=int, default=42,
                        help='random state, default is 42')
    parser.add_argument('--processors', type=int, default=8,
                        help='num of processors, default is 8')
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE,
                        help=f'chunk size for load large data, default is {CHUNKSIZE}')
    parser.add_argument('--inner-join', help="if inner join for merge data, default is outer join", action='store_true')
    parser.add_argument('--chrs', nargs='+', help='chromosomes used', default=None)
    parser.add_argument('--interactive', help="if output to console as interactive mode, quit use q/Q",
                        action='store_true')
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

    sklearn.show_versions()

    if args.inner_join:
        how_merge = 'inner'
    else:
        how_merge = 'outer'

    if args.m in nanome_model_dict:
        infn = os.path.join(xgboost_mode_base_dir, nanome_model_dict[args.m])
    else:
        infn = args.m
    logger.info(f"XGBoost model file: {infn}")
    model_cls = joblib.load(infn)
    try:
        logger.debug(f"Model info: xgboost_cls={model_cls}")
        logger.debug(f"best_params={model_cls.best_params_}")
    except:
        logger.debug(f"WARNING: print params encounter problem")

    feature_out = list(READS_COLUMN_LIST)
    chr_col_in_feature = args.feature_readids_col_order[1]

    if 'basic' in args.model_specific.lower():
        feature_out += top3_tools

    ## load tools results
    logger.info(f"Load tools results")
    if args.i is None:
        # input 3 tools files
        if args.nanopolish is None and args.megalodon is None and args.deepsignal is None:
            raise Exception("All tools input is None, not support.")

        toolDict = {
            'megalodon': args.megalodon,
            'nanopolish': args.nanopolish,
            'deepsignal': args.deepsignal,
        }

        datadf = None
        for tool, infn in toolDict.items():
            if infn is not None:
                df1 = load_tool_read_level_unified_as_df(
                    infn, toolname=tool,
                    filterChrSet=set(args.chrs) if args.chrs is not None else None,
                    chunksize=args.chunksize)
            else:
                df1 = dummy_empty_tool_df(tool)
            if datadf is None:
                datadf = df1
            else:
                datadf = datadf.merge(df1, on=READS_COLUMN_LIST, how=how_merge)
    else:
        # combined joined preds db tsv as input
        dflist = []
        for infn in tqdm(args.i):
            if args.chrs is not None and len(args.chrs) >= 1:
                iter_df = pd.read_csv(infn, header=0, index_col=False, sep=",", iterator=True,
                                      chunksize=args.chunksize)
                if args.chrs is not None:
                    datadf1 = pd.concat([chunk[chunk['Chr'].isin(set(args.chrs))] for chunk in iter_df])
                else:
                    datadf1 = pd.concat([chunk for chunk in iter_df])
            else:
                datadf1 = pd.read_csv(infn, header=0, sep=',', index_col=False)

            datadf1 = datadf1[READS_COLUMN_LIST + top3_tools]
            datadf1.dropna(subset=top3_tools, inplace=True, how='all')
            datadf1.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)
            dflist.append(datadf1)
        datadf = pd.concat(dflist)
        datadf = datadf[READS_COLUMN_LIST + top3_tools]
        datadf.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)

    logger.debug(f"datadf={datadf}")

    ## load features (DNAseq)
    if 'basic_w_seq' in args.model_specific.lower():
        # include DNA seq feature
        logger.info(f"Load DNA seq feature")
        feature_out += ['Seq']
        if args.feature is None:
            seqDF = dummy_empty_seq_df()
        else:
            col_read = args.feature_readids_col + [args.feature_seq_col]
            col_order = args.feature_readids_col_order + [len(args.feature_readids_col_order)]
            logger.debug(f"col_read={col_read}, col_order={col_order}")
            iter_df = pd.read_csv(args.feature, header=None, index_col=False, usecols=col_read, sep="\t",
                                  iterator=True,
                                  chunksize=args.chunksize)
            if args.chrs is not None:
                seqDF = pd.concat(
                    [chunk[chunk.iloc[:, chr_col_in_feature].isin(args.chrs)].iloc[:, col_order] for chunk in iter_df])
            else:
                seqDF = pd.concat([chunk.iloc[:, col_order] for chunk in iter_df])
            seqDF.columns = ['ID', 'Chr', 'Pos', 'Strand', 'Seq']
        # note the extracted feature file is 0-based, and read-level unified is 1-based
        # align the feature file into read-level unified file
        seqDF['Pos'] = seqDF['Pos'] + 1
        datadf = datadf.merge(seqDF, on=READS_COLUMN_LIST, how='left')

    datadf.drop_duplicates(subset=READS_COLUMN_LIST, inplace=True)
    datadf.reset_index(drop=True, inplace=True)
    logger.debug(datadf)

    if len(datadf) < 1:
        ## Build empty outputs with header only
        key_list = feature_out + ['Prediction', 'Prob_methylation']
        empty_frame = {keystr: [] for keystr in key_list}
        outdf = pd.DataFrame(empty_frame)
        outdf.to_csv(args.o, sep='\t', index=False)
        logger.info(f"make no predictions, bye bye!!!")
        sys.exit(0)

    ## Output read-level, site-level distributions
    logger.debug(f"Read stats: total={len(datadf):,}")

    sitedf = datadf[SITES_COLUMN_LIST].drop_duplicates()
    logger.debug(f"Site stats: total={len(sitedf):,}")
    sitedf = None

    logger.debug(f"Start predict by {args.model_specific}......")
    model_pred(datadf, model_cls, args.o, interactive=args.interactive)
    logger.info(f"### Done for model {args.model_specific} from {args.m} predict")
