#!/usr/bin/env python3
# @Author   : Emma Wadee
# @FileName : model_predict.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome
import argparse

import joblib
import pandas as pd

from nanome.weight_learning.utils import str2bool, prepare_one_hot_sequence, prepare_performer

CHUNKSIZE = 500000

'''
python model_predict.py \
    --importmodelfile "/pod/2/li-lab/wadee/progress/nanome/7-18/models/seq_weight_xgboost_model.joblib" \
    --dsname "chr1_testing" \
    --outdir "/pod/2/li-lab/wadee/progress/nanome/7-21/results" \
    --megalodon_file "/fastscratch/liuya/nanome/temp/ReadLevel-NA12878_NANOME/NA12878_Megalodon-perRead-score.sort.tsv.gz" \
    --deepsignal_file "/fastscratch/liuya/nanome/temp/ReadLevel-NA12878_NANOME/NA12878_DeepSignal-perRead-score.sort.tsv.gz" \
    --nanopolish_file "/fastscratch/liuya/nanome/temp/ReadLevel-NA12878_NANOME/NA12878_Nanopolish-perRead-score.sort.tsv.gz" \
    --feature_file "/fastscratch/liuya/nanome/temp/na12878_deepsignal2/NA12878_CHR1_deepsignal2_features_combine.tsv.gz" \
    --chrs "chr1" --model_type "XGBoost" --one_hot_sequence true 

'''


def parse_arguments():
    parser = argparse.ArgumentParser("")

    parser.add_argument('--importmodelfile', type=str, required=True,
                        help=f'model file, existing model list...')
    parser.add_argument('--dsname', type=str, required=True,
                        help='dataset name')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output file name')
    # parser.add_argument('-t', nargs='+', help='tools used for prediction, default is None',
    # default=None) ###change???

    parser.add_argument('--megalodon_file', type=str, required=False,
                        help="Megalodon prediction", default=None)
    parser.add_argument('--deepsignal_file', type=str, required=False,
                        help="DeepSignal prediction", default=None)
    parser.add_argument('--nanopolish_file', type=str, required=False,
                        help="Nanopolish prediction", default=None)

    parser.add_argument('--base_model_type', type=str, default="XGBoost", choices=["RF", "XGBoost"], required=False,
                        help="model type, Random Forest or XGBoost")
    parser.add_argument('--specific_model_type', type=str, default="base",
                        choices=["base", "weight", "seq_weight", "xgboost_seq_weight"], required=False,
                        help="model type, Random Forest or XGBoost")
    parser.add_argument("--one_hot_sequence", type=str2bool, nargs='?',
                        const=True, default=False, required=False,
                        help="Save model to file")  # add conditional parameter requirement? Like if true, then require output file

    parser.add_argument('--feature_file', type=str, required=True,
                        help="Feature extraction file")

    parser.add_argument('--random-state', type=int, default=42,
                        help='random state, default is 42')
    # parser.add_argument('--processors', type=int, default=8,
    # help='num of processors, default is 8')
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE,
                        help=f'chunk size for load large data, default is {CHUNKSIZE}')  ###change to chromosome???
    parser.add_argument('--nrows', type=int, default=None,
                        help=f'Number of rows for testing')  ###change to chromosome???
    # parser.add_argument('--dropna', help="if drop NA values before prediction", action='store_true')
    # parser.add_argument('--tsv-input', help="if input is tsv for tools' read-level format, or else is combined input",
    # action='store_true')
    parser.add_argument('--chrs', nargs='+', help='chromosomes used', default=None)
    # parser.add_argument('--interactive', help="if output to console as interactive mode, quit use q/Q",
    # action='store_true')
    # parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    print(args)

    '''
    3 tools' read level files, and one feature files
    feature_file = "/projects/li-lab/yang/dl_model/NA12878_data/prediction_input_human_sample/TestData_deepsignal2_feature_combine.tsv.gz"
    megalodon_file = "/projects/li-lab/yang/dl_model/NA12878_data/prediction_input_human_sample/TestData_Megalodon-perRead-score.tsv.gz"
    deepsignal_file = "/projects/li-lab/yang/dl_model/NA12878_data/prediction_input_human_sample/TestData_DeepSignal-perRead-score.tsv.gz"  
    nanopolish_file = "/projects/li-lab/yang/dl_model/NA12878_data/prediction_input_human_sample/TestData_Nanopolish-perRead-score.tsv.gz"

    '''

    toolList = ['megalodon', 'deepsignal', 'nanopolish']
    trainList = toolList  # used for adding k-mer to training feature

    ###1. input three tools' scores as DF1, DF2, DF3 --> merge with ReadID, chr, position, and strand using outer join, allowing NA missing values or some columns

    if args.megalodon_file is None and args.deepsignal_file is None and args.nanopolish_file is None:
        raise Exception("No performers included. Exiting...")

    # reads in by chunksize, renames column name
    megalodon = prepare_performer(args.megalodon_file, args.chunksize, args.nrows, args.chrs, "megalodon")
    nanopolish = prepare_performer(args.nanopolish_file, args.chunksize, args.nrows, args.chrs, "nanopolish")
    deepsignal = prepare_performer(args.deepsignal_file, args.chunksize, args.nrows, args.chrs, "deepsignal")

    ###2. save merged DF_score (e.g., the columns: ID, Chr, Pos, Strand, Nanopolish, DeepSignal, Megalodon)
    print("Merging scores...")
    df_score = pd.merge(pd.merge(megalodon, deepsignal, on=['ID', 'Chr', 'Pos', 'Strand'], how='outer'), nanopolish,
                        on=['ID', 'Chr', 'Pos', 'Strand'], how='outer')

    ###3. Input sequence feature Dataframe as DF_seq, only keep ReadID, chr, position, strand, sequence column
    print("Reading feature extraction file...")

    df_seq = None
    for chunk in pd.read_csv(args.feature_file, header=None, index_col=False, sep="\t", iterator=True,
                             chunksize=args.chunksize, nrows=args.nrows):
        chunk.columns = ['Chr', 'Pos', 'Strand', 'Pos_in_Strand', 'ID', 'Read_Strand', 'seq', 'signal_means',
                         'signal_stds', 'signal_lens', 'raw_signals', 'methyl_label']
        chunk.drop(['Pos_in_Strand', 'Read_Strand', 'signal_means', 'signal_stds', 'signal_lens', 'raw_signals',
                    'methyl_label'], axis=1, inplace=True)

        if args.chrs is not None: chunk = chunk[chunk['Chr'].isin(args.chrs)]  # need???

        chunk['Pos'] = chunk['Pos'] + 1  # make 1-based

        if df_seq is None:
            df_seq = chunk
        else:
            df_seq = df_seq.add(chunk)

    print(df_seq)

    print("Merging df_seq and df_score...")
    ###4. merge DF_score with DF_seq using outer join, allowing missing values. Save to DF_pred
    df_pred = pd.merge(df_seq, df_score, on=['ID', 'Chr', 'Pos', 'Strand'], how='outer')
    print(df_pred)

    ###5. Prepare one hot encoding and deal with missing values (fill with 0 for RF, no change for XGBoost)
    if args.one_hot_sequence:
        print("One hot encoding...")
        df_pred['k_mer_filled'] = df_pred['seq']
        df_pred['k_mer_filled'] = df_pred['k_mer_filled'].fillna("N" * 17)
        df_pred = prepare_one_hot_sequence(df_pred, 'k_mer_filled')
        trainList += list(range(0, 17))  # undo hardcode length of sequence!!!!

    print(df_pred.columns)

    if args.base_model_type == "RF":
        print("Filling NAs...")
        original_pred = df_pred[toolList]
        df_pred[trainList] = df_pred[trainList].fillna(0)

    ###6. Load classifier model
    # model = "/pod/2/li-lab/wadee/progress/nanome/7-18/models/seq_weight_xgboost_model.joblib"
    print("Loading model file...")
    clf = joblib.load(args.importmodelfile)

    try:
        print(f"Model info: classifier={clf}")
        print(f"best_params={clf.best_params_}")
    except:
        print(f"WARNING: print params encounter problem")

    print("Predicting...")
    ###7. Make predictions on DF_Pred
    X_pred = df_pred[trainList].to_numpy()
    y_pred = clf.predict(X_pred)
    y_score = clf.predict_proba(X_pred)[:, 1]

    ###8. Output the predictions, DF_Pred + two columns (probability, pred_class)
    # print(df_pred[["megalodon", "nanopolish", "deepsignal"]])
    print("Preparing save...")
    df_pred['Prediction'] = y_pred
    df_pred['Prob_methylation'] = y_score
    if args.base_model_type == "RF": df_pred[
        toolList] = original_pred  # revert filled in values to original NA values???
    if args.one_hot_sequence:
        df_pred.drop(list(range(0, 17)), axis=1, inplace=True)
        df_pred.drop('K-Mer-Hot', axis=1, inplace=True)
        df_pred.drop('k_mer_filled', axis=1, inplace=True)

    print("Saving...")
    # outfn = os.path.join(args.outdir, f'{args.dsname}_{args.specific_model_type}_prediction.tsv.gz')
    outfn = f'{args.outdir}/{args.dsname}_{args.specific_model_type}_prediction.tsv.gz'
    df_pred.to_csv(outfn, index=False, sep='\t', compression="gzip")  # compression won't work!!!!
