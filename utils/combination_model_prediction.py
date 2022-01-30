#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 18:03:39 2020

@author: akanksha
"""

import argparse
import warnings
from functools import reduce

import joblib
import numpy as np
import pandas as pd
import sklearn
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

parser = argparse.ArgumentParser(description='Prediction from combined models for the reads.')

parser.add_argument('--methodsfile', '-i', type=str, required=True,
                    help='TSV file containing name and path of the method output tsv file. The output tsv file from the method should be in the format [ID,Pos,Strand,Score]. Can be compressed in gz.')

parser.add_argument('--model', '-m', choices=["default", "optimized"], required=True, type=str,
                    help='which model to select from default RF or optimized RF with max_depth 3 and n_estimator 10')

parser.add_argument('--output', '-o', type=str, required=True,
                    help='Where to store the outputs')

parser.add_argument('--model-base-dir', '-b', type=str, required=True,
                    help='Where is the model dir', default='.')

options = parser.parse_args()


def mod_file(data_file_path):
    data_file = pd.read_csv(data_file_path, header=0, sep="\t")
    name = data_file_path.split("/")[-1].split(".")[0]
    data_file['Pos'] = data_file['Pos'].astype(np.int64)
    data_file.drop_duplicates(subset=['Chr', "ID", "Pos", "Strand"], inplace=True)  # add chr
    data_file.reset_index(inplace=True, drop=True)

    df_id_strand = data_file[["ID", 'Chr', 'Pos', "Strand"]]
    df_id_strand.drop_duplicates(subset=["ID", 'Chr', 'Pos', "Strand"], inplace=True)

    # mask = data_file.index[data_file.Strand == "-"].tolist() # select negative strand
    # data_file["Pos"][mask] = data_file["Pos"][mask] - 1 # merge to positive strand
    data_file.drop(["Strand"], axis=1, inplace=True)
    data_file.rename(columns={"Score": name}, inplace=True)
    data_file.reset_index(inplace=True, drop=True)

    return (data_file, df_id_strand)


def main(mp, combine_file, combine_id_strand):
    loaded_model = joblib.load(open(mp, 'rb'))
    print(f"model file: {mp}")
    print(f"loaded_model={loaded_model}", flush=True)

    X = combine_file[combine_file.columns[3:]]  # 2:
    X = sklearn.preprocessing.MinMaxScaler().fit_transform(X)

    prediction = pd.DataFrame(loaded_model.predict(X))  ##
    prediction.rename(columns={0: "Prediction"}, inplace=True)

    prediction_prob = pd.DataFrame(loaded_model.predict_proba(X))
    prediction_prob = prediction_prob[[1]]
    prediction_prob.rename(columns={1: "Prob_methylation"}, inplace=True)

    final_output = pd.concat([combine_file[combine_file.columns[:3]], prediction, prediction_prob], axis=1)
    final_output.dropna(inplace=True)
    final_output['Pos'] = final_output['Pos'].astype(np.int64)
    final_output['Prediction'] = final_output['Prediction'].astype(int)
    outdf = final_output.merge(combine_id_strand, how='left', on=['ID', 'Chr', 'Pos'])
    outdf.dropna(inplace=True)
    outdf.to_csv(options.output, header=True, index=False, sep='\t')
    print(f"save to {options.output}", flush=True)


if __name__ == '__main__':
    METEORE_Dir = options.model_base_dir
    df_file = pd.read_csv(options.methodsfile, header=None, sep='\t')
    if options.model == "default":
        fillval = "default"
    else:
        fillval = "max_depth_3_n_estimator_10"
    modelname = '_'.join(df_file[0])
    mp = f'{METEORE_Dir}/saved_models/rf_model_' + fillval + '_' + modelname + '.model'
    dfs = []
    dfs_id_strand = []
    for i in df_file[1]:
        df, df_id_strand = mod_file(i)
        dfs.append(df)
        dfs_id_strand.append(df_id_strand)
    combine_file = reduce(lambda left, right: pd.merge(left, right, how='inner', on=["ID", "Chr", "Pos"]),
                          dfs)
    combine_file.drop_duplicates(subset=["ID", "Chr", "Pos"], inplace=True)
    combine_file.replace([np.inf, -np.inf], np.nan, inplace=True)
    combine_file.dropna(inplace=True)
    combine_file.reset_index(inplace=True, drop=True)

    # print(f"combine_file={combine_file}")
    # print(f"combine_file.columns={combine_file.columns}", flush=True)

    # add strand, no need to combine strand since report read level outputs
    combine_id_strand = pd.concat(dfs_id_strand)
    combine_id_strand.drop_duplicates(inplace=True)

    if len(combine_file) > 0:
        main(mp, combine_file, combine_id_strand)
