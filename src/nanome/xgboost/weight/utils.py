import argparse
import joblib
import pandas as pd
import numpy as np
import logging as logger

base2onehot_dna = {'A' : '1000',
                   'C' : '0100',
                   'G' : '0010',
                   'T' : '0001',
                   'N' : '0000',
                   'a': '1000',
                   'c': '0100',
                   't': '0001',
                   'g': '0100',
                   'n': '0000'
                   }


def one_hot_encode(sequence):
    """
    https://ziweipan.notion.site/Transfer-learning-for-methylation-classifier-39366960aaa94cafa3cb8cfc82bf265d
    Function to convert DNA sequence with One-hot encoder method.
    param sequence: The DNA sequence.
    """
    vector = []
    for i, base in enumerate(str(sequence)):
        vector.append(base2onehot_dna[base])
    return vector

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

'''
from sklearn.preprocessing import OneHotEncoder

def one_hot_encode(sequence):
    enc=OneHotEncoder(max_categories=5)
    #split_seq = np.array(list(sequence)).reshape(1,-1)
    return enc.fit_transform(sequence)
'''

def load_xgboost_model(infn):
    logger.debug(f"Model file: {infn}")
    xgboost_cls = joblib.load(infn)
    try:
        logger.debug(f"Model info: xgboost_cls={xgboost_cls}")
        logger.debug(f"best_params={xgboost_cls.best_params_}")
    except:
        logger.debug(f"WARNNING: print params encounter problem, may due to scikit-learn version conflicts")
    return xgboost_cls

def prepare_for_RF(df):
    pass

#might not want to hardcode db[['k-mer']]
def prepare_one_hot_sequence(db, column):
    db[['K-Mer-Hot']] = db[[column]].applymap(one_hot_encode) 
    db = db.join(pd.DataFrame(db['K-Mer-Hot'].to_list()))
    return db

def prepare_performer(file, chunksize, nrows = None, chrs = None, name = 'Score'):
    print(f"Reading {name}...")

    if file is None: 
        df = pd.DataFrame({'Chr': pd.Series(dtype='str'), "ID": pd.Series(dtype='str'),
                               "Pos": pd.Series(dtype='int'), "Strand": pd.Series(dtype='str'),
                               name: pd.Series(dtype='float')})
    else:
        df = pd.read_csv(file, sep='\t', header=0, iterator = True, chunksize = chunksize, nrows = nrows)
        print(df)
        if chrs is not None:
                #df_seq = pd.concat([chunk[chunk['Chr'].isin(set(args.chrs))] for chunk in iter_df]) #don't need because feature_file is split into chromosomes??
                ###Change here??####
                #pd.concat([chunk[chunk['Chr'].isin(["chrY"])] for chunk in megalodon])
                df = pd.concat([chunk[chunk['Chr'].isin(chrs)] for chunk in df])
                

        else: 
            df = pd.concat([chunk for chunk in df])

        df = df.rename(columns = {'Score': name})
        df.replace([np.inf, -np.inf], 0, inplace=True)
        #print(np.isinf(df[name]).values.sum())

    print(df)
    return df

