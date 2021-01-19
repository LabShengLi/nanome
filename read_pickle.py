import pickle

import pandas as pd
from sklearn import metrics


def read_pickle_file(file):
    with open(file, 'rb') as f:
        ret = pickle.load(f)

    ytrue = ret['DeepSignal_true']
    ypred = ret['DeepSignal_pred']

    fpr, tpr, threshold = metrics.roc_curve(ytrue, ypred)
    df = pd.DataFrame(dict(fpr=fpr, tpr=tpr))
    return df
