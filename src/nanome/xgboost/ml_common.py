import os
from pathlib import Path

import joblib
import math

from nanome.common.global_config import logger
from nanome.common.global_settings import EPSLONG

# consensus model on top of top3 tools
top3_tools = ['megalodon', 'nanopolish', 'deepsignal']

region_order = ['Genome-wide', 'Discordant', 'Concordant', 'Singleton', 'Nonsingleton']

# use k_mer for DNAseq features
k_mer = 17

TRUTH_LABEL_COLUMN = 'Truth_label'
SITES_COLUMN_LIST = ["Chr", "Pos", "Strand"]
## ID  Chr Pos Strand
READS_COLUMN_LIST = ['ID'] + SITES_COLUMN_LIST

# each type of tools' input features
tool_feature_dict = {'basic': top3_tools,
                     'megalodon_deepsignal': ['megalodon', 'deepsignal']}

# XGBoost model default dir
xgboost_mode_base_dir = os.path.join(Path(__file__).parent, 'trained_model')

# XGBoost model name shortcuts
nanome_model_dict = {
    "NANOME2T": 'NANOME_NA12878_train1.0_megalodon_deepsignal_XGBoostNA2T_model.pkl',
    "NANOME3T": 'NANOME_NA12878_train1.0_nanopolish_megalodon_deepsignal_XGBoostNA3T_niter10_model.pkl',
    "xgboost_basic": 'NA12878_chr1_xgboost_basic_model.pkl',
    "xgboost_basic_w": 'NA12878_chr1_xgboost_basic_w_model.pkl',
    "xgboost_basic_w_seq": 'NA12878_chr1_xgboost_basic_w_seq_model.pkl',
}

# XGBoost model tools' order for input
nanome_model_tool_list_dict = {
    "NANOME2T": ['megalodon', 'deepsignal'],
    "NANOME3T": ['nanopolish', 'megalodon', 'deepsignal'],
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

## new suggested
gridcv_xgboost_params = {
    'learning_rate': [0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
    'n_estimators': [5, 10, 20, 40, 60, 80, 100],
    'max_depth': [3, 6, 9],
    'subsample': [0.6, 0.8, 1],
    'colsample_bytree': [0.6, 0.8, 1],
    'reg_alpha': [0.5, 0, 2, 5],
    'reg_lambda': [0.5, 1, 2],
}

## used currently
gridcv_xgboost_params = {
    'learning_rate': [0.01, 0.05, 0.1, 0.2],
    'n_estimators': [5, 10, 20, 40, 60],
    'max_depth': [3, 6, 9],
    'subsample': [0.6, 0.8, 1],
    'colsample_bytree': [0.6, 0.8, 1],
    'reg_alpha': [0.5, 0, 2, 5],
    'reg_lambda': [0.5, 1, 2],
}

default_rf_params = {
    'n_jobs': -1
}

gridcv_rf_params = {
    'n_estimators': [5, 10, 20, 40, 80, 100],
    'max_depth': [2, 3, 6, 9, 11],
    'oob_score': [True, False],
    'min_samples_leaf': [1, 5, 10, 50, 100]
}


def prob_to_llr_2(meth_prob):
    """
    METEORE, NANOME manner
    Args:
        meth_prob:

    Returns:

    """
    return math.log2((meth_prob + EPSLONG) / (1 - meth_prob + EPSLONG))


def prob_to_llr_e(meth_prob):
    """
    Megalodon manner:
    prob_to_llr_2(0.8)
    1.999945900626566
    prob_to_llr_e(0.8)
    1.3862568622917248

    Args:
        meth_prob:

    Returns:

    """
    return math.log((meth_prob + EPSLONG) / (1 - meth_prob + EPSLONG))


def load_meteore_model(infn):
    return joblib.load(open(infn, 'rb'))


def load_xgboost_model(infn):
    logger.debug(f"Model file: {infn}")
    xgboost_cls = joblib.load(infn)
    try:
        logger.debug(f"Model info: xgboost_cls={xgboost_cls}")
        logger.debug(f"best_params={xgboost_cls.best_params_}")
    except:
        logger.debug(f"WARNNING: print params encounter problem, may due to scikit-learn version conflicts")
    return xgboost_cls
