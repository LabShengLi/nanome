import os

import joblib

from nanocompare.global_config import logger

TRUTH_LABEL_COLUMN = 'Truth_label'
SITES_COLUMN_LIST = ["Chr", "Pos", "Strand"]
READS_COLUMN_LIST = ['ID'] + SITES_COLUMN_LIST

meteore_dir = '/projects/li-lab/yang/tools/METEORE'
meteore_deepsignal_megalodon_model = os.path.join(meteore_dir, 'saved_models',
                                                  'rf_model_max_depth_3_n_estimator_10_deepsignal_megalodon.model')
meteore_deepsignal_megalodon_tool = ['deepsignal', 'megalodon']

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
