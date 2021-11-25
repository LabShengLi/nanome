import joblib

from nanocompare.global_config import logger


def load_xgboost_model(infn):
    logger.debug(f"Model file: {infn}")
    xgboost_cls = joblib.load(infn)
    try:
        logger.debug(f"Model info: xgboost_cls={xgboost_cls}")
        logger.debug(f"best_params={xgboost_cls.best_params_}")
    except:
        logger.debug(f"WARNNING: print params encounter problem, may due to scikit-learn version conflicts")
    return xgboost_cls
