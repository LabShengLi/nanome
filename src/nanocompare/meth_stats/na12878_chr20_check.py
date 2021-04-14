"""
Check DeepMod performance on chr 20 of NA12878,
to see if results as what they reported?


Results:

Running on Linux
2020-02-14 22:46:33,092 - [Universal_meth_stats_evaluation.py:867] - INFO: ###	importPredictions_DeepMod SUCCESS: 1931856 methylation calls mapped to 1931856 CpGs from /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/NA12878_ultra_long_DeepMod_calls/NA12878_ultra_long_DeepMod_calls.C.chr20.C.bed file
2020-02-14 22:46:35,501 - [Universal_meth_stats_evaluation.py:867] - INFO: ###	importPredictions_DeepMod SUCCESS: 1370586 methylation calls mapped to 1370586 CpGs from /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/NA12878_ultra_long_DeepMod_calls/NA12878_ultra_long_DeepMod_calls.C.chr20.C.atLeast10.bed file

2020-02-14 22:46:40,864 - [na12878_chr20_check.py:122] - INFO: minCov=1
2020-02-14 22:46:42,313 - [na12878_chr20_check.py:42] - INFO: cutoff=0.90, Unmeth=221177, Meth=314429
2020-02-14 22:46:43,607 - [na12878_chr20_check.py:97] - INFO:
cr=              precision    recall  f1-score   support
          5C       0.63      0.98      0.77    107305
         5mC       0.99      0.80      0.89    310430
    accuracy                           0.85    417735
   macro avg       0.81      0.89      0.83    417735
weighted avg       0.90      0.85      0.85    417735
2020-02-14 22:46:44,241 - [na12878_chr20_check.py:101] - INFO:
roc_auc=0.970588961631317

2020-02-14 22:46:49,068 - [na12878_chr20_check.py:122] - INFO: minCov=5
2020-02-14 22:46:49,992 - [na12878_chr20_check.py:42] - INFO: cutoff=0.90, Unmeth=154669, Meth=259434
2020-02-14 22:46:50,688 - [na12878_chr20_check.py:97] - INFO:
cr=              precision    recall  f1-score   support
          5C       0.60      0.99      0.74     75550
         5mC       1.00      0.80      0.89    256602
    accuracy                           0.85    332152
   macro avg       0.80      0.90      0.82    332152
weighted avg       0.90      0.85      0.86    332152
2020-02-14 22:46:51,100 - [na12878_chr20_check.py:101] - INFO:
roc_auc=0.9750021335190482

2020-02-14 22:46:55,318 - [na12878_chr20_check.py:122] - INFO: minCov=10
2020-02-14 22:46:55,911 - [na12878_chr20_check.py:42] - INFO: cutoff=0.90, Unmeth=90907, Meth=174100
2020-02-14 22:46:56,377 - [na12878_chr20_check.py:97] - INFO:
cr=              precision    recall  f1-score   support
          5C       0.57      0.99      0.72     44232
         5mC       1.00      0.81      0.89    172343
    accuracy                           0.84    216575
   macro avg       0.78      0.90      0.81    216575
weighted avg       0.91      0.84      0.86    216575
2020-02-14 22:46:56,647 - [na12878_chr20_check.py:101] - INFO:
roc_auc=0.9776205628818496
"""

from sklearn.metrics import classification_report, roc_auc_score

from lilab.tcga.global_tcga import *

from nanocompare.read_level_eval import importGroundTruth_from_Encode, importPredictions_DeepMod3

bgENCFF279HCL_fn = '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/NA12878/ENCFF279HCL.chr20.bed'

bgENCFF835NTC_fn = '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/NA12878/ENCFF835NTC.chr20.bed'

bgJoin_fn = '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/ENCFF279HCL_and_ENCFF835NTC_combined.chr20.DeepMod_style.bed'

deepmod_fn = '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/NA12878_ultra_long_DeepMod_calls/NA12878_ultra_long_DeepMod_calls.C.chr20.C.bed'

deepmod_atlease10_fn = '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/NA12878_ultra_long_DeepMod_calls/NA12878_ultra_long_DeepMod_calls.C.chr20.C.atLeast10.bed'


def construct_meth_unmeth_bgtruth(t1, t2, cutoff=0.9):
    """
    Construct fully-methylated and unmethylated CpGs from two replicates.
    :param t1:  gb truth 1
    :param t2:  gb truth 2
    :param cutoff: >= cutoff means fully-methylated
    :return:    fully-methylated and unmethylated CpG sites
    """
    meth = dict()
    unmeth = dict()
    notused = dict()
    for cpg in t1.keys():
        if cpg not in t2:
            continue

        t1m = t1[cpg]
        t2m = t2[cpg]
        if t1m >= cutoff and t2m >= cutoff:
            meth[cpg] = 1.0
        elif t1m == 0.0 and t2m == 0.0:
            unmeth[cpg] = 0.0
        else:
            notused[cpg] = -1.0

    logger.info(f'cutoff={cutoff:.2f}, Unmeth={len(unmeth)}, Meth={len(meth)}')
    ret = meth.copy()
    ret.update(unmeth)
    return ret


def construct_ytrue_ypred(truth, ontCall, cutoff=0.5):
    """
    Using methylation percentage >= cutoff as positive cases
    :param truth:
    :param ontCall:
    :return:
    """
    ytrue = []
    ypred = []
    for cpg in truth.keys():
        if cpg not in ontCall:
            continue

        if truth[cpg] >= 1 - 1e-6:
            ytrue.append(1)
        else:
            ytrue.append(0)

        if ontCall[cpg] >= cutoff:
            ypred.append(1)
        else:
            ypred.append(0)

    return ytrue, ypred


# TODO: if needed
def construct_ytrue_ypred_each_read(truth, ontCall):
    pass


def construct_ytrue_ypred_for_rocauc(truth, ontCall):
    ytrue = []
    ypred = []
    for cpg in truth.keys():
        if cpg not in ontCall:
            continue

        if truth[cpg] >= 1 - 1e-6:
            ytrue.append(1)
        else:
            ytrue.append(0)

        ypred.append(ontCall[cpg])

    return ytrue, ypred


def report_classification_results(truth, ontCall):
    """
    Based on truth and ontCall, report classification results
    :param truth:
    :param ontCall:
    :return:
    """
    yt, yp = construct_ytrue_ypred(truth, ontCall)
    logger.debug(f'yt={len(yt)}, yp={len(yp)}')

    target_names = ['5C', '5mC']
    cr = classification_report(yt, yp, target_names=target_names)
    logger.info(f'\ncr={cr}')

    yt, yp = construct_ytrue_ypred_for_rocauc(truth, ontCall)
    roc_auc = roc_auc_score(yt, yp)
    logger.info(f'\nroc_auc={roc_auc}')

    pass


def main():
    # [1,5,10]

    deepmodOntCall1 = importPredictions_DeepMod3(deepmod_fn)

    deepmodOntCall_atleast10 = importPredictions_DeepMod3(deepmod_atlease10_fn)

    logger.debug(f'deepmodOntCall1 = {len(deepmodOntCall1)}, deepmodOntCall_atleast10 = {len(deepmodOntCall_atleast10)}')

    for minCov in [1, 5, 10]:
        bgTruth1 = importGroundTruth_from_Encode(bgENCFF279HCL_fn, covCutt=minCov)
        logger.debug(f'bgTruth1 = {len(bgTruth1)}')

        bgTruth2 = importGroundTruth_from_Encode(bgENCFF835NTC_fn, covCutt=minCov)
        logger.debug(f'bgTruth2 = {len(bgTruth2)}')

        logger.info(f'minCov={minCov}')

        # truth = construct_meth_unmeth_gbtruth(bgTruth1, bgTruth2, cutoff=1.0)
        # report_classification_results(truth, deepmodOntCall1)

        truth = construct_meth_unmeth_bgtruth(bgTruth1, bgTruth2, cutoff=0.9)
        report_classification_results(truth, deepmodOntCall1)

    return


if __name__ == '__main__':
    set_log_info_level()
    main()
    pass
