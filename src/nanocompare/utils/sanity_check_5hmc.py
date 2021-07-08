"""
Tool for pre-processing results

"""
import gzip

from scipy.stats import pearsonr, wilcoxon

from nanocompare.eval_common import import_bgtruth, import_call, readLevelToSiteLevelWithCov
from nanocompare.global_config import *


def output_cpg_set_0base(cpg_set, outfn):
    outf = gzip.open(outfn, 'wt')
    for key in cpg_set:
        # key is (chr, start, strand)
        if isinstance(bgtruthDict[key], list) or isinstance(bgtruthDict[key], tuple):
            cov = bgtruthDict[key][1]
        else:
            cov = None
        if cov is None:
            outf.write(f"{key[0]}\t{key[1] - 1}\t{key[1]}\t.\t.\t{key[2]}\n")
        else:
            outf.write(f"{key[0]}\t{key[1] - 1}\t{key[1]}\t.\t.\t{key[2]}\t{cov}\n")
    outf.close()
    logger.info(f"save to {outfn}")


def correlation_calls(call1, call2, keySet):
    v1 = []
    v2 = []
    for key in keySet:
        if isinstance(call1[key], list) or isinstance(call1[key], tuple):
            vv1 = call1[key][0]
        else:
            vv1 = call1[key]

        if isinstance(call2[key], list) or isinstance(call2[key], tuple):
            vv2 = call2[key][0]
        else:
            vv2 = call2[key]
        v1.append(vv1)
        v2.append(vv2)
    coe, pvalue = pearsonr(v1, v2)
    return coe, pvalue


def wilcoxon_test_calls(call1, call2, keySet):
    v1 = []
    v2 = []
    for key in keySet:
        if isinstance(call1[key], list) or isinstance(call1[key], tuple):
            vv1 = call1[key][0]
        else:
            vv1 = call1[key]

        if isinstance(call2[key], list) or isinstance(call2[key], tuple):
            vv2 = call2[key][0]
        else:
            vv2 = call2[key]
        v1.append(vv1)
        v2.append(vv2)
    stats, pvalue = wilcoxon(v1, v2)
    return stats, pvalue
    pass


if __name__ == '__main__':
    set_log_debug_level()

    mlmlFilename = "/pod/2/li-lab/Nanopore_compare/data/APL_5hmC_BSseq/APL.mlml.addstrand.bscov.oxbscov.selected.bed.gz"
    mlmlDict = import_bgtruth(mlmlFilename, encode="5hmc_ziwei", covCutoff=1, includeCov=True)
    logger.info(f"mlmlDict={len(mlmlDict):,}")

    bgtruthFilename = "/projects/li-lab/Nanopore_compare/data/APL/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.convert.add.strand.tsv.gz"
    # bgtruthDict is key -> (meth_freq, cov)
    bgtruthDict = import_bgtruth(bgtruthFilename, encode="bismark", covCutoff=1, includeCov=True)
    logger.info(f"bgtruthDict={len(bgtruthDict):,}")

    newBgtruthDict = {}
    joinKeys = set(bgtruthDict.keys()).intersection(set(mlmlDict.keys()))
    for key in joinKeys:
        newBgtruthDict[key] = bgtruthDict[key]
    logger.info(f"newBgtruthDict={len(newBgtruthDict):,}")

    # megalodonDict: key-> (meth_freq, cov)
    megalodonFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.megalodon.per_read.sorted.bed.gz"
    megalodonDict = import_call(megalodonFilename, encode="Megalodon.ZW")
    megalodonDict = readLevelToSiteLevelWithCov(megalodonDict, minCov=3, toolname="Megalodon")

    # nanopolishDict: key-> (meth_freq, cov)
    nanopolishFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.nanopolish.methylation_calls.combine.tsv.gz"
    nanopolishDict = import_call(nanopolishFilename, encode="Nanopolish")
    nanopolishDict = readLevelToSiteLevelWithCov(nanopolishDict, minCov=3, toolname="Nanopolish")

    logger.info(
        f"newBgtruthDict={type(newBgtruthDict)}, megalodonDict={type(megalodonDict)}, nanopolishDict={type(nanopolishDict)}")
    logger.info(
        f"newBgtruthDict={len(newBgtruthDict)}, megalodonDict={len(megalodonDict)}, nanopolishDict={len(nanopolishDict)}")

    joinSet = set(megalodonDict.keys()).intersection(set(nanopolishDict.keys())).intersection(set(newBgtruthDict.keys()))
    coe, pvalue = correlation_calls(megalodonDict, nanopolishDict, joinSet)
    logger.info(f"Megalodon vs. Nanopolish COE={coe:.3f}, pvalue={pvalue:e}, CpGs={len(joinSet):,}")

    for callDict, callName in zip([megalodonDict, nanopolishDict], ['Megalodon', 'Nanopolish']):
        logger.debug(f"Analyze tool={callName}")
        joinSet = set(newBgtruthDict.keys()).intersection(set(callDict.keys()))

        coe, pvalue = correlation_calls(newBgtruthDict, callDict, joinSet)
        logger.info(f"BS-seq vs. {callName} COE={coe:.3f}, pvalue={pvalue:e}, CpGs={len(joinSet):,}")

        ret_diff_get_40 = set()
        ret_diff_get_20 = set()
        ret_cons_let_05 = set()
        ret_cons_let_10 = set()
        for key in joinSet:
            if isinstance(newBgtruthDict[key], list) or isinstance(newBgtruthDict[key], tuple):
                bg_methfreq = newBgtruthDict[key][0]
            else:
                bg_methfreq = newBgtruthDict[key]
            tool_methfreq = callDict[key][0]
            if abs(tool_methfreq - bg_methfreq) >= 0.4:
                ret_diff_get_40.add(key)
            if abs(tool_methfreq - bg_methfreq) >= 0.2:
                ret_diff_get_20.add(key)
            if abs(tool_methfreq - bg_methfreq) <= 0.05:
                ret_cons_let_05.add(key)
            if abs(tool_methfreq - bg_methfreq) <= 0.10:
                ret_cons_let_10.add(key)

        logger.info(
            f"tool={callName}: ret_diff_get_40={len(ret_diff_get_40)}, ret_diff_get_20={len(ret_diff_get_20)}, ret_cons_let_05={len(ret_cons_let_05)}, ret_cons_let_10={len(ret_cons_let_10)}")

        stats_let_05, pv_let05 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_cons_let_05)
        stats_let_10, pv_let10 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_cons_let_10)
        stats_get_20, pv_get20 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_diff_get_20)
        stats_get_40, pv_get40 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_diff_get_40)

        logger.info(
            f"tool={callName} wilcoxon test:\n <5% = ({stats_let_05}, {pv_let05});\n <10% = ({stats_let_10},{pv_let10});\n >20% = ({stats_get_20}, {pv_get20});\n >40% = ({stats_get_40}, {pv_get40})\n")

        outfn = f"APL.bgtruth.{callName}.regions.diff.get.40.bed.gz"
        output_cpg_set_0base(ret_diff_get_40, os.path.join(pic_base_dir, outfn))

        outfn = f"APL.bgtruth.{callName}.regions.diff.get.20.bed.gz"
        output_cpg_set_0base(ret_diff_get_20, os.path.join(pic_base_dir, outfn))

        outfn = f"APL.bgtruth.{callName}.regions.cons.let.05.bed.gz"
        output_cpg_set_0base(ret_cons_let_05, os.path.join(pic_base_dir, outfn))

        outfn = f"APL.bgtruth.{callName}.regions.cons.let.10.bed.gz"
        output_cpg_set_0base(ret_cons_let_10, os.path.join(pic_base_dir, outfn))

    logger.info("sanity_check 5hmc DONE")
