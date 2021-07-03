"""
Tool for pre-processing results

"""
import gzip

from nanocompare.eval_common import import_bgtruth, import_call, readLevelToSiteLevelWithCov
from nanocompare.global_config import *


def output_cpg_set_0base(cpg_set, outfn):
    outf = gzip.open(outfn, 'wt')
    for key in cpg_set:
        # key is (chr, start, strand)
        outf.write(f"{key[0]}\t{key[1] - 1}\t{key[1]}\t.\t.\t{key[2]}\n")
    outf.close()
    logger.info(f"save to {outfn}")


if __name__ == '__main__':
    set_log_debug_level()

    bgtruthFilename = "/pod/2/li-lab/Nanopore_compare/data/APL_5hmC_BSseq/APL.cov5.mlml.addstrand.selected.bed.gz"
    # bgtruthDict is key -> meth_freq
    bgtruthDict = import_bgtruth(bgtruthFilename, encode="5hmc_ziwei", covCutoff=None, includeCov=False)

    # megalodonDict: key-> (meth_freq, cov)
    megalodonFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.megalodon.per_read.sorted.bed.gz"
    megalodonDict = import_call(megalodonFilename, encode="Megalodon.ZW")
    megalodonDict = readLevelToSiteLevelWithCov(megalodonDict, minCov=3, toolname="Megalodon")

    # nanopolishDict: key-> (meth_freq, cov)
    nanopolishFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.nanopolish.methylation_calls.combine.tsv.gz"
    nanopolishDict = import_call(nanopolishFilename, encode="Nanopolish")
    nanopolishDict = readLevelToSiteLevelWithCov(nanopolishDict, minCov=3, toolname="Nanopolish")

    logger.info(f"bgtruthDict={type(bgtruthDict)}, megalodonDict={type(megalodonDict)}, nanopolishDict={type(nanopolishDict)}")
    logger.info(f"bgtruthDict={len(bgtruthDict)}, megalodonDict={len(megalodonDict)}, nanopolishDict={len(nanopolishDict)}")
    for callDict, callName in zip([megalodonDict, nanopolishDict], ['Megalodon', 'Nanopolish']):
        logger.debug(f"Analyze tool={callName}")
        joinSet = set(bgtruthDict.keys()).intersection(set(callDict.keys()))

        ret_diff_get_40 = set()
        ret_diff_get_20 = set()
        ret_cons_let_05 = set()
        ret_cons_let_10 = set()
        for key in joinSet:
            bg_methfreq = bgtruthDict[key]
            tool_methfreq = callDict[key][0]
            if abs(tool_methfreq - bg_methfreq) >= 0.4:
                ret_diff_get_40.add(key)
            if abs(tool_methfreq - bg_methfreq) >= 0.2:
                ret_diff_get_20.add(key)
            if abs(tool_methfreq - bg_methfreq) <= 0.05:
                ret_cons_let_05.add(key)
            if abs(tool_methfreq - bg_methfreq) <= 0.10:
                ret_cons_let_10.add(key)

        logger.info(f"tool={callName}: ret_diff_get_40={len(ret_diff_get_40)}, ret_diff_get_20={len(ret_diff_get_20)}, ret_cons_let_05={len(ret_cons_let_05)}, ret_cons_let_10={len(ret_cons_let_10)}")
        outfn = f"APL.bgtruth.{callName}.regions.diff.get.40.bed.gz"
        output_cpg_set_0base(ret_diff_get_40, os.path.join(pic_base_dir, outfn))

        outfn = f"APL.bgtruth.{callName}.regions.diff.get.20.bed.gz"
        output_cpg_set_0base(ret_diff_get_20, os.path.join(pic_base_dir, outfn))

        outfn = f"APL.bgtruth.{callName}.regions.cons.let.05.bed.gz"
        output_cpg_set_0base(ret_cons_let_05, os.path.join(pic_base_dir, outfn))

        outfn = f"APL.bgtruth.{callName}.regions.cons.let.10.bed.gz"
        output_cpg_set_0base(ret_cons_let_10, os.path.join(pic_base_dir, outfn))

        pass

    logger.info("sanity_check 5hmc DONE")
