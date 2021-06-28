#!/usr/bin/env python3

"""
Generate site-level methylation results for TSS analysis in Nanocompare paper.
"""
import argparse

from nanocompare.eval_common import *
from nanocompare.global_settings import get_tool_name


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Site level correlation analysis')
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces', required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name1>", required=True)
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--beddir', type=str, help="base dir for bed files", default=None)  # need perform performance evaluation before, then get concordant, etc. bed files, like '/projects/li-lab/yang/results/2021-04-01'
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=' ')
    parser.add_argument('--processors', type=int, help="running processors", default=8)
    parser.add_argument('--min-bgtruth-cov', type=int, help="cutoff of coverage in bg-truth", default=1)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff of coverage in nanopore calls", default=1)
    parser.add_argument('--baseFormat', type=int, help="base format after imported", default=1)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', action='store_true')
    parser.add_argument('--using-cache', action='store_true')
    parser.add_argument('--output-meteore', action='store_true')
    return parser.parse_args()


def output_dict_to_bed(dictCalls, outfn, sep='\t'):
    """
    Assume dictCalls are key->value, key=(chr, 123, +), value=[(freq, cov), ...], note is 1-based format
    Output is format: 0-based format with tab-sep for analysis

    ==============================
    chr start end . . freq  cov

    :param dictCalls:
    :param outfn:
    :return:
    """
    ## Change to gzip output
    # with open(outfn, 'w') as outf:
    with gzip.open(outfn, 'wt') as outf:
        for key in dictCalls:
            strlist = [key[0], str(key[1] - 1), str(key[1]), '.', '.', key[2], str(dictCalls[key][0]), str(dictCalls[key][1])]
            outf.write(sep.join(strlist) + '\n')
    logger.debug(f'Output for TSS analysis: {outfn}')


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()

    if args.output_meteore:
        enable_cache = False
        using_cache = False
        output_meteore = True
    else:
        # cache function same with read level
        enable_cache = args.enable_cache
        using_cache = args.using_cache
        output_meteore = False

    # runid is always like 'MethCorr-K562_WGBS_2Reps', remove first word as RunPrefix like K562_WGBS_2Reps
    RunPrefix = args.runid.replace('TSS-', '')

    ## Now, we are exporting >=1 cov results, the cov = 3,5, etc. can be later used before plotting
    # tool coverage cutoff 1, or 3, 5
    minToolCovCutt = args.toolcov_cutoff
    minToolCovCutt = 1

    # bgtruth coverage cutoff 1, or 5, 10  --min-bgtruth-cov
    bgtruthCutt = args.min_bgtruth_cov
    bgtruthCutt = 1

    # load into program format 0-base or 1-base
    baseFormat = args.baseFormat
    # Currently we only use 1-base start format, for BED of singletons, non-singletons are use 1-base format
    baseFormat = 1

    # output csv seperator: , or tab
    sep = args.sep

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    # Add logging files also to result output dir
    # add_logging_file(os.path.join(out_dir, 'run-results.log'))

    logger.debug(args)

    logger.info(f'\n\n####################\n\n')

    if not output_meteore:
        # we import multiple (1 or 2) replicates and join them
        encode, fnlist = args.bgtruth.split(':')
        fnlist = fnlist.split(';')

        if len(fnlist) > 2:
            raise Exception(f'Currently only support bgtruth with upto two, but found more: {fnlist}')

        logger.debug(f'BGTruth fnlist={fnlist}, encode={encode}')

        bgTruthList = []
        for fn in fnlist:
            if len(fn) == 0:  # incase of input like 'bismark:/a/b/c;'
                continue
            # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
            bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True, using_cache=using_cache, enable_cache=enable_cache)
            bgTruthList.append(bgTruth1)

        # Combine one/two replicates, using cutoff=1 or 5
        bgTruth = combineBGTruthList(bgTruthList, covCutoff=bgtruthCutt)

        outfn = os.path.join(out_dir, f'{RunPrefix}.tss.bgtruth.cov{bgtruthCutt}.bed.gz')

        logger.info(f'Combined BS-seq data (cov>={bgtruthCutt}), all methylation level sites={len(bgTruth):,}')

        output_dict_to_bed(bgTruth, outfn)

    logger.info(f'\n\n####################\n\n')

    callfn_dict = defaultdict()  # callname -> filename

    # callname -> [call0, call1], call0 is no-filter results, call1 is filter by cutoff, and convert to [meth-freq, meth-cov] results.
    callresult_dict = defaultdict()
    loaded_callname_list = []

    for callstr in args.calls:
        callencode, callfn = callstr.split(':')

        if len(callfn) == 0:
            continue

        callname = get_tool_name(callencode)
        callfn_dict[callname] = callfn

        # We do now allow import DeepMod.Cluster for read level evaluation
        if callencode == 'DeepMod.C':
            raise Exception(f'{callencode} is not allowed for site level evaluation, please use DeepMod.Cluster file here')

        loaded_callname_list.append(callname)

        if output_meteore:
            outfn = os.path.join(out_dir, f"{args.dsname}_{callname}-METEORE-perRead-score.tsv.gz")
        else:
            outfn = None

        ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache, include_score=False, siteLevel=True, saveMeteore=output_meteore, outfn=outfn)

        if not output_meteore:
            ontCallWithCov = readLevelToSiteLevelWithCov(ontCall, minCov=minToolCovCutt, toolname=callname)
            callresult_dict[callname] = ontCallWithCov
            outfn = os.path.join(out_dir, f'{RunPrefix}.tss.{callname}.cov{minToolCovCutt}.bed')
            output_dict_to_bed(ontCallWithCov, outfn)

    logger.info(f'\n\n####################\n\n')
    logger.info("TSS DONE")
