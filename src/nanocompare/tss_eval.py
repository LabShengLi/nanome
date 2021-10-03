#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : tss_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

"""
Export read/site-level methylation results for TSS/METEORE analysis in nanome paper.
"""
import argparse
from multiprocessing import Manager

from nanocompare.eval_common import *
from nanocompare.global_settings import get_tool_name, save_done_file, nanome_version


def import_and_save_meteore(callfn, callencode, outfn):
    import_call(callfn, callencode, baseFormat=1, enable_cache=False, using_cache=False,
                include_score=False, siteLevel=False, save_unified_format=True, outfn=outfn,
                filterChr=args.chrs)


def output_calldict_to_unified_bed_as_0base(dictCalls, outfn, sep='\t'):
    """
    Assume dictCalls are key->value, key=(chr, 123, +), value=[(freq, cov), ...], note is 1-based format
    Output is format: 0-based start, 1-base end format with tab-sep for TSS analysis

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
            strlist = [key[0], str(key[1] - 1), str(key[1]), '.', '.', key[2], str(dictCalls[key][0]),
                       str(dictCalls[key][1])]
            outf.write(sep.join(strlist) + '\n')
    logger.debug(f'Output for TSS analysis: {outfn}')


def import_and_save_site_level(callfn, callname, callencode, minToolCovCutt, outfn, ns):
    ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache,
                          using_cache=using_cache, include_score=False, siteLevel=True, filterChr=args.chrs)

    ontCallWithCov = readLevelToSiteLevelWithCov(ontCall, minCov=minToolCovCutt, toolname=callname)
    ns.callsites[callname] = len(ontCallWithCov)
    output_calldict_to_unified_bed_as_0base(ontCallWithCov, outfn)


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='tss_eval (NANOME)',
        description='Export read/site level methylation results of all nanopore tools in nanome paper')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{nanome_version}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces',
                        required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name2>",
                        default=None)
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default='\t')
    parser.add_argument('--processors', type=int, help="running processors", default=1)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--output-unified-format', help="True:output read level results(1-based start), False: output site-level results(0-based start)", action='store_true')
    parser.add_argument('--chrs', nargs='+', help='chromosome list',
                        default=humanChrSet)
    parser.add_argument('--tagname', type=str, help="output unified file tagname", default=None)
    parser.add_argument('--debug', help="if output debug info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    if args.debug:
        set_log_debug_level()
    else:
        set_log_info_level()

    logger.debug(f"args={args}")

    if args.output_unified_format:  ## if output METEORE format, must read directly
        enable_cache = False
        using_cache = False
    else:
        enable_cache = args.enable_cache
        using_cache = args.using_cache

    # runid is always like 'MethCorr-K562_WGBS_2Reps', remove first word as RunPrefix like K562_WGBS_2Reps
    RunPrefix = args.runid.replace('TSS-', '')

    minToolCovCutt = 1
    bgtruthCutt = 1

    # Currently we only use 1-base start format, for BED of singletons, non-singletons are use 1-base format
    baseFormat = 1

    # output csv seperator: , or tab
    sep = args.sep

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    logger.debug(args)
    logger.debug(f'\n\n####################\n\n')

    if args.output_unified_format:
        ## Output read-level unified format with 1-based start
        logger.info(f"We are outputing each tool's unified results for read-level, same as METEORE format")
        input_list = []
        for callstr in args.calls:
            callencode, callfn = callstr.split(':')
            if len(callfn) == 0:
                continue
            callname = get_tool_name(callencode)
            # Consider tools have read-level outputs, except for DeepMod
            if callname not in ['Nanopolish', 'Megalodon', 'DeepSignal', 'Guppy', 'Tombo']:
                continue
            outfn = os.path.join(out_dir,
                                 f"{args.dsname}_{callname}{f'-{args.tagname}' if args.tagname else ''}-perRead-score.tsv.gz")

            input_list.append((callfn, callencode, outfn,))
        with Pool(processes=args.processors) as pool:
            pool.starmap(import_and_save_meteore, input_list)

        save_done_file(out_dir)
        logger.info("### Unified format output DONE")
        sys.exit(0)

    ## Convert into 0-based format site level bed CpG files, used for TSS plot
    if args.bgtruth:
        logger.debug("We are generating bed CpG results for BG-Truth")
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
            bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True,
                                      using_cache=using_cache, enable_cache=enable_cache)
            bgTruthList.append(bgTruth1)

        # Combine one/two replicates, using cutoff=1 or 5
        bgTruth = combineBGTruthList(bgTruthList, covCutoff=bgtruthCutt)
        outfn = os.path.join(out_dir, f'{args.dsname}_BSseq-perSite-cov{bgtruthCutt}.bed.gz')
        logger.debug(f'Combined BS-seq data (cov>={bgtruthCutt}), all methylation level sites={len(bgTruth):,}')
        output_calldict_to_unified_bed_as_0base(bgTruth, outfn)
        logger.debug(f'\n\n####################\n\n')

    logger.debug("We are outputing bed CpG results for each tool")

    mgr = Manager()
    ns = mgr.Namespace()

    # callname -> # of sites
    ns.callsites = mgr.dict()

    input_list = []
    for callstr in args.calls:
        # logger.info(f'\n\n####################\n\n')
        callencode, callfn = callstr.split(':')
        if len(callfn) == 0:
            continue
        callname = get_tool_name(callencode)

        outfn = os.path.join(out_dir, f'{args.dsname}_{callname}-perSite-cov{minToolCovCutt}.bed.gz')
        input1 = (callfn, callname, callencode, minToolCovCutt, outfn, ns,)
        input_list.append(input1)

    with Pool(processes=args.processors) as pool:
        pool.starmap(import_and_save_site_level, input_list)

    for key in ns.callsites.keys():
        logger.debug(f"tool={key}, sites={ns.callsites[key]}")
    save_done_file(out_dir)
    logger.info("TSS bed file generation DONE")
