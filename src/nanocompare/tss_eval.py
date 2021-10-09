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

from nanocompare.eval_common import *
from nanocompare.global_settings import get_tool_name, save_done_file, nanome_version


def import_and_save_read_level(callfn, callencode, outfn):
    """
    Export to unified read-level output format for each tool
    Args:
        callfn:
        callencode:
        outfn:

    Returns:

    """
    import_call(callfn, callencode, baseFormat=1, filterChr=args.chrSet,
                include_score=False, siteLevel=False, save_unified_format=True, outfn=outfn,
                enable_cache=False, using_cache=False)


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
    with gzip.open(outfn, 'wt') as outf:
        for key in dictCalls:
            strlist = [key[0], str(key[1] - 1), str(key[1]), '.', '.', key[2], str(dictCalls[key][0]),
                       str(dictCalls[key][1])]
            outf.write(sep.join(strlist) + '\n')
    logger.debug(f'Output for TSS analysis: {outfn}')


def import_and_save_site_level(callfn, callname, callencode, minToolCovCutt, outfn):
    """
    Output to unified site-level format for each tool
    Args:
        callfn:
        callname:
        callencode:
        minToolCovCutt:
        outfn:

    Returns:

    """
    ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache,
                          using_cache=using_cache, include_score=False, siteLevel=True, filterChr=args.chrSet,
                          cache_dir=ds_cache_dir)

    ontCallWithCov = readLevelToSiteLevelWithCov(ontCall, minCov=minToolCovCutt, toolname=callname)
    ontcall_tools_dict[callname] = len(ontCallWithCov)
    output_calldict_to_unified_bed_as_0base(ontCallWithCov, outfn)


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='tss_eval (NANOME)',
                                     description='Export read/site level methylation results of all nanopore tools in nanome paper')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{nanome_version}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix/output dir name", required=True)
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces',
                        required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name2>",
                        default=None)
    parser.add_argument('--read-level-format',
                        help="if true, it will output read level results (1-based start), else it will output site-level results (0-based start, 1-based end)",
                        action='store_true')
    parser.add_argument('--sep', type=str, help="seperator for output csv file, default is tab character", default='\t')
    parser.add_argument('--processors', type=int, help="running processors, default is 1", default=1)
    parser.add_argument('-o', type=str, help="output base dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--cache-dir', type=str,
                        help=f'cache dir used for loading calls/bs-seq (speed up running), default is {global_cache_dir}',
                        default=global_cache_dir)
    parser.add_argument('--chrSet', nargs='+', help='chromosome list, default is human chromosome chr1-22, X and Y',
                        default=humanChrSet)
    parser.add_argument('--tagname', type=str, help="output unified file's tagname", default=None)
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    logger.debug(f"args={args}")
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()

    ## Set cache dir for each dataset
    if args.enable_cache or args.using_cache:
        ds_cache_dir = os.path.join(args.cache_dir, args.dsname)
        # os.makedirs(ds_cache_dir, exist_ok=True)
    else:
        ds_cache_dir = None

    if args.read_level_format:  ## if output read-level format, must read directly from raw file
        enable_cache = False
        using_cache = False
    else:
        enable_cache = args.enable_cache
        using_cache = args.using_cache

    tool_cutoff = 1
    bs_cutoff = 1

    # Currently we only use 1-base start format, for BED of singletons, non-singletons are use 1-base format
    baseFormat = 1

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    logger.debug(args)
    logger.debug(f'\n\n####################\n\n')

    if args.read_level_format:
        ## Output read-level unified format with 1-based start for ONT tools
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

        executor = ThreadPoolExecutor(max_workers=args.processors)
        for arg in input_list:
            executor.submit(import_and_save_read_level, *arg)
        executor.shutdown()

        save_done_file(out_dir)
        logger.info(f"Memory report: {get_current_memory_usage()}")
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
                                      using_cache=using_cache, enable_cache=enable_cache, cache_dir=ds_cache_dir)
            bgTruthList.append(bgTruth1)

        # Combine one/two replicates, using cutoff=1 or 5
        combine_bsdata = combineBGTruthList(bgTruthList, covCutoff=bs_cutoff)

        # Clean up bgTruthList
        del bgTruthList

        outfn = os.path.join(out_dir, f'{args.dsname}_BSseq-perSite-cov{bs_cutoff}.bed.gz')
        logger.debug(f'Combined BS-seq data (cov>={bs_cutoff}), all methylation level sites={len(combine_bsdata):,}')
        output_calldict_to_unified_bed_as_0base(combine_bsdata, outfn)

        # Clean up bgTruth, not used anymore
        del combine_bsdata

        logger.debug(f"Memory report: {get_current_memory_usage()}")
        logger.debug(f'\n\n####################\n\n')

    logger.debug("We are outputing bed CpG results for each tool")

    # callname -> # of sites
    ontcall_tools_dict = dict()

    input_list = []
    for callstr in args.calls:
        callencode, callfn = callstr.split(':')
        if len(callfn) == 0:
            continue
        callname = get_tool_name(callencode)

        outfn = os.path.join(out_dir, f'{args.dsname}_{callname}-perSite-cov{tool_cutoff}.bed.gz')
        input1 = (callfn, callname, callencode, tool_cutoff, outfn,)
        input_list.append(input1)

    executor = ThreadPoolExecutor(max_workers=args.processors)
    for arg in input_list:
        executor.submit(import_and_save_site_level, *arg)
    executor.shutdown()

    for key in ontcall_tools_dict.keys():
        logger.debug(f"tool={key}, sites={ontcall_tools_dict[key]}")
    save_done_file(out_dir)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    logger.info("TSS unified format generation DONE")
