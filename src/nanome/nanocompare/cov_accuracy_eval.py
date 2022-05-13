#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : site_level_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Generate site-level methylation correlation results in nanome paper.
"""

import argparse
import random

import pybedtools
from scipy import stats
from scipy.stats import PearsonRConstantInputWarning, pearsonr
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

from nanome.common.eval_common import *
from nanome.common.global_settings import get_tool_name, NANOME_VERSION, load_genome_annotation_config, save_done_file


def get_nsites_in_regions(callSet, bedfn, tagname):
    subset = filter_cpgkeys_using_bedfile(callSet, bedfn)
    ret = {tagname: len(subset)}
    return ret


def correlation_report_on_regions(corr_infn, bed_tuple_list, dsname=None, runid=None, outdir=None,
                                  large_mem=False,
                                  enable_base_detection_bedfile=enable_base_detection_bedfile,
                                  enable_cache=False, using_cache=False):
    """
    Calculate Pearson's correlation coefficient at different regions.
    :param corr_infn:
    :param beddir:
    :param dsname:
    :param outdir:
    :return:
    """
    global progress_bar_global_site
    progress_bar_global_site = tqdm(total=len(bed_tuple_list))
    progress_bar_global_site.set_description(f"MT-PCC-{dsname}-regions")

    executor = ThreadPoolExecutor(max_workers=args.processors)
    all_task = []
    ## input: df, bed_tuple
    ## return: list of dict [{tool1's pcc}, {toolk's pcc}]
    for bed_tuple in bed_tuple_list:
        future = executor.submit(compute_pcc_at_region, corr_infn, bed_tuple)
        future.add_done_callback(update_progress_bar_site_level)
        all_task.append(future)
    executor.shutdown()
    progress_bar_global_site.close()

    ret_list = []  # list of dict for dataframe
    for future in all_task:
        ret_l1 = future.result()
        if ret_l1 is None:
            continue
        ret_list.extend(ret_l1)

    # logger.info(dataset)
    outdf = pd.DataFrame(ret_list)
    logger.debug(outdf)

    outfn = os.path.join(outdir, f'{runid}_{dsname}.corrdata.coe.pvalue.each.regions.xlsx')
    outdf.to_excel(outfn)
    logger.debug(f'save to {outfn}')
    return outdf


def save_meth_corr_data(callresult_dict, bgTruth, reportCpGSet, outfn):
    """
    Save meth freq and cov results into csv file
    :param callresult_dict:
    :param bgTruth:
    :param reportCpGSet:
    :param outfn:
    :return:
    """
    # outfile = open(outfn, 'w')
    outfile = gzip.open(outfn, 'wt')

    header_list = ['chr', 'start', 'end', 'BGTruth_freq', 'BGTruth_cov', 'strand']

    # Output header
    for name in loaded_callname_list:
        header_list.extend([f'{name}_freq', f'{name}_cov'])
    outfile.write(sep.join(header_list))
    outfile.write("\n")

    for cpg in reportCpGSet:
        if baseFormat == 0:
            end = cpg[1] + 1
        elif baseFormat == 1:
            end = cpg[1]

        # Ouput ground truth
        row_list = [cpg[0], str(cpg[1]), str(end), f'{bgTruth[cpg][0]:.3f}', str(bgTruth[cpg][1]), cpg[2]]

        # Output each tool results
        for name in loaded_callname_list:
            if cpg in callresult_dict[name]:  # if cpg is in tool results
                row_list.extend([f'{callresult_dict[name][cpg][0]:.3f}', f'{callresult_dict[name][cpg][1]}'])
            else:  # if cpg key is not exist, we use NA as ''
                row_list.extend(['', ''])
        outfile.write(sep.join(row_list))
        outfile.write("\n")
    outfile.close()
    logger.debug(f"save to {outfn}\n")


def update_progress_bar_site_level(*a):
    """
    Update progress for multiprocessing
    :param a:
    :return:
    """
    global progress_bar_global_site
    progress_bar_global_site.update()


def get_num_intersect(callBed, region_bed, bedfn="", tagname=None):
    ret_bed = intersect_bed_regions(callBed, region_bed, bedfn=bedfn)
    if len(ret_bed) < 1:
        logger.error(
            f"Found 0 intersect CPG for tool with {tagname} region, callBed={callBed}, region_bed={region_bed}")
    ret = {tagname: len(ret_bed)}
    return ret


def compute_pcc_at_region(corr_infn, bed_tuple):
    infn, tagname, coord_bed = bed_tuple
    logger.debug(f'tagname={tagname}, coord_fn={infn}')
    if not args.large_mem and tagname != genome_wide_tagname and coord_bed is None:  # load on demand
        eval_coord_bed = get_region_bed_tuple(
            infn, enable_base_detection_bedfile=enable_base_detection_bedfile,
            enable_cache=args.enable_cache, using_cache=args.using_cache,
            cache_dir=ds_cache_dir)[2]
    else:  # large memory, or genome wide - None
        eval_coord_bed = coord_bed

    if tagname != genome_wide_tagname and eval_coord_bed is None:
        logger.debug(f"Region name={tagname} is not found, not compute PCC")
        return None

    df = pd.read_csv(corr_infn)
    newdf = filter_corrdata_df_by_bedfile(df, eval_coord_bed, infn)
    if newdf is None:
        logger.debug(f"Found intersection 0 CPGs for tagname={tagname}, no report for PCC")
        return None

    # Computer COE and pvalue
    newdf = newdf.filter(regex='_freq$', axis=1)
    ret_list = []
    for i in range(1, len(newdf.columns)):
        toolname = str(newdf.columns[i]).replace('_freq', '')
        try:  # too few samples will fail
            # with warnings.catch_warnings(): # not function
            warnings.filterwarnings('ignore', category=PearsonRConstantInputWarning)
            coe, pval = stats.pearsonr(newdf.iloc[:, 0], newdf.iloc[:, i])
        except:
            coe, pval = None, None

        # report to dataset
        ret = {
            'dsname': dsname,
            'Tool': toolname,
            'Location': tagname,
            '#Bases': len(newdf),
            'COE': coe,
            'p-value': pval
        }
        ret_list.append(ret)
    logger.debug(f"tagname={tagname}, pcc_return={ret_list}")
    return ret_list


def down_sample(call, k=5):
    ret_call = {}
    for key in call:
        ret_call[key] = random.sample(call[key], k)
    return ret_call


def pcc_mse_evaluation(call, bgTruth, bed_tuple=None, tool=None, down_sample_cov=None):
    joinedCPG = set(call.keys()).intersection(set(bgTruth.keys()))
    tagname = bed_tuple[1]
    if bed_tuple[2]:  # join bed with CPGs
        logger.info(f"Perform filter with bed region={tagname}")
        joinedCPG_bed = calldict2bed(joinedCPG)
        intersect_bed = intersect_bed_regions(joinedCPG_bed, bed_tuple[2], bed_tuple[0])
        joinedCPG = set(bedtxt2dict(intersect_bed).keys())
    logger.debug(f"Joined CPG={len(joinedCPG):,}, tagname={tagname}")
    if (len(joinedCPG)) < 1:
        return None, None

    cpg_list = []
    freq_list = []
    for cpg in joinedCPG:
        cpg_list.append(cpg)
        freq_list.append(bgTruth[cpg][0])
    cpg_df = pd.DataFrame({'cpg': cpg_list, 'freq': freq_list})

    ## Bin is [1,0.1], (0.1, 0.2], ..., (0.9, 1], class are 1-10
    cpg_df['bin_cat'] = pd.cut(cpg_df['freq'], bins=np.linspace(0, 1, args.bin_num + 1),
                               labels=list(range(1, args.bin_num + 1)), include_lowest=True)
    logger.debug(cpg_df)
    logger.debug(cpg_df['bin_cat'].value_counts())

    bg_freq = []
    pd_freq = []
    bin_bg_freq = defaultdict(list)  # category list
    bin_pd_freq = defaultdict(list)
    for index, row in cpg_df.iterrows():
        cpg = row['cpg']
        freq1 = row['freq']
        bin_cat = row['bin_cat']  # 1--10
        freq2 = sum(call[cpg]) / len(call[cpg])

        ## No bin results
        bg_freq.append(freq1)
        pd_freq.append(freq2)

        ## bin cat results
        bin_bg_freq[bin_cat].append(freq1)
        bin_pd_freq[bin_cat].append(freq2)

        bin_bg_freq[0].append(freq1)
        bin_pd_freq[0].append(freq2)

    logger.debug(f"bin_bg_freq key= {bin_bg_freq.keys()}, bin_pd_freq key={bin_pd_freq.keys()}")

    try:
        pcc, pcc_pvalue = pearsonr(bg_freq, pd_freq)
    except:
        pcc, pcc_pvalue = 0.0, 0.0
    mse = mean_squared_error(bg_freq, pd_freq)
    r2 = r2_score(bg_freq, pd_freq)
    mae = mean_absolute_error(bg_freq, pd_freq)

    ## No bin results
    ret = {
        'dsname': args.dsname,
        'Tool': tool, 'Location': tagname, 'Down-sample-cov': down_sample_cov,
        'Bases': len(joinedCPG),
        'PCC': pcc,
        'P-value': pcc_pvalue,
        'MAE': mae,
        'MSE': mse,
        'RMSE': np.sqrt(mse),
        'R2': r2,
    }

    ## Note, 0 - whole sets, 1-10, 10 bin results
    ret_bin_cat = []
    for k in range(0, args.bin_num + 1):
        bg_freqk = bin_bg_freq[k]
        pd_freqk = bin_pd_freq[k]
        try:
            pcc, pcc_pvalue = pearsonr(bg_freqk, pd_freqk)
            mse = mean_squared_error(bg_freqk, pd_freqk)
            rmse = np.sqrt(mse)
            r2 = r2_score(bg_freqk, pd_freqk)
            mae = mean_absolute_error(bg_freqk, pd_freqk)
        except:
            pcc = pcc_pvalue = mse = r2 = mae = rmse = None
        retk = {
            'dsname': args.dsname,
            'Tool': tool, 'Location': tagname, 'Down-sample-cov': down_sample_cov,
            'Bases': len(bg_freqk),
            'Bin-cat': k,
            'PCC': pcc,
            'P-value': pcc_pvalue,
            'MAE': mae,
            'MSE': mse,
            'RMSE': rmse,
            'R2': r2,

        }
        ret_bin_cat.append(retk)

    return ret, ret_bin_cat


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='cov_accuracy_eval (NANOME)',
                                     description='Coverage and performance evaluation')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix/output folder name, such as MethCorr-DS_WGBS_2reps",
                        required=True)
    parser.add_argument('--calls', nargs='+',
                        help='all ONT call results <tool-name>:<file-name> seperated by spaces, tool-name can be Nanopolish, Megalodon, DeepSignal, Guppy, Tombo, METEORE, DeepMod',
                        required=True)
    parser.add_argument('--bgtruth', type=str,
                        help="background truth file <encode-type>:<file-name1>;<file-name1>, encode-type can be 'encode' or 'bismark'",
                        default=None)
    parser.add_argument('--genome-annotation', type=str,
                        help='genome annotation dir, contain BED files',
                        default=None)
    parser.add_argument('--beddir', type=str,
                        help="base dir for concordant/discordant BED files generated by read-level analysis, make sure the dsname is same",
                        default=None)
    parser.add_argument('--min-bgtruth-cov', type=int, help="cutoff for coverage in bg-truth, default is >=5",
                        default=5)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff for coverage in nanopore tools, default is >=3",
                        default=3)
    parser.add_argument('--seed', type=int, help="random seed, default is 53",
                        default=53)
    parser.add_argument('--bin-num', type=int, help="bin number for bin analysis, default is 10",
                        default=10)
    parser.add_argument('--bin-analysis', help="bin number for bin analysis, default is 10",
                        action='store_true')
    parser.add_argument('--startk', type=int, help="start coverage for sampling, default is 5", default=5)
    parser.add_argument('--stepk', type=int, help="step for coverage sampling, default is 5", default=5)

    parser.add_argument('--chrSet', nargs='+', help='chromosome list, default is human chr1-22, X and Y',
                        default=HUMAN_CHR_SET)
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=',')
    parser.add_argument('--processors', type=int, help="number of processors used, default is 1", default=1)
    parser.add_argument('-o', type=str, help=f"output base dir, default is {pic_base_dir}", default=pic_base_dir)
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--bedtools-tmp', type=str, help=f'bedtools temp dir, default is {global_temp_dir}',
                        default=global_temp_dir)
    parser.add_argument('--cache-dir', type=str,
                        help=f'cache dir used for loading calls/bs-seq(speed up running), default is {global_cache_dir}',
                        default=global_cache_dir)
    parser.add_argument('--disable-bed-check',
                        help="if disable auto-checking the 0/1 base format for genome annotations",
                        action='store_true')
    parser.add_argument('--mpi',
                        help="if using multi-processing/threading for evaluation, it can speed-up but need more memory",
                        action='store_true')
    parser.add_argument('--config', help="if print out config file for genome annotation", action='store_true')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()

    random.seed(args.seed)
    dsname = args.dsname
    ## Set tmp dir for bedtools, each process use a bed tmp dir
    ## because the tmp dir files may be cleaned by the end of the process
    bed_temp_dir = os.path.join(args.bedtools_tmp, dsname)
    os.makedirs(bed_temp_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(bed_temp_dir)

    ## Set cache dir for each dataset
    if args.enable_cache or args.using_cache:
        ds_cache_dir = os.path.join(args.cache_dir, dsname)
        # os.makedirs(ds_cache_dir, exist_ok=True)
    else:
        ds_cache_dir = None

    # cache function same with read level
    enable_cache = args.enable_cache
    using_cache = args.using_cache

    # runid is always like 'MethCorr-K562_WGBS_2Reps', remove first word as RunPrefix like K562_WGBS_2Reps
    RunPrefix = args.runid.replace('CovAccuracy-', '')

    # tool coverage cutoff 1, or 3, 5
    minToolCovCutt = args.toolcov_cutoff

    # bgtruth coverage cutoff 1, or 5, 10  --min-bgtruth-cov
    bgtruthCutt = args.min_bgtruth_cov

    # We import and report use 1-base start format
    baseFormat = 1

    # output csv seperator: , or tab
    sep = args.sep

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    # Add logging files also to result output dir
    add_logging_file(os.path.join(out_dir, 'run-results.log'))
    logger.debug(args)

    if args.config:
        load_genome_annotation_config(verbose=True)
    logger.debug(f'\n\n####################\n\n')

    if args.bgtruth:
        # we import multiple (1 or 2) replicates and join them
        encode, fnlist = args.bgtruth.split(':')
        fnlist = fnlist.split(';')

        if len(fnlist) > 2:
            raise Exception(f'Currently only support the number of bgtruth upto two, but found more: {fnlist}')

        logger.debug(f'BGTruth fnlist={fnlist}, encode={encode}')

        bgTruthList = []
        for fn in fnlist:
            if len(fn) == 0:  # incase of input like 'bismark:/a/b/c;'
                continue
            # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
            bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True,
                                      filterChr=args.chrSet,
                                      using_cache=using_cache, enable_cache=enable_cache, cache_dir=ds_cache_dir)
            bgTruthList.append(bgTruth1)

        # Combine one/two replicates, using cutoff=1 or 5
        bgTruth = combineBGTruthList(bgTruthList, covCutoff=bgtruthCutt)

        logger.info(f'Combined BS-seq data (cov>={bgtruthCutt}), all methylation level sites={len(bgTruth):,}')
        logger.debug(f"Memory report: {get_current_memory_usage()}")
        logger.debug(f'\n\n####################\n\n')
    else:
        bgTruth = None

    callfn_dict = defaultdict()  # callname -> filename

    # callname -> [call0, call1], call0 is no-filter results, call1 is filter by cutoff, and convert to [meth-freq, meth-cov] results.
    callresult_dict_cov3 = defaultdict()
    loaded_callname_list = []

    for callstr in args.calls:
        callencode, callfn = callstr.split(':')

        if len(callfn) == 0:
            continue

        callname = get_tool_name(callencode)
        callfn_dict[callname] = callfn

        loaded_callname_list.append(callname)

        # For site level evaluation, only need (freq, cov) results, no score needed. Especially for DeepMod, we must import as freq and cov format from DeepMod.Cluster encode
        # Do not filter bgtruth, because we use later for overlapping (without bg-truth)
        call1 = import_call(callfn, callencode, baseFormat=baseFormat, filterChr=args.chrSet,
                            enable_cache=enable_cache, using_cache=using_cache,
                            include_score=False, siteLevel=False, cache_dir=ds_cache_dir)
        call3 = filter_cpg_dict_by_cov(call1, coverage=minToolCovCutt)
        callresult_dict_cov3[callname] = call3
        logger.info(f"Tool={callname}, cov=1, CPG={len(call1):,}; cov={minToolCovCutt}, CPG={len(call3):,}; ")
        logger.info(f"Memory report: {get_current_memory_usage()}")

    logger.debug(f'\n\n####################\n\n')

    annot_dir = args.genome_annotation if args.genome_annotation is not None else '.'
    regions_full_filepath = [None] + [os.path.join(annot_dir, cofn) for cofn in region_filename_dict.keys()]
    region_bed_list = [(infn, get_region_tagname(infn), None,)
                       for infn in regions_full_filepath]

    ## Downsampleing to different coverage for tools: 5,10,15,20,25
    dataset = []
    dataset_bin_cat = []

    for bed_tuple in region_bed_list:
        infn, tagname, coord_bed = bed_tuple
        if tagname != genome_wide_tagname and coord_bed is None:  # load on demand
            eval_coord_bed = get_region_bed_tuple(
                infn, enable_base_detection_bedfile=enable_base_detection_bedfile,
                enable_cache=args.enable_cache, using_cache=args.using_cache,
                cache_dir=ds_cache_dir)[2]
            if eval_coord_bed is None:
                continue
        else:
            eval_coord_bed = coord_bed
        eval_coord_bed_tuple = (infn, tagname, eval_coord_bed,)
        for down_sample_cov in range(args.startk, minToolCovCutt + 1, args.stepk):
            for tool in callresult_dict_cov3:
                call = callresult_dict_cov3[tool]
                # down sample the orginal 25X call into 5, 10, etc.
                down_sample_call = down_sample(call, down_sample_cov)
                logger.debug(f"After downsample to cov={down_sample_cov}, Tool={tool}")
                ret, ret_bin_cat_list = pcc_mse_evaluation(down_sample_call, bgTruth, eval_coord_bed_tuple, tool=tool,
                                                           down_sample_cov=down_sample_cov)
                if not ret:
                    continue
                logger.debug(f"Tool={tool}, tagname={tagname}, down-sample cov={down_sample_cov}, ret={ret}")
                logger.debug(f"\n\nbin_cat results: ret_bin_cat_list={ret_bin_cat_list}")
                dataset.append(ret)
                dataset_bin_cat += ret_bin_cat_list

    outdf = pd.DataFrame(dataset)
    logger.info(outdf)

    outdf2 = pd.DataFrame(dataset_bin_cat)
    logger.info(outdf2)

    outfn = os.path.join(out_dir,
                         f"Coverage_vs_accuracy_{args.dsname}_bgtruth{bgtruthCutt}_toolcov{minToolCovCutt}.csv")
    outdf.to_csv(outfn)
    logger.info(f"save to {outfn}")

    outfn2 = os.path.join(out_dir,
                          f"Coverage_vs_accuracy_{args.dsname}_bgtruth{bgtruthCutt}_toolcov{minToolCovCutt}_bincat{args.bin_num}.csv")
    outdf2.to_csv(outfn2)
    logger.info(f"save to {outfn2}")

    save_done_file(out_dir)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    logger.info("### Coverage accuracy analysis DONE")
