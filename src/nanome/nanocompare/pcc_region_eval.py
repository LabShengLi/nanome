#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : pcc_region_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Evaluate PCC at different regions in nanome paper
"""

import argparse

import pybedtools
from scipy import stats
from sklearn.metrics import mean_squared_error

from nanome.common.eval_common import *
from nanome.common.global_settings import NANOME_VERSION, load_genome_annotation_config, save_done_file


def update_progress_bar_site_level(*a):
    """
    Update progress for multiprocessing
    :param a:
    :return:
    """
    global progress_bar_global_site
    progress_bar_global_site.update()


def correlation_report_on_regions(corr_infn, bed_tuple_list, dsname=None, outdir=None,
                                  mpi=True):
    """
    Calculate Pearson's correlation coefficient at different regions.
    :param corr_infn:
    :param beddir:
    :param dsname:
    :param outdir:
    :return:
    """
    df = pd.read_csv(corr_infn, index_col=False)
    logger.info(f"df={df}")

    ## Correlation matrix
    # Report correlation matrix for joined results
    df2 = df.filter(regex='_freq$', axis=1)
    cordf = df2.corr()
    # Count CpGs
    num_join_with_bsseq = [len(df2.iloc[:, 0])]
    for k in range(1, len(df2.columns)):
        num_join_with_bsseq.append(df2.iloc[:, k].notna().sum())
    cordf = pd.concat(
        [cordf, pd.Series(num_join_with_bsseq, index=cordf.index).rename('CpGs_with_BSseq')],
        axis=1)

    cordf.columns = [col_name.replace('_freq', '') for col_name in cordf.columns]
    cordf.index = [idx_name.replace('_freq', '') for idx_name in cordf.index]
    logger.debug(f'Correlation matrix is:\n{cordf}')
    corr_outfn = os.path.join(outdir,
                              f"{args.runid}.{os.path.basename(args.i).replace('.csv.gz', '')}.correlation.matrix.xlsx")
    cordf.to_excel(corr_outfn)
    logger.debug(f"save to {corr_outfn}")

    ## Evaluate on regions
    # if df.isnull().values.any():
    #     df.fillna('.', inplace=True)

    global progress_bar_global_site
    progress_bar_global_site = tqdm(total=len(bed_tuple_list))
    progress_bar_global_site.set_description(f"MT-PCC-{dsname}-regions")

    if mpi:
        executor = ThreadPoolExecutor(max_workers=args.processors)
    else:
        executor = ThreadPoolExecutor(max_workers=1)
    all_task = []

    for bed_tuple in bed_tuple_list:
        future = executor.submit(compute_pcc_at_region, df, bed_tuple)
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

    outfn = os.path.join(outdir, f"{args.runid}.{os.path.basename(args.i).replace('.csv.gz', '')}.pcc.regions.xlsx")
    outdf.to_excel(outfn)
    logger.debug(f'save to {outfn}')
    return outdf


def compute_pcc_at_region(df, bed_tuple):
    """
    Compute PCC of input DF with respect to bed region tuple
    Args:
        df:
        bed_tuple:

    Returns:

    """
    infn, tagname, coord_bed = bed_tuple
    logger.debug(f'tagname={tagname}, coord_fn={infn}')
    if not args.large_mem and tagname != genome_wide_tagname and coord_bed is None:  # load on demand
        eval_coord_bed = get_region_bed_tuple(
            infn, enable_base_detection_bedfile=not args.disable_bed_check,
            enable_cache=args.enable_cache, using_cache=args.using_cache,
            cache_dir=ds_cache_dir)[2]
    else:  # large memory, or genome wide - None
        eval_coord_bed = coord_bed

    if tagname != genome_wide_tagname and eval_coord_bed is None:
        logger.debug(f"Region name={tagname} is not found, not compute PCC")
        return None

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
            # warnings.filterwarnings('ignore', category=PearsonRConstantInputWarning)
            mask_notna = newdf.iloc[:, i].notna().values
            coe, pval = stats.pearsonr(newdf.iloc[mask_notna, 0], newdf.iloc[mask_notna, i].astype(np.float64))
            mse = mean_squared_error(newdf.iloc[mask_notna, 0], newdf.iloc[mask_notna, i].astype(np.float64))
        except:  ## May be pearsonr function failed
            coe, pval = None, None
            mse = None

        # report to dataset
        ret = {
            'dsname': dsname,
            'Tool': toolname,
            'Location': tagname,
            '#Bases': newdf.iloc[:, i].notna().sum(),
            'MSE': mse,
            'COE': coe,
            'p-value': pval
        }
        ret_list.append(ret)
    logger.debug(f"tagname={tagname}, pcc_return={ret_list}")
    return ret_list


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='pcc_region_eval (NANOME)',
                                     description='Site-level PCC correlation at different genomic regions')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('-i', type=str, help="input freq file for BS-seq and tools", required=True)
    parser.add_argument('--runid', type=str,
                        help="running prefix/output folder name, such as PCCRegion-Dataset_WGBS_2Reps",
                        required=True)
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
    parser.add_argument('--chrSet', nargs='+', help='chromosome list, default is human chr1-22, X and Y',
                        default=HUMAN_CHR_SET)
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=',')
    parser.add_argument('--processors', type=int, help="number of processors used, default is 1", default=1)

    parser.add_argument('-o', type=str, help=f"output base dir, default is {pic_base_dir}", default=pic_base_dir)
    parser.add_argument('--report-no-join', help="if output no-join report also", action='store_true')
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--bedtools-tmp', type=str, help=f'bedtools temp dir, default is {global_temp_dir}',
                        default=global_temp_dir)
    parser.add_argument('--cache-dir', type=str,
                        help=f'cache dir used for loading calls/bs-seq(speed up running), default is {global_cache_dir}',
                        default=global_cache_dir)
    parser.add_argument('--large-mem', help="if using large memory (>100GB) for speed up", action='store_true')
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

    dsname = args.dsname
    ## Set tmp dir for bedtools, each process use a bed tmp dir
    ## because the tmp dir files may be cleaned by the end of the process
    bed_temp_dir = os.path.join(args.bedtools_tmp, f"{dsname}_pcc_region")
    os.makedirs(bed_temp_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(bed_temp_dir)

    ## Set cache dir for each dataset
    if args.enable_cache or args.using_cache:
        ds_cache_dir = os.path.join(args.cache_dir, f"{dsname}_pcc_region")
        # os.makedirs(ds_cache_dir, exist_ok=True)
    else:
        ds_cache_dir = None

    # cache function same with read level
    enable_cache = args.enable_cache
    using_cache = args.using_cache

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

    ## load region bed list
    logger.debug("Create region bed list firstly, take times......")

    # Evaluated all region filename lists, bed objects
    # assume all files are located in args.genome_annotation dir
    annot_dir = args.genome_annotation if args.genome_annotation is not None else '.'

    # regions_full_filepath = [os.path.join(annot_dir, cofn) for cofn in narrowCoordNameList[1:]] + \
    #                         [os.path.join(annot_dir, cofn) for cofn in cg_density_coord_name_list] + \
    #                         [os.path.join(annot_dir, cofn) for cofn in rep_coord_name_list]
    # region file path from genome-wide, singletons, to genic/intergenic, cg-density, and repetitive, the concordant and discoradnt wil be discovered later
    regions_full_filepath = [None] + [os.path.join(annot_dir, cofn) for cofn in region_filename_dict.keys()]

    if args.large_mem:  # load all in memory
        region_bed_list = get_region_bed_pairs_list_mp(
            regions_full_filepath,
            processors=args.processors,
            enable_base_detection_bedfile=not args.disable_bed_check,
            enable_cache=args.enable_cache,
            using_cache=args.using_cache,
            cache_dir=ds_cache_dir)
        logger.info(f"Memory report: {get_current_memory_usage()}")
    else:  # load bed coord later
        region_bed_list = [(infn, get_region_tagname(infn), None,)
                           for infn in regions_full_filepath]

    if args.beddir:  # add concordant and discordant region coverage if needed
        logger.debug(f'We are finding Concordant and Discordant BED file at basedir={args.beddir}')
        concordantFileName = find_bed_filename(basedir=args.beddir,
                                               pattern=f'*{args.dsname}*.concordant.bed.gz')
        concordant_bed = get_region_bed(concordantFileName) if concordantFileName is not None else None

        discordantFileName = find_bed_filename(basedir=args.beddir,
                                               pattern=f'*{args.dsname}*.discordant.bed.gz')
        discordant_bed = get_region_bed(discordantFileName) if discordantFileName is not None else None
    else:
        concordant_bed = None
        discordant_bed = None
    ## Add concordant/discordant if possible
    if concordant_bed is not None:
        region_bed_list += [(concordantFileName, 'Concordant', concordant_bed,)]
    if discordant_bed is not None:
        region_bed_list += [(discordantFileName, 'Discordant', discordant_bed,)]

    logger.debug(f"Evaluated on regions: {region_bed_list}")

    eval_genomic_context_tuple = region_bed_list

    # file like: Meth_corr_plot_data_joined-TestData_RRBS_2Reps-bsCov1-minToolCov1-baseFormat1.sorted.csv.gz
    dsname = args.dsname

    logger.info(f"Start report PCC in different genomic regions based on file={args.i}, dsname={args.dsname}")
    correlation_report_on_regions(
        args.i, bed_tuple_list=eval_genomic_context_tuple, dsname=args.dsname,
        outdir=out_dir, mpi=args.mpi)
    logger.debug(f"Memory report: {get_current_memory_usage()}")

    save_done_file(out_dir)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    logger.info("### PCC region evaluation DONE")
