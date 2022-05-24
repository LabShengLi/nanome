#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : read_level_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
This is a new proposed evaluation criteria:
we join all exactly same predictions (read-id, loci) of tools at read-level firstly, then evaluate the per-read/per-site performance.
"""

import argparse
import os.path

import pybedtools
from scipy import stats
from sklearn.metrics import accuracy_score, roc_auc_score, mean_squared_error

from nanome.common.eval_common import *
from nanome.common.global_settings import NANOME_VERSION, save_done_file, \
    region_filename_dict

SITE_COLUMNS = ['Chr', 'Pos', 'Strand']
READ_COLUMNS = ['ID'] + SITE_COLUMNS

DTYPE = {'ID': str, 'Chr': str, 'Pos': int, 'Strand': str,
         'Freq': float, 'Coverage': int}


def prob_to_llr_2(meth_prob):
    """
    METEORE, NANOME manner
    Args:
        meth_prob:

    Returns:

    """
    return math.log2((meth_prob + EPSLONG) / (1 - meth_prob + EPSLONG))


def prob_to_llr_e(meth_prob):
    """
    Megalodon manner:
    prob_to_llr_2(0.8)
    1.999945900626566
    prob_to_llr_e(0.8)
    1.3862568622917248

    Args:
        meth_prob:

    Returns:

    """
    return math.log((meth_prob + EPSLONG) / (1 - meth_prob + EPSLONG))


DEFAULT_CUTOFF_MEGALODON = prob_to_llr_e(0.8)
DEFAULT_CUTOFF_NANOPOLISH = 2


def update_progress_bar_join_preds_eval(*a):
    """
    Update progress for multithreading/multiprocessing
    :param a:
    :return:
    """
    global progress_bar_global_join_preds
    progress_bar_global_join_preds.update()


def report_read_level_performance(df, toolList=None, dsname=None, outdir=None, location="Genome-wide", if_save=False):
    """
    Report read-level performance on joined predictions
    Args:
        df:
        toolList:
        outdir:
        location:
        if_save:

    Returns:

    """
    logger.debug("Evaluate on data...")
    logger.debug(f"Total predictions:{len(df):,}")

    numSites = len(df.drop_duplicates(subset=SITE_COLUMNS))
    logger.debug(f"Total CpGs:{numSites:,}")

    y_truth = df['Freq'].apply(freq_to_label).astype(int)

    dataset = defaultdict(list)
    for tool in toolList:
        y_score = df[tool]
        y_pred = df[tool].apply(lambda x: 1 if x > 0 else 0)
        if len(y_score) != len(y_truth):
            raise Exception(f"Not assert for length of y_truth and y_score")

        accuracy = accuracy_score(y_truth, y_pred)
        precision = precision_score(y_truth, y_pred)
        recall = recall_score(y_truth, y_pred)
        f1 = f1_score(y_truth, y_pred)
        roc_auc = roc_auc_score(y_truth, y_score)
        dataset['Dataset'].append(dsname)
        dataset['Tool'].append(tool)
        dataset['Location'].append(location)
        dataset['Accuracy'].append(accuracy)
        dataset['Precision'].append(precision)
        dataset['Recall'].append(recall)
        dataset['F1'].append(f1)
        dataset['ROC_AUC'].append(roc_auc)
        dataset['#Bases'].append(numSites)
        dataset['#Pred'].append(len(y_truth))

    outdf = pd.DataFrame.from_dict(dataset)
    logger.debug(outdf)

    if if_save:
        outfn = os.path.join(outdir, f'{args.dsname}_read_level_joined_preds_perf_{location}.csv')
        outdf.to_csv(outfn)
        logger.info(f"save to {outfn}")
    return outdf


def load_read_unified_df(infn, chrSet, cutoff=0.0, test=None, chunksize=1000000):
    """
    Load one read-level data into DF, with chr, and cutoff filter
    Args:
        infn:
        chrSet:
        cutoff:
        test:
        chunksize:

    Returns:

    """
    logger.debug(f"Start load file:{infn}, chrs={chrSet}")
    iter_df = pd.read_csv(infn, header=0, index_col=False, sep="\t", iterator=True,
                          chunksize=chunksize, nrows=test)
    df = pd.concat([chunk[chunk['Chr'].isin(chrSet)] for chunk in iter_df])
    df.drop_duplicates(subset=READ_COLUMNS, inplace=True)
    if cutoff > EPSLONG:  # has cutoff, such as 2, 2.5
        df = df[df['Score'].abs >= cutoff]
    return df


def joined_preds_for_all_tools(callDict, chrFilter, dsname="dsname", bsseq_df=None, test=None,
                               is_save=False, outdir=None, chunksize=1000000):
    """
    Join the unified read-level for all tools by chromosome
    Args:
        callDict:
        chrFilter:
        dsname:
        bsseq_df:
        test:
        is_save:
        outdir:

    Returns:

    """
    chrs = set(chrFilter)

    dfall = None
    for toolName in callDict:
        df1 = load_read_unified_df(callDict[toolName], chrs, test=test,
                                   chunksize=chunksize)
        df1.rename(columns={'Score': toolName}, inplace=True)
        logger.debug(f"df1={df1}")
        if dfall is None:
            dfall = df1
        else:
            dfall = dfall.merge(df1, on=READ_COLUMNS, how='inner')

    if bsseq_df is not None:
        dfall = dfall.merge(bsseq_df, on=SITE_COLUMNS, how='inner')
    logger.debug(f"dfall={dfall}")
    if is_save and outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        chrListStr = '_'.join(chrs)
        outfn = os.path.join(outdir, f'{dsname}_T{len(callDict.keys())}_{chrListStr}_join_preds_eval_db.csv.gz')
        dfall.to_csv(outfn, index=False)
        logger.debug(f"save to {outfn}")

        ## sort db
        sort_per_read_csv_file(outfn, outfn.replace('csv.gz', 'sort.csv.gz'), deduplicate=True)
        os.remove(outfn)
        outfn = outfn.replace('csv.gz', 'sort.csv.gz')
        logger.debug(f"sort and save to {outfn}")

    return True


def find_join_preds_bgtruth_as_df(dsname, chrs, db_dir, cutoffDict, dtype=None, certainSites=False,
                                  fully_meth_threshold=1.0,
                                  sel_cols=None,
                                  is_save=False,
                                  outdir=None):
    """
    Find all files for joined predictions, and create the DF
    Args:
        dsname:
        chrs:
        db_dir:
        dtype:
        fully_meth_threshold:
        sel_cols:
        is_save:
        outdir:

    Returns:

    """
    find_files = []
    for chr in chrs:
        fnlist = glob.glob(os.path.join(db_dir, '**', f'{dsname}_T*_{chr}_join_preds_eval_db.sort.csv.gz'),
                           recursive=True)
        if len(fnlist) < 1:
            logger.warn(f"Not found {chr} for {dsname} at {db_dir}")
        elif len(fnlist) > 1:
            logger.warn(
                f"Found more than one file for {chr} for {dsname} at {db_dir}: {fnlist}, we use only first one.")
        if len(fnlist) >= 1:
            find_files.append(fnlist[0])
    logger.debug(f"find_files={find_files}, len={len(find_files)}")

    dflist = []
    for infn in find_files:
        try:
            df1 = pd.read_csv(infn, header=0, index_col=None, dtype=dtype)
        except:
            logger.error(f"File read fail: {infn}")
            continue
        if certainSites:
            df1 = df1[(df1['Freq'] <= EPSLONG) | (df1['Freq'] >= fully_meth_threshold - EPSLONG)]
        if sel_cols is not None:
            df1 = df1[sel_cols]
        ## DeepSignal may contain na, drop these rows
        df1.replace([np.inf, -np.inf], np.nan, inplace=True)
        df1.dropna(inplace=True)

        ## apply tool cutoff
        for toolName in cutoffDict.keys():
            cutoff = cutoffDict[toolName]
            if cutoff > EPSLONG:
                df1 = df1[df1[toolName].abs() >= cutoff]

        dflist.append(df1)
    df = pd.concat(dflist)
    df.reset_index(drop=True, inplace=True)
    if is_save:
        outdirf = os.path.join(outdir, 'eval_db')
        os.makedirs(outdirf, exist_ok=True)
        if certainSites:
            tagname = "certain_cpgs"
        else:
            tagname = "all_cpgs"
        outfn = os.path.join(outdirf, f'{dsname}_eval_preds_{tagname}_db.csv.gz')
        df.to_csv(outfn, index=False)
    return df


def filter_preds_df_by_bedfile(df, coord_bed, coord_fn):
    """
    Filter lines in correlation data, within coordinate BED file
    :param df:
    :param coord_fn:
    :return:
    """
    if coord_bed is None:  # No need to filter for genome wide
        return df
    # suppress warnings /home/liuya/anaconda3/envs/nanocompare/lib/python3.6/site-packages/pybedtools/bedtool.py:3287: UserWarning: Default names for filetype bed are:
    # ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    # but file has 20 fields; you can supply custom names with the `names` kwarg
    #   % (self.file_type, _names, self.field_count()))
    warnings.filterwarnings('ignore', category=UserWarning)

    ## In order to support NaN in dataframe
    ## Need replace NaN to ., then for sort, or else will encounter error: Differing number of BED fields encountered at line: 1339
    ## Ref: NA values https://github.com/daler/pybedtools/issues/257
    bed_of_df = BedTool.from_dataframe(df).sort()

    if not isinstance(coord_bed, pybedtools.BedTool) or not isinstance(bed_of_df, pybedtools.BedTool):
        logger.error(f"bed_of_df={type(bed_of_df)}, coord_bed={type(coord_bed)}")
        raise Exception(f"The bed object is not correct")

    bed_of_intersect = intersect_bed_regions(bed_of_df, coord_bed, coord_fn)
    ## logger.debug(f"bed_of_df={len(bed_of_df)}, bed_of_intersect={len(bed_of_intersect)}")
    if len(bed_of_intersect) > 0:
        ### retdf = bed_of_intersect.to_dataframe().replace('.', np.NaN)
        ### MUST set column names here, or else there is one row missing!!!
        retdf = bed_of_intersect.to_dataframe(names=df.columns)
        ## logger.debug(f"retdf={len(retdf)}")
    else:
        retdf = None
    return retdf


def create_region_bed_list(genome_annotation, beddir, dsname):
    """
    Create a region tuple list for genomic region evaluation
    Args:
        genome_annotation:
        beddir:
        dsname:

    Returns:

    """
    annot_dir = genome_annotation if genome_annotation is not None else '.'
    regions_full_filepath = [None] + [os.path.join(annot_dir, cofn) for cofn in region_filename_dict.keys()
                                      if os.path.exists(os.path.join(annot_dir, cofn))]
    region_bed_list = [(infn, get_region_tagname(infn), None,)
                       for infn in regions_full_filepath]
    if beddir:  # add concordant and discordant region coverage if needed
        logger.debug(f'We are finding Concordant and Discordant BED file at basedir={beddir}')
        concordantFileName = find_bed_filename(basedir=beddir,
                                               pattern=f'*{dsname}*.concordant.bed.gz')
        # concordant_bed = get_region_bed(concordantFileName) if concordantFileName is not None else None

        discordantFileName = find_bed_filename(basedir=beddir,
                                               pattern=f'*{dsname}*.discordant.bed.gz')
        # discordant_bed = get_region_bed(discordantFileName) if discordantFileName is not None else None
        ## Add concordant/discordant if possible
        if concordantFileName is not None:
            region_bed_list += [(concordantFileName, 'Concordant', None,)]
        if discordantFileName is not None:
            region_bed_list += [(discordantFileName, 'Discordant', None,)]

    logger.debug(f"Evaluated on regions: {region_bed_list}, len={len(region_bed_list)}")
    return region_bed_list


def read_to_site_df(read_df, toolList, dsname, tool_cov, if_save=False, outdir=None):
    """
    Group read level CpGs, freq for tools, and bs-seq
    Args:
        read_df:

    Returns:

    """
    df = read_df.copy()
    for tool in toolList:
        df[tool] = df[tool].apply(lambda x: 1 if x > 0 else 0)

    agg_func = {tool: 'sum' for tool in toolList}
    agg_func.update({'Freq': 'first', 'ID': 'count'})
    site_df = df[["Chr", "Pos", "Strand", "ID"] + toolList + ["Freq"]].groupby(by=["Chr", "Pos", "Strand"]).agg(
        agg_func)
    site_df.rename(columns={"ID": "Reads"}, inplace=True)

    for tool in toolList:
        site_df[tool] = site_df[tool] / site_df['Reads']
    site_df = site_df[site_df['Reads'] >= tool_cov]
    site_df.reset_index(inplace=True)
    logger.debug(site_df)

    if if_save:
        outdirf = os.path.join(outdir, 'eval_db')
        os.makedirs(outdirf, exist_ok=True)
        outfn = os.path.join(outdirf, f"{dsname}_eval_sites_db.csv.gz")
        site_df.to_csv(outfn, index=False)
        logger.info(f"save to {outfn}")
    return site_df


def report_site_level_performance(site_df, toolList, dsname=None, location="Genome-wide", if_save=False, outdir=None):
    """
    Report site level PCC and MSE
    Args:
        site_df:
        toolList:
        location:
        if_save:
        outdir:

    Returns:

    """
    dataset = defaultdict(list)
    for tool in toolList:
        coe, pval = stats.pearsonr(site_df['Freq'], site_df[tool])
        mse = mean_squared_error(site_df['Freq'], site_df[tool])
        dataset['Dataset'].append(dsname)
        dataset['Tool'].append(tool)
        dataset['Location'].append(location)
        dataset['PCC'].append(coe)
        dataset['P_value'].append(pval)
        dataset['MSE'].append(mse)
        dataset['#Bases'].append(len(site_df))
        dataset['#Pred'].append(site_df['Reads'].sum())

    site_perf_df = pd.DataFrame.from_dict(dataset)
    logger.debug(f"site_perf_df={site_perf_df}")

    if if_save:
        outfn = os.path.join(outdir, f"{args.dsname}_site_level_joined_preds_perf_{location}.csv")
        site_perf_df.to_csv(outfn, index=False)
        logger.info(f"save to {outfn}")

    return site_perf_df


def eval_read_level_at_region(infn, regionName, read_df1, dsname, callDict):
    """
    Used for submit to multi-threading evaluation at a region for read-level
    Args:
        infn:
        regionName:
        read_df1:
        dsname:
        callDict:

    Returns:

    """
    if regionName != genome_wide_tagname:
        eval_coord_bed = get_region_bed_tuple(infn)[2]
        if eval_coord_bed is None:
            logger.warn(f"Region name={regionName} is not found, not compute read-level")
            return None
    if regionName == genome_wide_tagname:
        eval_df = read_df1
    else:
        eval_df = filter_preds_df_by_bedfile(read_df1, eval_coord_bed, infn)
    if eval_df is None:
        logger.warn(f"No interset rows for {regionName}")
        return None
    df1 = report_read_level_performance(eval_df, toolList=list(callDict.keys()), dsname=dsname,
                                        location=regionName)
    return df1


def eval_site_level_at_region(infn, regionName, site_df1, dsname, callDict):
    """
    Used for submit to multi-threading evaluation at a region for site-level
    Args:
        infn:
        regionName:
        site_df1:
        dsname:
        callDict:

    Returns:

    """
    if regionName != genome_wide_tagname:
        eval_coord_bed = get_region_bed_tuple(infn)[2]
        if eval_coord_bed is None:
            logger.warn(f"Region name={regionName} is not found, not compute read-level")
            return None
    if regionName == genome_wide_tagname:
        eval_df = site_df1
    else:
        eval_df = filter_preds_df_by_bedfile(site_df1, eval_coord_bed, infn)

    if eval_df is None:
        logger.warn(f"No interset rows for {regionName}")
        return None
    df1 = report_site_level_performance(eval_df, toolList=list(callDict.keys()), dsname=dsname,
                                        location=regionName)
    return df1


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='join_preds_eval (NANOME)',
                                     description='Performance evaluation for joined predictions by all tools')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str,
                        help="running prefix/output folder name, such as JoinPredsEval-Dataset_WGBS_2Reps",
                        required=True)
    parser.add_argument('--calls', nargs='+',
                        help='all ONT call results <tool-name>:<unified-per-read-file> seperated by space, tool-name can be any name for the file, the file MUST be encoded in unified read-level format by NANOME definition',
                        required=True)
    parser.add_argument('--bs-seq-bed', type=str,
                        help="background truth file of sorted BED unified format",
                        default=None)
    parser.add_argument('--genome-annotation', type=str,
                        help='genome annotation dir, contain BED files',
                        default=None)
    parser.add_argument('--min-bgtruth-cov', type=int, help="min bg-truth coverage cutoff, default is %(default)s",
                        default=5)
    parser.add_argument('--toolcov-site', type=int,
                        help="cutoff for coverage in nanopore tools of site level evaluation, default is >=%(default)s",
                        default=3)
    parser.add_argument('--chunksize', type=int, help="min bg-truth coverage cutoff, default is %(default)s",
                        default=1000000)
    parser.add_argument('--processors', type=int, help="number of processors used, default is %(default)s", default=1)
    parser.add_argument('--chrs', nargs='+', type=str, help='chromosome list, default is human chr1-22, X and Y',
                        default=HUMAN_CHR_SET)
    parser.add_argument('--score-cutoff', nargs='+', type=float,
                        help='cutoff of LLR score for listed tool in --calls, default is 0 for all/not specified',
                        default=None)
    parser.add_argument('--fully-meth-threshold', type=float, default=1.0,
                        help='fully methylated threshold (e.g., 0.9), default is 1.0')
    parser.add_argument('-o', type=str, help=f"output base dir, default is {pic_base_dir}", default=pic_base_dir)
    parser.add_argument('--bedtools-tmp', type=str, help=f'bedtools temp dir, default is %(default)s',
                        default=global_temp_dir)
    parser.add_argument('--disable-bed-check',
                        help="if disable auto-checking the 0/1 base format for genome annotations",
                        action='store_true')
    parser.add_argument('--config', help="if print out config file for genome annotation", action='store_true')
    parser.add_argument('--test', type=int, help="only test for small lines, default is None", default=None)
    parser.add_argument('--beddir', type=str, help="concordant and discordant bed file find base dir", default=None)
    parser.add_argument('--dbdir', type=str, help="specify the db file base dir", default=None)
    parser.add_argument('--skip-join-preds', help="assumes the inputs has been there, you can skip this step",
                        action='store_true')
    parser.add_argument('--skip-read-eval', help="if skip read level eval",
                        action='store_true')
    parser.add_argument('--skip-site-eval', help="if skip site level eval",
                        action='store_true')
    parser.add_argument('--region-report', help="if report results at genomic regions",
                        action='store_true')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()

    logger.debug(f"args={args}")

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)

    bed_temp_dir = os.path.join(args.bedtools_tmp, f"{args.dsname}_join_preds")
    os.makedirs(bed_temp_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(bed_temp_dir)

    ## Infer call and cutoff for all tools
    callDict = {}
    for callStr in args.calls:
        toolName, toolFile = callStr.strip().split(':')
        callDict[toolName] = toolFile
    cutoffDict = {}
    for ind, toolName in enumerate(callDict.keys()):
        if args.score_cutoff is not None and ind < len(args.score_cutoff):  # preset cutoff
            cutoffDict[toolName] = args.score_cutoff[ind]
        else:  # default values
            if toolName.lower() == 'nanopolish':
                cutoffDict[toolName] = DEFAULT_CUTOFF_NANOPOLISH
            elif toolName.lower() == 'megalodon':
                cutoffDict[toolName] = DEFAULT_CUTOFF_MEGALODON
            else:
                cutoffDict[toolName] = 0.0
    logger.debug(f"callDict={callDict}, cutoffDict={cutoffDict}")

    global progress_bar_global_join_preds

    if args.skip_join_preds:
        logger.info(
            f"Assume you have generated all joined preds db at:{out_dir} or {args.dbdir}, make sure this is correct if you skip make predictions joined db by --skip-join-preds.")
    else:
        logger.info("Join preds of bs-seq with tool, take times...")
        ## Step 1: load bs-seq, tool read-level unified inputs, make joined preds DF
        ## Load bs-seq
        if args.bs_seq_bed is not None:
            logger.info(f"Start load bs-seq from file:{args.bs_seq_bed}")
            bsseq_df = pd.read_csv(args.bs_seq_bed, sep='\t', header=None, index_col=False).iloc[:, [0, 2, 5, 6, 7]]
            bsseq_df.columns = SITE_COLUMNS + ['Freq', 'Coverage']
            bsseq_df.drop_duplicates(subset=SITE_COLUMNS, inplace=True)
            bsseq_df = bsseq_df[bsseq_df['Coverage'] >= args.min_bgtruth_cov]
            logger.debug(f"bsseq_df={bsseq_df}")
        else:
            bsseq_df = None

        ## join the preds for all tools
        logger.info(f"Start load tool files for tool:{callDict.keys()}")
        executor = ThreadPoolExecutor(max_workers=args.processors)
        all_task_future = []

        progress_bar_global_join_preds = tqdm(total=len(args.chrs))
        progress_bar_global_join_preds.set_description(f"MergePreds-{args.dsname}")

        for chr in args.chrs:
            future = executor.submit(joined_preds_for_all_tools, callDict, chrFilter=[chr],
                                     bsseq_df=bsseq_df, dsname=args.dsname,
                                     test=args.test, is_save=True, chunksize=args.chunksize,
                                     outdir=os.path.join(out_dir, 'joined_db'))
            future.add_done_callback(update_progress_bar_join_preds_eval)
            all_task_future.append(future)
        executor.shutdown()
        progress_bar_global_join_preds.close()
        logger.info(f"Memory report: {get_current_memory_usage()}")
        logger.info(f"### DONE for make joined prediction db.")

    ## Step 2: read level eval on joined preds
    new_dtype = dict(DTYPE)
    for toolName in callDict:
        new_dtype.update({toolName: float})
    logger.debug(f"new_dtype={new_dtype}")

    if not args.skip_read_eval or not args.skip_site_eval:
        if args.region_report:
            logger.debug("Create region bed list firstly, take times......")
            region_bed_list = create_region_bed_list(args.genome_annotation, args.beddir, args.dsname)

    if args.skip_read_eval:
        logger.info("You do not need read level eval, we will skip this section.")
    else:
        read_df = find_join_preds_bgtruth_as_df(args.dsname, args.chrs,
                                                out_dir if args.dbdir is None else args.dbdir, cutoffDict,
                                                dtype=new_dtype, certainSites=True,
                                                fully_meth_threshold=args.fully_meth_threshold,
                                                is_save=True,
                                                outdir=out_dir)
        logger.debug(f"read_df={read_df}")
        report_read_level_performance(read_df, toolList=list(callDict.keys()), dsname=args.dsname, outdir=out_dir,
                                      if_save=True)
        logger.info(f"Memory report: {get_current_memory_usage()}")
        logger.info(f"### DONE for read-level eval for all preds")

        if args.region_report:
            logger.info(f"Start report read-level at genomic regions, take time...")
            read_df1 = read_df.drop('ID', axis=1)
            read_df1['End'] = read_df1['Pos']
            read_df1['Gene'] = '.'
            read_df1['Score'] = '.'
            read_df1 = read_df1[
                ['Chr', 'Pos', 'End', 'Gene', 'Score', 'Strand'] +
                list(callDict.keys()) + ['Freq', 'Coverage']]
            logger.debug(f"read_df1={read_df1}")

            ## preparing multi-threading
            progress_bar_global_join_preds = tqdm(total=len(region_bed_list))
            progress_bar_global_join_preds.set_description(f"Read-level-MT-{args.dsname}")
            all_task_future = []
            executor = ThreadPoolExecutor(max_workers=args.processors)

            ## submit multi-threading function
            for infn, regionName, bedObj in region_bed_list:
                future = executor.submit(
                    eval_read_level_at_region, infn, regionName, read_df1, args.dsname, callDict)
                future.add_done_callback(update_progress_bar_join_preds_eval)
                all_task_future.append(future)
            executor.shutdown()
            progress_bar_global_join_preds.close()

            ## collecting multi-threading results
            dflist = []
            for future in all_task_future:
                if future.result() is not None:
                    dflist.append(future.result())
            region_df = pd.concat(dflist)
            logger.debug(f"Read level region_df={region_df}")

            outfn = os.path.join(out_dir, f"{args.dsname}_read_level_joined_preds_perf_genomic_regions.csv")
            region_df.to_csv(outfn, index=False)
            logger.info(f"save to {outfn}")

            logger.info(f"Memory report: {get_current_memory_usage()}")
            logger.info(f"### DONE for read-level eval for genomic regions")
        read_df = read_df1 = None

    ## Step 3: site level eval on joined preds
    if args.skip_site_eval:
        logger.info("You do not need site level eval, we will skip this section")
    else:
        read_df = find_join_preds_bgtruth_as_df(args.dsname, args.chrs,
                                                out_dir if args.dbdir is None else args.dbdir,
                                                cutoffDict,
                                                dtype=new_dtype,
                                                is_save=True,
                                                outdir=out_dir)
        site_df = read_to_site_df(read_df, list(callDict.keys()), args.dsname, args.toolcov_site, if_save=True,
                                  outdir=out_dir)

        report_site_level_performance(site_df, list(callDict.keys()), dsname=args.dsname, if_save=True, outdir=out_dir)
        logger.info(f"Memory report: {get_current_memory_usage()}")
        logger.info(f"### DONE for site-level eval for all preds")

        if args.region_report:
            logger.info(f"Start report site-level at genomic regions, take time...")
            site_df1 = site_df.copy()
            site_df1['End'] = site_df1['Pos']
            site_df1['Gene'] = '.'
            site_df1['Score'] = '.'
            site_df1 = site_df1[
                ['Chr', 'Pos', 'End', 'Gene', 'Score', 'Strand'] +
                list(callDict.keys()) + ['Freq', 'Reads']]
            logger.debug(f"site_df1={site_df1}")

            ## preparing multi-threading
            progress_bar_global_join_preds = tqdm(total=len(region_bed_list))
            progress_bar_global_join_preds.set_description(f"Site-level-MT-{args.dsname}")
            all_task_future = []
            executor = ThreadPoolExecutor(max_workers=args.processors)

            ## submit multi-threading function
            for infn, regionName, bedObj in region_bed_list:
                future = executor.submit(
                    eval_site_level_at_region, infn, regionName, site_df1, args.dsname, callDict)
                future.add_done_callback(update_progress_bar_join_preds_eval)
                all_task_future.append(future)
            executor.shutdown()
            progress_bar_global_join_preds.close()

            ## collecting multi-threading results
            dflist = []
            for future in all_task_future:
                if future.result() is not None:
                    dflist.append(future.result())
            region_df = pd.concat(dflist)
            logger.debug(f"Site level region_df={region_df}")

            outfn = os.path.join(out_dir, f"{args.dsname}_site_level_joined_preds_perf_genomic_regions.csv")
            region_df.to_csv(outfn, index=False)
            logger.info(f"save to {outfn}")

            logger.info(f"Memory report: {get_current_memory_usage()}")
            logger.info(f"### DONE for site-level eval for genomic regions")

    save_done_file(out_dir)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    logger.info("### Join predictions performance analysis DONE.")
