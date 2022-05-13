#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : meth_stats_tool.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Tool for pre-processing results
"""
import argparse
import glob
import gzip
import sys
from collections import defaultdict
from multiprocessing import Pool

import h5py
import numpy as np
import pandas as pd
from Bio import SeqIO
from ont_fast5_api.fast5_interface import get_fast5_file
from tqdm import tqdm

from nanome.common.eval_common import load_tombo_df, load_deepmod_df, get_dna_base_from_reference, \
    load_sam_as_strand_info_df, load_nanopolish_df
from nanome.common.global_config import *
from nanome.common.global_settings import HUMAN_CHR_SET


def add_strand_info_for_nanopolish(
        nanopolish_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.methylation_calls.tsv',
        sam_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.sam'):
    """
    No need for new nanopolish output
    Combine the nanopolish output tsv results with strand-info from SAM files. This will add last column as strand-info.

    This is due to original nanopolish output results contain no strand-info, we are going to solve this problem.

    Return results columns are:
     [(0, 'chromosome'), (1, 'start'), (2, 'end'), (3, 'read_name'), (4, 'log_lik_ratio'), (5, 'log_lik_methylated'), (6, 'log_lik_unmethylated'), (7, 'num_calling_strands'), (8, 'num_cpgs'), (9, 'sequence'), (10, 'strand-info')]


    :param nanopolish_fn: nanopolish file name
    :param sam_fn: SAM file name for strand-info
    :return:
    """
    if args.i is not None:
        nanopolish_fn = args.i

    if args.ibam is not None:
        sam_fn = args.ibam

    df2 = load_sam_as_strand_info_df(infn=sam_fn)
    df1 = load_nanopolish_df(infn=nanopolish_fn)

    df = df1.merge(df2, left_on='read_name', right_on='read-name', how='left')
    df = df.drop('read-name', axis=1)
    logger.info(df)
    logger.info(list(enumerate(df.columns)))

    if len(df1) != len(df):
        raise Exception(
            "We found the read-name of Nanopolish results is not mapped all to SAM/BAM file, please check if the BAM file is used for Nanopolish")

    # df = df.iloc[:, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]

    outfn = os.path.join(pic_base_dir,
                         f'{os.path.splitext(os.path.basename(nanopolish_fn))[0]}-nanopolish-strand-info.tsv')
    df.to_csv(outfn, sep='\t', index=False)
    logger.info(f'save to {outfn}')
    return df


def sanity_check_get_dna_seq(chrstr):
    """
    Check 0-based start, input as 'chr1:123'
    :param chrstr:
    :return:
    """

    chr, start = chrstr.strip().split(':')
    start = int(start)

    show_arrow = ''.join(['~'] * 5 + ['↑'] + ['~'] * 5)

    ret = get_dna_base_from_reference(chr, start, ref_fasta=ref_fasta)
    logger.info(f'chr={chr}, start={start}\nSEQ={ret}\nPOS={show_arrow}')


def filter_noncg_sites_ref_seq(df, tagname, ntask=1, ttask=1, num_seq=5, chr_col=0, start_col=1, strand_col=5,
                               toolname='tombo'):
    """
    Filter out rows that are non-CG patterns in Tombo results, reference sequence is based on BAM files

    from SAM to BAM (with index) script is as follows:

    samtools view -S -b K562.sam > K562.bam
    samtools sort -o K562.sorted.bam K562.bam
    samtools index K562.sorted.bam

    :param tombo_fn:
    :param sam_fn:
    :return:
    """

    chrs = df.iloc[:, chr_col].unique()
    chrs = np.sort(chrs)
    logger.info(chrs)
    logger.info(len(chrs))

    all_list = list(range(len(df)))
    cpg_pattern_index = subset_of_list(all_list, ntask, ttask)

    # sel_chrs = subset_of_list(chrs, ntask, ttask)
    # logger.info(sel_chrs)
    # df = df[df[0].isin(sel_chrs)]
    df = df.iloc[cpg_pattern_index, :]
    logger.info(df)

    rep_chr = df.iloc[0, chr_col]

    seq_col = []
    cpg_pattern_index = []

    print_first = True
    for index, row in tqdm(df.iterrows()):
        if print_first:
            logger.info(f"index={index}, row={row}")
            print_first = False
        chr = row[chr_col]
        start = int(row[start_col])
        strand_info = row[strand_col]

        # ret = get_dna_sequence_from_samfile(chr, start, start + num_seq, samfile)  # may return None, if no sequence at all reads

        ret = get_dna_base_from_reference(chr, start, num_seq=num_seq, ref_fasta=ref_fasta)
        seq_col.append(ret)

        if toolname == 'tombo':
            if ret[5:7] == 'CG':
                cpg_pattern_index.append(index)
        elif toolname == 'deepmod':
            if strand_info == '+':
                if ret[5:7] == 'CG':
                    cpg_pattern_index.append(index)
            elif strand_info == '-':
                if ret[4:6] == 'CG':
                    cpg_pattern_index.append(index)

    # TODO: using ret if it is CG pattern, or will remove later

    # logger.info(f'chr={chr}, start={start}, strand={strand_info}, ret={ret}')
    # if index > 10000:
    #     break
    df['sequence'] = seq_col

    logger.debug(f'before filter:{len(df)}, after non-CG filter:{len(cpg_pattern_index)}')
    df = df.loc[cpg_pattern_index, :]

    # tagname is like 'K562.tombo.perReadsStats.combine'
    # then outfn is like 'K562.tombo.perReadsStats.combine-with-seq-info-t001-chr1.tsv'
    outfn = os.path.join(args.o, f'{tagname}-with-seq-info-n{ntask}-t{ttask:03d}-{rep_chr}.tsv')
    df.to_csv(outfn, sep='\t', header=False, index=False)
    logger.info(f"save to {outfn}")


def filter_noncg_sites_ref_seq_mpi(df, tagname, ntask=1, ttask=1, num_dna_seq=5, chr_col=0, start_col=1, strand_col=5,
                                   toolname='tombo', print_first=False):
    """
    MPI version
    invoke like: res = p.apply_async(testFunc, args=(2, 4), kwds={'calcY': False})
                    or pool.apply_async(test, (t,), dict(arg2=5))

    Filter out rows that are non-CG patterns in Tombo results, reference sequence is based on BAM files

    :param tombo_fn:
    :param sam_fn:
    :return:
    """

    rep_chr = df.iloc[0, chr_col]

    seq_col = []
    only_cpg_pattern_index = []

    for index, row in df.iterrows():
        if print_first:
            logger.info(f"index={index}, row={row}")
            print_first = False
        chr = row[chr_col]
        start = int(row[start_col])
        strand_info = row[strand_col]

        ret = get_dna_base_from_reference(chr, start, num_seq=num_dna_seq, ref_fasta=ref_fasta)
        seq_col.append(ret)

        if toolname == 'tombo':
            if ret[5:7] == 'CG':
                only_cpg_pattern_index.append(index)
        elif toolname in ['deepmod', 'deepmod-read-level']:
            if strand_info == '+':
                if ret[5:7] == 'CG':
                    only_cpg_pattern_index.append(index)
            elif strand_info == '-':
                if ret[4:6] == 'CG':
                    only_cpg_pattern_index.append(index)

    df['sequence'] = seq_col

    # logger.debug(f'Subprocess [{ttask}:{ntask}] finished, before filter:{len(df)}, after non-CG filter:{len(only_cpg_pattern_index)}')
    df = df.loc[only_cpg_pattern_index, :]

    # tagname is like 'K562.tombo.perReadsStats.combine'
    # then outfn is like 'K562.tombo.perReadsStats.combine-with-seq-info-t001-chr1.tsv'
    # outfn = os.path.join(args.o, f'{tagname}-with-seq-info-n{ntask}-t{ttask:03d}-{rep_chr}.tsv')
    # df.to_csv(outfn, sep='\t', header=False, index=False)
    # logger.info(f"save to {outfn}")
    logger.info(f"Finished of subprocess {ttask}:{ntask}")
    return df


def filter_noncg_sites_for_tombo(
        tombo_fn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.tombo_perReadsStats.bed',
        sam_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.sorted.bam', ntask=1, ttask=1, num_seq=5):
    if args.i is not None:
        tombo_fn = args.i

    df = load_tombo_df(infn=tombo_fn)
    basefn = os.path.basename(tombo_fn)
    basename = os.path.splitext(basefn)[0]
    filter_noncg_sites_ref_seq(df=df, tagname=basename, ntask=ntask, ttask=ttask, num_seq=num_seq)


def convert_bismark_add_strand_and_seq(indf, outfn):
    """
    Check start pointer, if point to CG's C, it is positive strand, or else, it is reverse strand
    Note: input file is 1-based start, we also output to a 1-based format that is compatable to our Bismark import functions.
    :param indf:
    :param outf:
    :param report_num:
    :return:
    """
    logger.debug(f'Start add strand and seq to bismark cov file, total len={len(indf)}')

    outf = gzip.open(outfn, 'wt')

    for index, row in tqdm(indf.iterrows(), total=len(indf), desc='Bismark_cov'):
        # if report_num and index % report_num == 0:
        #     logger.debug(f'processed index={index}')
        chr = row['chr']
        start = int(row['start'])  # Keep raw 1-based format of bismark results
        ret = get_dna_base_from_reference(chr, start - 1, ref_fasta=ref_fasta)
        if ret[5] == 'C':  # strand is +
            strand = '+'
        elif ret[5] == 'G':
            strand = '-'
        else:
            raise Exception(f'We can not identify this bg-truth file with non-CG results, such as row={row}')

        outstr = '\t'.join([chr, str(start), strand, str(row['mcount']), str(row['ccount']), ret[4:7]])
        outf.write(f'{outstr}\n')
    outf.close()
    logger.info(f'save to {outfn}')

    logger.debug(f'Finish add strand info task')


def convert_bismark_cov_to_gw_format(df, tagname=None):
    """
    Save adding strand info and dna seq format, which is in same format of Bismark Genome-wide output files

    Input format:
    chr17	61127	61127	100	1	0
    chr17	61139	61139	100	1	0
    chr17	61193	61193	100	1	0

    Output format:
    chr17	61127	+	1	0	GCG
    chr17	61139	+	1	0	ACG
    chr17	61193	+	1	0	CCG

    :param df:
    :return:
    """
    basefn = os.path.basename(args.i)
    basename = os.path.splitext(basefn)[0]

    outfn = os.path.join(args.o,
                         f'{basename.replace(".gz", "")}.{"" if tagname is None else tagname}.convert.add.strand.tsv.gz')
    convert_bismark_add_strand_and_seq(df, outfn)


def filter_noncg_sites_mpi(df, ntask=300, toolname='tombo'):
    """
    MPI version of filter out non-CG patterns
    :return:
    """
    basefn = os.path.basename(args.i)
    basename = os.path.splitext(basefn)[0]

    all_list = list(range(len(df)))

    # Store each sub-process return results
    df_list = []
    with Pool(processes=args.processors) as pool:
        for epoch in range(ntask):
            cpg_pattern_index = subset_of_list(all_list, ntask, epoch + 1)
            seldf = df.iloc[cpg_pattern_index, :]

            if toolname == 'tombo':
                df_list.append(pool.apply_async(filter_noncg_sites_ref_seq_mpi, (seldf, basename, ntask, epoch + 1)))
            elif toolname == 'deepmod':
                df_list.append(pool.apply_async(filter_noncg_sites_ref_seq_mpi, (seldf, basename, ntask, epoch + 1),
                                                dict(chr_col=0, start_col=1, strand_col=5, toolname='deepmod')))
            elif toolname == 'deepmod-read-level':
                df_list.append(pool.apply_async(filter_noncg_sites_ref_seq_mpi, (seldf, basename, ntask, epoch + 1),
                                                dict(chr_col=0, start_col=1, strand_col=5,
                                                     toolname='deepmod-read-level')))
            else:
                raise Exception(f"{toolname} is no valid.")
        pool.close()
        pool.join()

    # Combine df
    logger.debug("Start to combine all results")
    df_list = [df1.get() for df1 in df_list]
    retdf = pd.concat(df_list)
    logger.debug(retdf)

    ## Note: original   input=K562.tombo.perReadsStats.combine.tsv
    ##                  output=K562.tombo.perReadsStatsOnlyCpG.combine.tsv

    if toolname == 'tombo':
        basefn = basefn.replace("perReadsStats", "perReadsStatsOnlyCG").replace("combined", "combine")
    elif toolname == 'deepmod':
        ## Example: HL60.deepmod.C.combined.tsv
        basefn = basefn.replace(".C.", ".C_OnlyCG.").replace("combined", "combine")
    else:
        raise Exception(f"{toolname} is no valid.")

    outfn = os.path.join(args.o, f'{basefn}')
    retdf.to_csv(outfn, sep='\t', index=False, header=False)
    logger.debug(f"Save to {outfn}")


def filter_noncg_sites_for_deepmod(
        deepmod_fn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.deepmod_combined.bed',
        sam_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.sorted.bam', ntask=1, ttask=1, num_seq=5):
    if args.i is not None:
        deepmod_fn = args.i

    df = load_deepmod_df(infn=deepmod_fn)
    basefn = os.path.basename(deepmod_fn)
    basename = os.path.splitext(basefn)[0]
    filter_noncg_sites_ref_seq(df=df, tagname=basename, ntask=ntask, ttask=ttask, num_seq=num_seq, chr_col=0,
                               start_col=1, strand_col=5, toolname='deepmod')


def subset_of_list(alist, n, t):
    """
    Subset of a list for multi-processing
        n=1 to 100
        t=1 to N
        return subset list of alist
    :param alist:
    :param n:
    :param t:
    :return:
    """
    if t < 1 or t > n:
        raise Exception(f't={t} is not accept, must be 1-N (include)')

    if n > len(alist):  # if n is bigger than all list, return only 1 for t<=len
        if t <= len(alist):
            return [alist[t - 1]]
        else:
            return None

    m = int(len(alist) / n)  # each task of a section of list

    start_index = int((t - 1) * m)
    if t == n:
        sublist = alist[start_index:]
    else:
        sublist = alist[start_index:start_index + m]
    # logger.debug(f'n={n}, t={t}, section={m}, index={start_index}:{start_index + m}')
    return sublist


def get_f5_readid_map(flist):
    f5_readid_map = defaultdict(str)
    for fn in flist:
        basename = os.path.basename(fn)
        with get_fast5_file(fn, mode="r") as f5:
            # if len(f5.get_reads()) >= 2:
            #     raise Exception(f'We can only deal with one read in fast5, but fn={fn}, contains {len(f5.get_reads())} multiple reads')
            for read in f5.get_reads():
                f5_readid_map[basename] = str(read.read_id)
    return f5_readid_map


def build_map_fast5_to_readid_mp(
        basedir='/fastscratch/liuya/nanocompare/K562-Runs/K562-DeepMod-N50/K562-DeepMod-N50-basecall', ntask=300):
    patfn = os.path.join(basedir, '**', '*.fast5')
    fast5_flist = glob.glob(patfn, recursive=True)

    logger.info(f'Total fast5 files: {len(fast5_flist)}')

    ret_list = []
    with Pool(processes=args.processors) as pool:
        for epoch in range(ntask):
            subflist = subset_of_list(fast5_flist, ntask, epoch + 1)
            ret_list.append(pool.apply_async(get_f5_readid_map, (subflist,)))
        pool.close()
        pool.join()
    logger.debug('Finish fast5 to read-id mapping')
    f5_readid_map = defaultdict(str)
    for ret in ret_list:
        f5_readid_map.update(ret.get())

    # for fn in fast5_flist[:]:
    #     # logger.debug(fn)
    #     basename = os.path.basename(fn)
    #
    #     with get_fast5_file(fn, mode="r") as f5:
    #         for read in f5.get_reads():
    #             # logger.debug(read.read_id)
    #             f5_readid_map[basename] = str(read.read_id)
    return f5_readid_map


def process_pred_detail_f5file(fn, f5_readid_map):
    """
    For each deepmod prediction results file, we analyze a df result of read-level results
    :param fn:
    :param f5_readid_map:
    :return:
    """

    f5_pred_key = '/pred/pred_0/predetail'
    dflist = []
    with h5py.File(fn, 'r') as mr:
        # m_pred = mr[f5_pred_key].value
        # logger.debug(m_pred)
        for name in mr['/pred']:
            # logger.debug(name)
            pred_num_key = f'/pred/{name}'
            f5file = os.path.basename(mr[pred_num_key].attrs['f5file'])
            mapped_chr = mr[pred_num_key].attrs['mapped_chr']
            mapped_strand = mr[pred_num_key].attrs['mapped_strand']

            # logger.debug(f'{pred_num_key}: chr={mapped_chr}, strand={mapped_strand}, f5file={f5file}')

            pred_detail_key = f'{pred_num_key}/predetail'
            # m_pred = mr[pred_detail_key].value
            m_pred = mr[pred_detail_key][()]
            m_pred = np.array(m_pred, dtype=[('refbase', 'U1'), ('readbase', 'U1'), ('refbasei', np.uint64),
                                             ('readbasei', np.uint64), ('mod_pred', np.int)])

            dataset = []
            for mi in range(len(m_pred)):
                if m_pred['refbase'][mi] not in ['C']:
                    continue
                if m_pred['refbase'][mi] in ['-', 'N', 'n']:
                    continue
                # if m_pred['readbase'][mi] == '-':
                #     continue

                # Filter non-CG patterns results
                ret = get_dna_base_from_reference(mapped_chr, int(m_pred['refbasei'][mi]), ref_fasta=ref_fasta)

                if mapped_strand == '+':
                    if ret[5:7] != 'CG':
                        continue
                elif mapped_strand == '-':
                    if ret[4:6] != 'CG':
                        continue

                if -0.1 < m_pred['mod_pred'][mi] - 1 < 0.1:
                    meth_indicator = 1
                else:
                    meth_indicator = 0
                # sp_options['4NA'][m_pred['refbase'][mi]][(cur_chr, cur_strand, int(m_pred['refbasei'][mi]) )][0] += 1
                ret = {'start': int(m_pred['refbasei'][mi]), 'pred': meth_indicator, 'base': m_pred['refbase'][mi],
                       'sequence': ret}
                dataset.append(ret)
            df = pd.DataFrame(dataset)

            if len(df) < 1:
                continue
            df['chr'] = str(mapped_chr)
            df['end'] = df['start'] + 1
            df['strand'] = str(mapped_strand)
            df['read-id'] = f5_readid_map[f5file]
            df = df[['chr', 'start', 'end', 'read-id', 'base', 'strand', 'sequence', 'pred']]
            # logger.info(df)
            dflist.append(df)

    sumdf = pd.concat(dflist)

    # logger.debug(f'Process pred detail file {fn} finished, total reads={len(sumdf)}.')
    return sumdf


def extract_deepmod_read_level_results_mp(
        basecallDir='/fastscratch/liuya/nanocompare/K562-Runs/K562-DeepMod-N50/K562-DeepMod-N50-basecall',
        methcallDir='/fastscratch/liuya/nanocompare/K562-Runs/K562-DeepMod-N50/K562-DeepMod-N50-methcall', ntask=50):
    f5_readid_map = build_map_fast5_to_readid_mp(basedir=basecallDir, ntask=ntask)
    # logger.debug(f5_readid_map)

    # dirname = '/fastscratch/liuya/nanocompare/K562-Runs/K562-DeepMod-N50/K562-DeepMod-N50-methcall/**/rnn.pred.detail.fast5.*'
    dirname = os.path.join(methcallDir, '**', 'rnn.pred.detail.fast5.*')
    fast5_flist = glob.glob(dirname, recursive=True)
    logger.info(f'Total deepmod fast5 files:{len(fast5_flist)}')

    dflist = []
    with Pool(processes=args.processors) as pool:
        for fn in fast5_flist[:]:
            # df = process_pred_detail_f5file(fn, f5_readid_map)
            dflist.append(pool.apply_async(process_pred_detail_f5file, (fn, f5_readid_map,)))
            # logger.debug(df)
            # logger.debug(df.iloc[1, :])
            # logger.debug(fn)
            # pass
        pool.close()
        pool.join()

    dflist = [df.get() for df in dflist]
    sumdf = pd.concat(dflist)
    logger.debug('Finish get df from deepmod fast5 predetail files')

    cpgDict = defaultdict(lambda: [0, 0])  # 0:cov, 1:meth-cov
    for index, row in sumdf.iterrows():
        chr = row['chr']
        start = row['start']
        strand = row['strand']
        basekey = (chr, start, strand)
        cpgDict[basekey][0] += 1
        if row['pred'] == 1:
            cpgDict[basekey][1] += 1
    logger.debug(f'CpG sites={len(cpgDict)}')

    dataset = []
    for site in cpgDict:
        ret = {'chr': site[0], 'start': site[1], 'end': site[1] + 1, 'base': 'C', 'cap-cov': cpgDict[site][0],
               'strand': site[2], 'no-use1': '', 'start1': site[1], 'end1': site[1] + 1, 'no-use2': '0,0,0',
               'cov': cpgDict[site][0], 'meth-freq': int(100 * cpgDict[site][1] / cpgDict[site][0]),
               'meth-cov': cpgDict[site][1]}
        dataset.append(ret)
    beddf = pd.DataFrame(dataset)
    beddf = beddf[
        ['chr', 'start', 'end', 'base', 'cap-cov', 'strand', 'no-use1', 'start1', 'end1', 'no-use2', 'cov', 'meth-freq',
         'meth-cov']]
    logger.debug('Finish bed df, extract all DONE.')

    return sumdf, beddf


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Multi-task')
    parser.add_argument("cmd", help="name of command: compute, combine, or gen-pixel-info")
    parser.add_argument('-n', type=int, help="the total number of tasks (1-27)", default=1)
    parser.add_argument('-t', type=int, help="the current task id (1-N)", default=1)
    parser.add_argument('-i', type=str, help="input file", default=None)
    parser.add_argument('-o', type=str, help="output dir or file", default=pic_base_dir)
    parser.add_argument('--o2', type=str, help="second output dir or file", default=None)
    parser.add_argument('--ibam', type=str, help="input bam/sam file", default=None)
    parser.add_argument('--basecallDir', type=str, help="basecallDir dir name", default=None)
    parser.add_argument('--methcallDir', type=str, help="methcallDir dir name", default=None)
    parser.add_argument('--processors', type=int, help="Number of processors", default=8)
    parser.add_argument('--mpi', action='store_true')
    parser.add_argument('--chrs', nargs='+', help='all chrs need to check', default=[])

    return parser.parse_args()


def output_bed_by_bin(bin_id):
    num_bins = 5
    density_col = 4
    output_cols = [0, 1, 2]
    bin_value = int(bin_id / num_bins * 100 + 1e-5)

    logger.info(f"start with bin_id={bin_id}, bin_value={bin_value}")

    ndf = df[df[density_col] == bin_value]
    ndf = ndf.iloc[:, output_cols]

    logger.info(f"start to save, df={len(df):,}, ndf={len(ndf):,}, for bin_value={bin_value}")
    outfn = os.path.join(args.o, f"hg38.gc5Base.bin{bin_value}.bed.gz")
    ndf.to_csv(outfn, sep='\t', header=False, index=False)
    logger.info(f"save to {outfn}")


def output_bed_by_bin2(infn, num_bins):
    inf = gzip.open(infn, 'rt')
    outf_list = []
    for bin_id in range(0, num_bins + 1):
        bin_value = int(bin_id / num_bins * 100 + 1e-5)
        outf_list.append(gzip.open(os.path.join(args.o, f"hg38.gc5Base.bin{bin_value}.bed.gz"), 'wt'))

    for row in tqdm(inf):
        tmp = row.strip().split("\t")
        density_col = 4
        bin_value = int(float(tmp[density_col]) + 1e-5)
        bin_id = bin_value // 20
        if bin_id not in range(0, num_bins + 1):
            logger.error(f"Error found: bin_value={bin_value}, bin_id={bin_id}, for row={row}")
            raise Exception(f"Error found: bin_value={bin_value}, bin_id={bin_id}, for row={row}")
        outf_list[bin_id].write(f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}\n")

    [outf.close for outf in outf_list]
    logger.info("Finished bin bed for gc density")


def save_tss_bed_for_5hmc(infn, outfn):
    logger.info(f"open infn={infn}")
    df = pd.read_csv(infn, sep='\t', header=None)
    logger.debug(df)

    df = df.iloc[:, [0, 1, 2, 4, 7]]
    df.columns = ['chr', 'start', 'end', '5hmc_level', 'strand']
    df['n1'] = '.'
    df['start'] = df['start'].astype(int) - 1
    df['end'] = df['end'].astype(int) - 1
    df['5hmc_level'] = df['5hmc_level'].astype(float)
    df = df[['chr', 'start', 'end', '5hmc_level', 'n1', 'strand']]

    logger.info(f"df['5hmc_level'] = {df['5hmc_level'].describe()}")
    logger.info(f"len(df['5hmc_level'] >= 1.0) = {(df.loc[:, '5hmc_level'] >= 1.0 - 1e-3).sum()}")

    df.to_csv(outfn, sep='\t', header=False, index=False)
    logger.info(f"save to {outfn}")
    pass


if __name__ == '__main__':
    set_log_debug_level()
    args = parse_arguments()
    logger.debug(args)

    ref_fasta = None
    if args.cmd in ['tombo-add-seq', 'deepmod-add-seq', 'deepmod-read-level', 'sanity-check-seq',
                    'bismark-convert']:  # These command will use reference genome
        ref_fn = '/projects/li-lab/reference/hg38/hg38.fasta'
        ref_fasta = SeqIO.to_dict(SeqIO.parse(open(ref_fn), 'fasta'))

    if args.cmd == 'tombo-add-seq':
        if args.mpi:
            logger.debug('in mpi mode')

            import multiprocessing

            logger.debug(
                "There are %d CPUs on this machine by multiprocessing.cpu_count()" % multiprocessing.cpu_count())

            df = load_tombo_df(infn=args.i)

            filter_noncg_sites_mpi(df)
        else:
            filter_noncg_sites_for_tombo(ntask=args.n, ttask=args.t)
    elif args.cmd == 'deepmod-add-seq':
        if args.mpi:
            logger.debug('in mpi mode')
            import multiprocessing

            logger.debug(
                "There are %d CPUs on this machine by multiprocessing.cpu_count()" % multiprocessing.cpu_count())

            df = load_deepmod_df(infn=args.i)
            filter_noncg_sites_mpi(df, toolname='deepmod')
        else:
            filter_noncg_sites_for_deepmod(ntask=args.n, ttask=args.t)
    elif args.cmd == 'nanopolish-add-strand':
        add_strand_info_for_nanopolish()
    elif args.cmd == 'sanity-check-seq':
        ## bash meth_stats_tool.sh sanity-check-seq --chrs chr4:10164 chr4:10298
        for chrstr in args.chrs:
            # logger.info(chrstr)
            sanity_check_get_dna_seq(chrstr)
    elif args.cmd == 'deepmod-read-level':
        ### Running bash:
        """
         sbatch meth_stats_tool_mpi.sh deepmod-read-level --basecallDir /fastscratch/liuya/nanocompare/K562-Runs/K562-DeepMod-N50/K562-DeepMod-N50-basecall --methcallDir /fastscratch/liuya/nanocompare/K562-Runs/K562-DeepMod-N50/K562-DeepMod-N50-methcall -o /fastscratch/liuya/nanocompare/deepmod-read-level1.tsv --o2 /fastscratch/liuya/nanocompare/deepmod-read-level1-extract-output.bed
        """

        sumdf, beddf = extract_deepmod_read_level_results_mp(basecallDir=args.basecallDir, methcallDir=args.methcallDir)
        logger.info(sumdf)
        logger.info(sumdf.iloc[1, :])
        logger.info(sumdf['chr'].unique())
        # outfn = os.path.join('/fastscratch/liuya/nanocompare/', 'deepmod-read-level.tsv')

        # Save read level results
        outfn = args.o
        sumdf.to_csv(outfn, sep='\t', index=False, header=False)
        logger.info(f'save to {outfn}')

        if args.o2:  # Save CpG base level results bed file for cluster module use
            outfn = args.o2
            beddf.to_csv(outfn, sep=' ', index=False, header=False)
            logger.info(f'save to {outfn}')
    elif args.cmd == 'bismark-convert':  # Convert non-strand info bismark to strand
        ## bash meth_stats_tool.sh bismark-convert -i /pod/2/li-lab/Ziwei/Nanopore_methyl_compare/result/APL_BSseq/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bed

        ## sbatch meth_stats_tool.sh bismark-convert -i /pod/2/li-lab/Ziwei/Nanopore_methyl_compare/result/APL_BSseq/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bed

        df = pd.read_csv(args.i, sep='\t', header=None, index_col=None)
        if len(df.columns) != 6:
            raise Exception(f"Can no recognize input file format for infn={args.i}, df={df}")
        df.columns = ['chr', 'start', 'end', 'freq100', 'mcount', 'ccount']

        if args.chrs is not None and len(args.chrs) >= 1:
            logger.debug(f"Filter by chrs={args.chrs}")
            df = df[df.chr.isin(args.chrs)].copy()
            df.reset_index(drop=True, inplace=True)
        logger.debug(df)
        convert_bismark_cov_to_gw_format(df, tagname='_'.join(args.chrs))
    elif args.cmd == 'gc-density-bed':
        # sbatch meth_stats_tool.sh gc-density-bed
        infn = "/projects/li-lab/yang/workspace/nano-compare/data/genome-annotation/hg38.gc5Base.bed.gz"
        output_bed_by_bin2(infn, num_bins=5)
        if True:
            sys.exit(0)

        df = pd.read_csv(infn, sep='\t', header=None)
        df.iloc[:, 4] = df.iloc[:, 4].astype(int)
        logger.debug(df)
        bin_list = list(range(1, 6))
        os.makedirs(args.o, exist_ok=True)
        with Pool(processes=args.processors) as pool:
            pool.map(output_bed_by_bin, bin_list)
    elif args.cmd == 'repetitive-bed':
        # sbatch meth_stats_tool.sh repetitive-bed
        # bash meth_stats_tool.sh repetitive-bed
        infn = "/projects/li-lab/yang/results/2021-07-01/hg38.repetitive.bed.gz"
        df = pd.read_csv(infn, sep='\t')
        df = df[df['genoName'].isin(HUMAN_CHR_SET)]
        df['n1'] = '.'
        df['n2'] = '.'
        logger.info(df)
        outfn = f"hg38.repetitive.rep_All.bed.gz"
        df[['genoName', 'genoStart', 'genoEnd', 'n1', 'n2', 'strand']].to_csv(os.path.join(args.o, outfn), sep='\t',
                                                                              header=False, index=False)

        region_dict = {
            "LINE": ["LINE"],
            "SINE": ["SINE"],
            "LTR": ["LTR"],
            "DNA": ["DNA"]
        }
        used_list = []
        for key in region_dict:
            logger.info(f"seperate {key}")
            used_list += region_dict[key]
            ndf = df[df['repClass'].isin(region_dict[key])]
            ndf = ndf[['genoName', 'genoStart', 'genoEnd', 'n1', 'n2', 'strand']]
            # logger.info(ndf)
            outfn = f"hg38.repetitive.rep_{key}.bed.gz"
            ndf.to_csv(os.path.join(args.o, outfn), sep='\t', header=False, index=False)
            logger.info(f"len={len(ndf)}, save to {outfn}")

        ## Output others
        ndf = df[~df['repClass'].isin(used_list)]
        ndf = ndf[['genoName', 'genoStart', 'genoEnd', 'n1', 'n2', 'strand']]
        # logger.info(ndf)
        outfn = f"hg38.repetitive.rep_Others.bed.gz"
        ndf.to_csv(os.path.join(args.o, outfn), sep='\t', header=False, index=False)
        logger.info(f"len={len(ndf)}, save to {outfn}")
    elif args.cmd == 'apl-5hmc-bed':
        # Extract TSS format BED file for 5hmC
        # convert 1-based to 0-based results, output 5hmc level
        # bash meth_stats_tool.sh apl-5hmc-bed
        # file will be later converted into BW file
        infn = "/pod/2/li-lab/Nanopore_compare/data/APL_5hmC_BSseq/APL.cov5.mlml.addstrand.selected.bed.gz"
        outfn = os.path.join(args.o, "APL.5hmc.tss.cov5.bed.gz")
        save_tss_bed_for_5hmc(infn, outfn)

        infn = "/pod/2/li-lab/Nanopore_compare/data/APL_5hmC_BSseq/APL.mlml.addstrand.selected.bed.gz"
        outfn = os.path.join(args.o, "APL.5hmc.tss.cov1.bed.gz")
        save_tss_bed_for_5hmc(infn, outfn)
        pass
    elif args.cmd == 'merge-basecall-summary':
        ## sbatch meth_stats_tool.sh merge-basecall-summary -i /projects/li-lab/yang/results/2021-07-17/NA12878_basecall_logs_output
        baseDir = args.i
        flist = glob.glob(os.path.join(baseDir, '**', '*sequencing_summary.txt'))
        logger.info(flist)
        logger.info(len(flist))
        dflist = []
        for fn in flist:
            df = pd.read_csv(fn, sep='\t')
            dflist.append(df)
        dfall = pd.concat(dflist)
        outfn = os.path.join(args.o, 'NA12878-allChrs-basecall.sequencing_summary.txt')
        dfall.to_csv(outfn, sep='\t', index=False)
        logger.info(f"save to {outfn}")
    else:
        raise Exception(f"Not support command={args.cmd}")

    logger.info("meth_stats_tool DONE")
