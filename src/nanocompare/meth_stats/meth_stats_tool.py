"""
Tool for pre-processing results

"""
import argparse
import sys

nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
sys.path.append(nanocompare_prj)
from nanocompare.meth_stats.meth_stats_common import load_tombo_df, load_deepmod_df, get_dna_sequence_from_reference, load_sam_as_strand_info_df, load_nanopolish_df

from multiprocessing import Pool

from tqdm import tqdm
import pandas as pd
from global_config import *

import os
import numpy as np

from Bio import SeqIO


def add_strand_info_for_nanopolish(nanopolish_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.methylation_calls.tsv', sam_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.sam'):
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
        raise Exception("We found the read-name of Nanopolish results is not mapped all to SAM/BAM file, please check if the BAM file is used for Nanopolish")

    # df = df.iloc[:, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]

    outfn = os.path.join(pic_base_dir, f'{os.path.splitext(os.path.basename(nanopolish_fn))[0]}-nanopolish-strand-info.tsv')
    df.to_csv(outfn, sep='\t', index=False)
    logger.info(f'save to {outfn}')
    return df


def sanity_check_get_dna_seq(chr='chr1', start_list=[11027, 11068, 11129, 10563, 10484, 10660]):
    for start in start_list:
        ret = get_dna_sequence_from_reference(chr, start, ref_fasta=ref_fasta)
        logger.info(f'chr={chr}, start={start}, ret={ret}')


def filter_noncg_sites_ref_seq(df, tagname, ntask=1, ttask=1, num_seq=5, chr_col=0, start_col=1, strand_col=5, toolname='tombo'):
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

        ret = get_dna_sequence_from_reference(chr, start, num_seq=num_seq, ref_fasta=ref_fasta)
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
    # then outfn is like 'K562.tombo.perReadsStats.combine-with-seq-info-n300-t001-chr1.tsv'
    outfn = os.path.join(args.o, f'{tagname}-with-seq-info-n{ntask}-t{ttask:03d}-{rep_chr}.tsv')
    df.to_csv(outfn, sep='\t', header=False, index=False)
    logger.info(f"save to {outfn}")


def filter_noncg_sites_ref_seq_mpi(df, tagname, ntask=1, ttask=1, num_dna_seq=5, chr_col=0, start_col=1, strand_col=5, toolname='tombo', print_first=False):
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

        ret = get_dna_sequence_from_reference(chr, start, num_seq=num_dna_seq, ref_fasta=ref_fasta)
        seq_col.append(ret)

        if toolname == 'tombo':
            if ret[5:7] == 'CG':
                only_cpg_pattern_index.append(index)
        elif toolname == 'deepmod':
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
    # then outfn is like 'K562.tombo.perReadsStats.combine-with-seq-info-n300-t001-chr1.tsv'
    # outfn = os.path.join(args.o, f'{tagname}-with-seq-info-n{ntask}-t{ttask:03d}-{rep_chr}.tsv')
    # df.to_csv(outfn, sep='\t', header=False, index=False)
    # logger.info(f"save to {outfn}")
    logger.info(f"Finished of subprocess {ttask}:{ntask}")
    return df


def filter_noncg_sites_for_tombo(tombo_fn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.tombo_perReadsStats.bed', sam_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.sorted.bam', ntask=1, ttask=1, num_seq=5):
    if args.i is not None:
        tombo_fn = args.i

    df = load_tombo_df(infn=tombo_fn)
    basefn = os.path.basename(tombo_fn)
    basename = os.path.splitext(basefn)[0]
    filter_noncg_sites_ref_seq(df=df, tagname=basename, ntask=ntask, ttask=ttask, num_seq=num_seq)


def filter_noncg_sites_for_tombo_mpi():
    """
    MPI version of filter out non-CG patterns
    :return:
    """
    ntask = 300

    basefn = os.path.basename(args.i)
    basename = os.path.splitext(basefn)[0]

    df = load_tombo_df(infn=args.i)
    all_list = list(range(len(df)))

    df_list = []
    with Pool(processes=args.processors) as pool:
        for epoch in range(ntask):
            cpg_pattern_index = subset_of_list(all_list, ntask, epoch + 1)
            seldf = df.iloc[cpg_pattern_index, :]
            # logger.info(seldf)

            df_list.append(pool.apply_async(filter_noncg_sites_ref_seq_mpi, (seldf, basename, ntask, epoch + 1)))

            # filter_noncg_sites_ref_seq_mpi(seldf, tagname=basename, ntask=ntask, ttask=epoch + 1)
        pool.close()
        pool.join()

    # Combine df
    logger.debug("Start to combine all results")
    df_list = [df1.get() for df1 in df_list]
    retdf = pd.concat(df_list)
    logger.debug(retdf)

    ## Note: original   input=K562.tombo.perReadsStats.combine.tsv
    ##                  output=K562.tombo.perReadsStatsOnlyCpG.combine.tsv
    basefn = basefn.replace("perReadsStats", "perReadsStatsOnlyCG")
    outfn = os.path.join(args.o, f'{basefn}')
    retdf.to_csv(outfn, sep='\t', index=False, header=False)
    logger.debug(f"Save to {outfn}")


def filter_noncg_sites_for_deepmod(deepmod_fn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.deepmod_combined.bed', sam_fn='/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.sorted.bam', ntask=1, ttask=1, num_seq=5):
    if args.i is not None:
        deepmod_fn = args.i

    df = load_deepmod_df(infn=deepmod_fn)
    basefn = os.path.basename(deepmod_fn)
    basename = os.path.splitext(basefn)[0]
    filter_noncg_sites_ref_seq(df=df, tagname=basename, ntask=ntask, ttask=ttask, num_seq=num_seq, chr_col=0, start_col=1, strand_col=5, toolname='deepmod')


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


def parse_arguments():
    """
    usage: volume_calculation.py [-h] [-n N] [-t T] [--show]
                             [--input INPUT [INPUT ...]] [--output OUTPUT]
                             [--dcm] [--single-scan] [--lu-score LU_SCORE]
                             [--le-score LE_SCORE]
                             cmd

    Volume calculation for lung and lesion

    positional arguments:
      cmd                   name of command: compute, combine, or gen-pixel-info

    optional arguments:
      -h, --help            show this help message and exit
      -n N                  the total number of tasks (1-27)
      -t T                  the current task id (1-N)
      --show                show prediction images if using this switch
      --input INPUT [INPUT ...]
                            the input dir that contains scanid of pic/dcm files
      --output OUTPUT       the input pic dir
      --dcm                 folders are scanid that containing DCM files if using
                            this switch
      --single-scan         folders are directly the scanid folder if using this
                            switch
      --lu-score LU_SCORE   the lung field detection score
      --le-score LE_SCORE   the lesion field detection score
    :return:
    """
    parser = argparse.ArgumentParser(description='Multi-task')
    parser.add_argument("cmd", help="name of command: compute, combine, or gen-pixel-info")
    parser.add_argument('-n', type=int, help="the total number of tasks (1-27)", default=1)
    parser.add_argument('-t', type=int, help="the current task id (1-N)", default=1)
    parser.add_argument('-i', type=str, help="input file", default=None)
    parser.add_argument('-o', type=str, help="output dir or file", default=pic_base_dir)
    parser.add_argument('--ibam', type=str, help="input bam/sam file", default=None)
    parser.add_argument('--processors', type=int, help="Number of processors", default=8)
    parser.add_argument('--mpi', action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    set_log_debug_level()
    args = parse_arguments()
    logger.debug(args)

    ref_fasta = None
    if args.cmd in ['tombo-add-seq', 'deepmod-add-seq', 'sanity-get-seq']:
        ref_fn = '/projects/li-lab/Ziwei/Nanopore/data/reference/hg38.fa'
        ref_fasta = SeqIO.to_dict(SeqIO.parse(open(ref_fn), 'fasta'))

    if args.cmd == 'tombo-add-seq':
        if args.mpi:
            logger.debug('in mpi mode')

            import multiprocessing

            logger.debug("There are %d CPUs on this machine by multiprocessing.cpu_count()" % multiprocessing.cpu_count())

            filter_noncg_sites_for_tombo_mpi()
        else:
            filter_noncg_sites_for_tombo(ntask=args.n, ttask=args.t)
    elif args.cmd == 'deepmod-add-seq':
        filter_noncg_sites_for_deepmod(ntask=args.n, ttask=args.t)
    elif args.cmd == 'nanopolish-add-strand':
        add_strand_info_for_nanopolish()
    elif args.cmd == 'sanity-get-seq':
        sanity_check_get_dna_seq()
    # samfile = pysam.AlignmentFile('/projects/li-lab/yang/results/12-09/K562.nanopolish/K562.sorted.bam', "rb")
    # #
    # chr = 'chr1'
    # start = 45834
    # strand_info = '+'
    # #
    # ret = get_dna_sequence_from_samfile(chr, start, start + 4, samfile)
    # logger.info(f'chr={chr}, start={start}, strand={strand_info}, ret={ret}')
