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

from nanome.common.eval_common import get_dna_base_from_reference
from nanome.common.global_config import *
from nanome.common.global_settings import HUMAN_CHR_SET


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

    if args.cmd == 'sanity-check-seq':
        ## bash meth_stats_tool.sh sanity-check-seq --chrs chr4:10164 chr4:10298
        for chrstr in args.chrs:
            # logger.info(chrstr)
            sanity_check_get_dna_seq(chrstr)
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
