#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : computeRawReadsCoverage.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Calculate nanopore raw fast5 reads coverage at each CpG sites based on genomic context for Table S1.
"""
import argparse
import glob
import gzip
import os
import re
from multiprocessing import Pool

import pandas as pd
import pybedtools
from pybedtools import BedTool
from tqdm import tqdm

from nanome.common.eval_common import get_dna_seq_from_reference, open_file_gz_or_txt, find_bed_filename, get_ref_fasta, \
    get_region_bed, intersect_bed_regions, get_region_tagname
from nanome.common.global_config import set_log_debug_level, logger, pic_base_dir, global_temp_dir, set_log_info_level
from nanome.common.global_settings import HUMAN_CHR_SET, datasets_order, reference_genome_hg38_fn, NANOME_VERSION, \
    region_filename_dict

# used for convert region bed cov to base level cov
rawReadDir = '/pod/2/li-lab/Nanopore_compare/data/Nanopore_cov'
rawReadDir = '/projects/li-lab/Nanopore_compare/nanome_paper_result/basecall_results/raw_fast5_reads_cpg_coverage'

# base level cov bed files (combined + - strand)
base_cov_dir = '/projects/li-lab/Nanopore_compare/nanome_paper_result/basecall_results/raw_fast5_reads_cpg_coverage'

# dir to newly concordant and discordant perf results bed dir
bedDir = '/projects/li-lab/yang/results/2021-07-08'


def convert_region_to_cpg_base(dsname):
    """
    Conver bed file of regions into CpG base bed files

    Assume infn is 0-based file, we print 1-based output. Because 1-based site level reports are easy to stats overlap with region files.
    :param infn:
    :return:
    """

    fnlist = glob.glob(os.path.join(rawReadDir, f'{dsname}*coverage.*strand.bed.gz'))
    logger.info(f'convert_region_to_cpg_base:{fnlist}')

    outfn = os.path.join(outdir, f'{dsname}.rawfast5.coverage.base.bed.gz')
    outfile = gzip.open(outfn, 'wt')

    print_first = True
    for infn in fnlist:
        logger.info(f'processing file={infn}')
        infile = open_file_gz_or_txt(infn)
        for row in tqdm(infile):
            tmp = row.strip().split(" ")
            if print_first:
                print_first = False
                logger.info(f'row={row}, tmp={tmp}')
            chr = tmp[0]
            start = int(tmp[1])
            end = int(tmp[2])
            cov = int(tmp[3])
            # TODO: get strand info, + strand link to CG's C, - strand link to CG's G
            strand = tmp[4]

            if chr not in HUMAN_CHR_SET:  # filter out non-human chrs
                continue

            # we want get seq at least two bases, evaluate 'CG' patterns
            # if end == start + 1:
            #     newend = end + 1
            # else:
            #     newend = end
            newend = end + 1

            dnastr = get_dna_seq_from_reference(chr, start, newend, ref_fasta=refFasta)
            for cpg in re.finditer("CG", dnastr):
                cpgstart = start + cpg.start()
                if strand == '+':  # point to CG's C, just report that position (0-based) into 1-based
                    outstart = cpgstart + 1
                elif strand == '-':  # negative strand, point to CG's G, based on seq of + strand
                    outstart = cpgstart + 2
                else:
                    raise Exception(f'strand={strand} is no suppport')

                out_txt = '\t'.join([chr, str(outstart), str(outstart), str(cov), '.', strand])
                outfile.write(f'{out_txt}\n')
                pass
        infile.close()
    outfile.close()


def preprocess_bed_to_cgp_base():
    dsname_list = ['HL60', 'K562', 'APL', 'NA19240']
    dsname_list = ['NA12878']

    for dsname in dsname_list:
        logger.info(f'dsname={dsname}')
        convert_region_to_cpg_base(dsname)


def count_sites_in_coord(readBed, coordfn, tagname, cutoff_list):
    """
    Counts how many CpG sites from readBed are also in range of coordinate file

    Sample readBed files:
    head HL60.rawfast5.coverage.base.bed
    chr1	10470	10470	8	.	-
    chr1	10472	10472	6	.	-
    chr1	10485	10485	4	.	-
    chr1	10490	10490	6	.	-

    :param readBed:
    :param coordfn:
    :param cutoff_list:
    :return:
    """
    ret = {}

    coordBed = get_region_bed(coordfn, enable_base_detection_bedfile=not args.disable_bed_check)
    if not coordBed:
        logger.debug(f"ERROR: not found BED region: tagname={tagname}, coordfn={coordfn}")
        return None

    intersectBed = intersect_bed_regions(readBed, coordBed, coordfn)

    covList = []
    for row in intersectBed:
        cov = int(row[3])
        covList.append(cov)
    covSeries = pd.Series(covList)
    for i, cutoff in enumerate(cutoff_list):
        ret.update({f'{tagname}.cutoff{cutoff}': (covSeries >= cutoff).sum()})

    return ret


def get_len_of_bedtool(bed_obj, cov_col=3, cutoff=1):
    """
    Get the number of rows (sites) of BedTool object
    :param bed_obj:
    :param cutoff:
    :return:
    """
    # covList = []
    # for row in bed_obj:
    #     covList.append(int(row[cov_col]))
    # covSeries = pd.Series(covList)
    df = bed_obj.to_dataframe()
    return (df.iloc[:, cov_col] >= cutoff).sum()


def save_dataset(dataset, tagname=""):
    ## dataframe for raw reads at each genomic regions
    df = pd.DataFrame(dataset)
    logger.info(df)

    column_list = ['dsname', 'total', 'total_cutoff'] + \
                  ["Singletons", "Non-singletons", "Promoters", "Exons", "Introns", "Intergenic", "CpG islands",
                   "CpG shores", "CpG shelves"] + \
                  ["CG_20", "CG_40", "CG_60", "CG_80", "CG_100"] + \
                  ["rep_SINE", "rep_LINE", "rep_LTR", "rep_DNA", "rep_Others"] + \
                  ['Concordant', 'Discordant', 'filename']

    for cutoff in cutoff_list:
        outdf = pd.concat([df.loc[:, ['dsname', 'filename', 'total']], df.filter(regex=f'.cutoff{cutoff}$', axis=1)],
                          axis=1)

        oldname_list = outdf.columns
        newname_list = [the_name.replace(f'.cutoff{cutoff}', "") for the_name in oldname_list]
        outdf.columns = newname_list

        ## Reorder columns based on list, if not exist, skit it
        selcols = [colname for colname in column_list if colname in outdf.columns]
        outdf = outdf[selcols]

        outdf['dsname'] = pd.Categorical(df['dsname'], datasets_order)
        outdf = outdf.sort_values(by='dsname', ascending=True)
        outfn = os.path.join(outdir,
                             f'raw.fast5.reads.cpg.coverage.across.regions.cutoff{cutoff}{tagname}.table.s1.xlsx')
        outdf.to_excel(outfn)
        logger.info(f'save to {outfn}')


def combine_na12878_coverage_bed():
    baseDir = "/fastscratch/liuya/nanocompare/NA12878-coverage"
    outfn = os.path.join(outdir, "NA12878-allChrs.coverage.bothstrand.bed.gz")
    outf = gzip.open(outfn, 'wt')
    for chrName in HUMAN_CHR_SET:
        logger.info(f"Processing chr={chrName}")
        flist = glob.glob(os.path.join(baseDir, f"NA12878-{chrName.upper()}.coverage.*.bed.gz"))
        logger.info(flist)

        for fn in flist:
            logger.info(fn)
            with gzip.open(fn, 'rt') as inf:
                for row in tqdm(inf):
                    tmp = row.strip().split(" ")
                    chr = tmp[0]
                    if chr != chrName:
                        continue
                    outf.write(f"{row}")

    outf.close()
    logger.info(f"save to {outfn}")


def report_raw_fast5_cpg_in_regions_table():
    """
    Generate Figure 2 C and D data: stats #cpgs of raw fast5 sequencing in each region
    :return:
    """
    dataset = []
    # distribution_dataset = defaultdict(list)
    for dsname in dsname_list:
        fnlist = glob.glob(os.path.join(args.base_cov_dir, f'{dsname}.rawfast5.coverage.base.bed.gz'))
        if len(fnlist) != 1:
            raise Exception(f"Not only one file, fnlist={fnlist}, for dsname={dsname}")

        fn = fnlist[0]
        logger.info(fn)

        rawReadBed = BedTool(fn).sort()
        logger.info(f'len(rawReadBed)={len(rawReadBed):,}')
        dataDict = {'dsname': dsname, 'filename': fn, 'total': len(rawReadBed)}

        # list of return for each parrallellized jobs
        retList = []

        ## we now use to_dataframe() for BedTool
        bed_df_cov_col = rawReadBed.to_dataframe().iloc[:, 3]

        ## total sites
        ret1 = {}
        for cutoff in cutoff_list:
            ret1.update({f'total_cutoff.cutoff{cutoff}': (bed_df_cov_col >= cutoff).sum()})

        logger.info(f"Evaluated regions: {eval_regions_file_path}")
        with Pool(processes=args.processors) as pool:
            for bedfn in eval_regions_file_path:
                tagname = get_region_tagname(bedfn)
                ret = pool.apply_async(count_sites_in_coord, (rawReadBed, bedfn, tagname,),
                                       kwds={'cutoff_list': cutoff_list})
                retList.append(ret)

            concordantFileName = find_bed_filename(basedir=bedDir,
                                                   pattern=f'*{dsname}*.concordant.bed')
            ret = pool.apply_async(count_sites_in_coord, (rawReadBed, concordantFileName, 'Concordant',),
                                   kwds={'cutoff_list': cutoff_list})
            retList.append(ret)

            discordantFileName = find_bed_filename(basedir=bedDir,
                                                   pattern=f'*{dsname}*.discordant.bed')
            ret = pool.apply_async(count_sites_in_coord, (rawReadBed, discordantFileName, 'Discordant',),
                                   kwds={'cutoff_list': cutoff_list})
            retList.append(ret)

            pool.close()
            pool.join()
        # Get each jobs's return results
        retList = [ret.get() for ret in retList]

        # Update returned results into dataDict: each line of results in the table
        dataDict.update(ret1)
        for ret in retList:
            if not ret:
                continue
            dataDict.update(ret)
        dataset.append(dataDict)
        logger.debug(dataDict)

        ## Temp output for real-time view
        save_dataset(dataset, tagname='tmp')

    save_dataset(dataset)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='computeRawReadsCoverage (NANOME)',
                                     description='compute coverage summary for raw reads in nanome paper')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument("cmd", help="name of command")
    parser.add_argument('--session-name', type=str, help='run name, default is "RawReadsCompute"',
                        default='RawReadsCompute')
    parser.add_argument('--genome-annotation', type=str,
                        help='genome annotation dir, contain BED files such as singleton, nonsingleton, etc.',
                        default=None)
    parser.add_argument('--reference-genome', type=str, help='reference genome file',
                        default=reference_genome_hg38_fn)
    parser.add_argument('--base-cov-dir', type=str, help='raw fast5 base coverage dir',
                        default=base_cov_dir)
    parser.add_argument('--bed-dir', type=str, help='bed dir for concordant and discordant',
                        default=bedDir)
    parser.add_argument('--processors', type=int, help="number of processors used, default is 8", default=8)
    parser.add_argument('--disable-bed-check',
                        help="if disable auto-checking the 0/1 base format for genome annotations",
                        action='store_true')
    parser.add_argument('--bedtools-tmp', type=str, help=f'bedtools temp dir, default is {global_temp_dir}',
                        default=global_temp_dir)
    parser.add_argument('-o', type=str, help=f"output base dir, default is {pic_base_dir}", default=pic_base_dir)
    parser.add_argument('--config', help="if print out config file for genome annotation", action='store_true')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()

    logger.debug(args)

    bed_temp_dir = os.path.join(args.bedtools_tmp, args.session_name)
    os.makedirs(bed_temp_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(bed_temp_dir)

    outdir = os.path.join(args.o, args.session_name)
    os.makedirs(outdir, exist_ok=True)
    logger.info(f"Output to dir: {outdir}")

    refFasta = None

    # Combine NA12878 all chrs bed regions
    if args.cmd == 'combine-na12878':
        combine_na12878_coverage_bed()

    # Generate bed sites file from region file using sequencing string
    if args.cmd == 'region-to-base':
        refFasta = get_ref_fasta(ref_fn=args.reference_genome)
        preprocess_bed_to_cgp_base()

    # Calculate raw fast5 cpg sites in each region
    if args.cmd == 'report-raw':
        dsname_list = ['HL60', 'K562', 'APL', 'NA19240', 'NA12878']
        # Nanopore reads coverage cutoff
        cutoff_list = [3]

        # list of evaluated regions (except for Genome-wide)
        # eval_regions = [os.path.join(args.genome_annotation, cofn) for cofn in narrowCoordNameList[1:]] + \
        #                [os.path.join(args.genome_annotation, cofn) for cofn in
        #                 cg_density_coord_name_list] + \
        #                [os.path.join(args.genome_annotation, cofn) for cofn in rep_coord_name_list]
        annot_dir = args.genome_annotation if args.genome_annotation is not None else '.'
        eval_regions_file_path = [os.path.join(annot_dir, cofn) for cofn in region_filename_dict.keys()]

        report_raw_fast5_cpg_in_regions_table()
    logger.info("### computeRawCoverage DONE")
