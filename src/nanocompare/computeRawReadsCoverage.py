#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : computeRawReadsCoverage.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

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

from nanocompare.eval_common import get_dna_seq_from_reference, open_file_gz_or_txt, find_bed_filename, get_ref_fasta, \
    get_region_bed, intersect_bed_regions
from nanocompare.global_config import set_log_debug_level, logger, pic_base_dir, temp_dir, data_base_dir
from nanocompare.global_settings import humanChrSet, location_filename_to_abbvname, \
    datasets_order, narrowCoordNameList, cg_density_coord_name_list, \
    rep_coord_name_list, referenceGenomeFile, nanome_version

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

    outfn = os.path.join(pic_base_dir, f'{dsname}.rawfast5.coverage.base.bed.gz')
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

            if chr not in humanChrSet:  # filter out non-human chrs
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

    coordBed = get_region_bed(coordfn)
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

        # outfn = os.path.join(pic_base_dir,
        #                      f'raw.fast5.reads.cpg.coverage.across.regions.cutoff{cutoff}.not.sorted{tagname}.xlsx')
        # outdf.to_excel(outfn)

        ## Reorder columns based on list, if not exist, skit it
        selcols = [colname for colname in column_list if colname in outdf.columns]
        outdf = outdf[selcols]

        outdf['dsname'] = pd.Categorical(df['dsname'], datasets_order)
        outdf = outdf.sort_values(by='dsname', ascending=True)
        outfn = os.path.join(pic_base_dir,
                             f'raw.fast5.reads.cpg.coverage.across.regions.cutoff{cutoff}{tagname}.table.s1.xlsx')
        outdf.to_excel(outfn)
        logger.info(f'save to {outfn}')


def combine_na12878_coverage_bed():
    baseDir = "/fastscratch/liuya/nanocompare/NA12878-coverage"
    outfn = os.path.join(pic_base_dir, "NA12878-allChrs.coverage.bothstrand.bed.gz")
    outf = gzip.open(outfn, 'wt')
    for chrName in humanChrSet:
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

        ret1 = {}
        for cutoff in cutoff_list:
            ret1.update({f'total_cutoff.cutoff{cutoff}': (bed_df_cov_col >= cutoff).sum()})

        logger.info(f"Evaluated regions: {eval_regions}")
        with Pool(processes=16) as pool:
            for bedfn in eval_regions:
                tagname = location_filename_to_abbvname[os.path.basename(bedfn)]
                ret = pool.apply_async(count_sites_in_coord, (rawReadBed, bedfn, tagname,),
                                       kwds={'cutoff_list': cutoff_list})
                retList.append(ret)

            concordantFileName = find_bed_filename(basedir=bedDir,
                                                   pattern=f'{dsname}*hg38_nonsingletons.concordant.bed')
            ret = pool.apply_async(count_sites_in_coord, (rawReadBed, concordantFileName, 'Concordant',),
                                   kwds={'cutoff_list': cutoff_list})
            retList.append(ret)

            discordantFileName = find_bed_filename(basedir=bedDir,
                                                   pattern=f'{dsname}*hg38_nonsingletons.discordant.bed')
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
            dataDict.update(ret)
        dataset.append(dataDict)
        logger.debug(dataDict)

        ## Temp output for real-time view
        save_dataset(dataset, tagname='tmp')

    save_dataset(dataset)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='comp_raw_read_cov (NANOME)',description='Plot and export data for Nanocompare paper.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{nanome_version}')
    parser.add_argument("cmd", help="name of command")
    parser.add_argument('--bedtools-tmp', type=str, help='bedtools temp dir', default=temp_dir)
    parser.add_argument('--genome-annotation', type=str, help='genome annotation dir',
                        default=os.path.join(data_base_dir, 'genome-annotation'))
    parser.add_argument('--reference-genome', type=str, help='reference genome file',
                        default=referenceGenomeFile)
    parser.add_argument('--base-cov-dir', type=str, help='raw fast5 base coverage dir',
                        default=base_cov_dir)
    parser.add_argument('--bed-dir', type=str, help='bed dir for concordant and discordant',
                        default=bedDir)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    set_log_debug_level()
    args = parse_arguments()
    logger.debug(args)

    os.makedirs(args.bedtools_tmp, exist_ok=True)
    pybedtools.helpers.set_tempdir(args.bedtools_tmp)

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
        # eval_regions = narrowCoordFileList[1:] + cg_density_file_list + rep_file_list

        eval_regions = [os.path.join(args.genome_annotation, cofn) for cofn in narrowCoordNameList[1:]] + \
                       [os.path.join(args.genome_annotation, cofn) for cofn in
                        cg_density_coord_name_list] + \
                       [os.path.join(args.genome_annotation, cofn) for cofn in rep_coord_name_list]

        report_raw_fast5_cpg_in_regions_table()
    logger.info("### computeRawCoverage DONE")
