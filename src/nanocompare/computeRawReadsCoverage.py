#!/usr/bin/env python3
"""
Calculate Nanopore raw reads coverage at each CpG sites based on raw fast5 files
"""
import glob
import os
import re
from multiprocessing import Pool

import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm

from nanocompare.eval_common import get_dna_seq_from_reference, open_file_gz_or_txt, find_bed_filename
from nanocompare.global_config import set_log_debug_level, logger, pic_base_dir
from nanocompare.global_settings import humanChrSet, narrowCoordFileList, location_filename_to_abbvname

# used for convert region bed cov to base level cov
rawReadDir = '/pod/2/li-lab/Nanopore_compare/data/Nanopore_cov'
# rawReadDir = '/projects/li-lab/yang/results/2021-04-05'

# base level cov bed files (combined + - strand)
baseReadCovDir = '/projects/li-lab/yang/results/2021-04-05'

# modify this dir to newly concordant and discordant perf results base dir
datasetBedDir = '/projects/li-lab/yang/results/2021-04-12'


def convert_region_to_cpg_base(dsname):
    """
    Conver bed file of regions into CpG base bed files

    Assume infn is 0-based file, we print 1-based output. Because 1-based site level reports are easy to stats overlap with region files.
    :param infn:
    :return:
    """

    fnlist = glob.glob(os.path.join(rawReadDir, f'{dsname}*coverage.positivestrand.bed'))
    fnlist += glob.glob(os.path.join(rawReadDir, f'{dsname}*coverage.negativestrand.bed'))

    logger.info(f'convert_region_to_cpg_base:{fnlist}')

    outfn = os.path.join(pic_base_dir, f'{dsname}.rawfast5.coverage.base.bed')
    outfile = open(outfn, 'w')

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
    # dsname_list = ['HL60', 'K562', 'APL']
    # fnlist = sorted(glob.glob(os.path.join(rawReadDir, '*.strand.bed')), key=os.path.getsize)
    # dsname_list = ['NA19240']

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
    covList = []
    for row in readBed:
        covList.append(int(row[3]))
    covSeries = pd.Series(covList)

    for i, cutoff in enumerate(cutoff_list):
        ret.update({f'total_sites_after_cutoff.cutoff{cutoff}': (covSeries >= cutoff).sum()})

    coordBed = BedTool(coordfn)
    coordBed = coordBed.sort()
    intersectBed = readBed.intersect(coordBed, u=True, wa=True)

    covList = []
    for row in intersectBed:
        cov = int(row[3])
        covList.append(cov)
    covSeries = pd.Series(covList)
    for i, cutoff in enumerate(cutoff_list):
        ret.update({f'{tagname}.cutoff{cutoff}': (covSeries >= cutoff).sum()})

    return ret


def report_table():
    """
    Generate Figure 2 C and D data
    :return:
    """
    dsname_list = ['HL60', 'K562', 'APL', 'NA19240']
    # dsname_list = ['HL60', 'K562', 'APL']
    # dsname_list = ['NA19240']
    # dsname_list = ['HL60']

    cutoff_list = [1, 3]

    dataset = []
    for dsname in dsname_list:
        fnlist = glob.glob(os.path.join(baseReadCovDir, f'{dsname}*.base.bed'))
        fn = fnlist[0]
        logger.info(fn)

        rawReadBed = BedTool(fn)
        rawReadBed = rawReadBed.sort()
        logger.info(f'len(rawReadBed)={len(rawReadBed):,}')
        dataDict = {'dsname': dsname, 'filename': fn, 'total': len(rawReadBed)}
        retList = []
        with Pool(processes=8) as pool:
            for bedfn in narrowCoordFileList[1:]:
                tagname = location_filename_to_abbvname[os.path.basename(bedfn)]
                ret = pool.apply_async(count_sites_in_coord, (rawReadBed, bedfn, tagname,), kwds={'cutoff_list': cutoff_list})
                retList.append(ret)

            concordantFileName = find_bed_filename(basedir=datasetBedDir, pattern=f'{dsname}*hg38_nonsingletons.concordant.bed')
            ret = pool.apply_async(count_sites_in_coord, (rawReadBed, concordantFileName, 'Concordant',), kwds={'cutoff_list': cutoff_list})
            retList.append(ret)

            discordantFileName = find_bed_filename(basedir=datasetBedDir, pattern=f'{dsname}*hg38_nonsingletons.discordant.bed')
            ret = pool.apply_async(count_sites_in_coord, (rawReadBed, discordantFileName, 'Discordant',), kwds={'cutoff_list': cutoff_list})
            retList.append(ret)

            pool.close()
            pool.join()
        retList = [ret.get() for ret in retList]
        for ret in retList:
            dataDict.update(ret)

        dataset.append(dataDict)
        logger.debug(dataDict)
    df = pd.DataFrame(dataset)
    logger.info(df)

    for cutoff in cutoff_list:
        outdf = pd.concat([df.loc[:, ['dsname', 'filename', 'total']], df.filter(regex=f'.cutoff{cutoff}$', axis=1)], axis=1)
        outfn = os.path.join(pic_base_dir, f'raw.fast5.reads.cpg.coverage.across.regions.cutoff{cutoff}.xlsx')
        outdf.to_excel(outfn)
        logger.info(f'save to {outfn}')


if __name__ == '__main__':
    set_log_debug_level()

    refFasta = None
    # refFasta = get_ref_fasta()
    # preprocess_bed_to_cgp_base()
    report_table()
